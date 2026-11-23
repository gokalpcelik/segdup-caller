"""CNV/depth calling: CNVscope wrapper and depth/segment caching.

`call_cnvscope`, `call_depth`, and `trim_to_region_vcf` are pure helpers
taking a Bam instance as a state-holder. `call_depth` mutates Bam's depth/
segment caches as part of its caching contract. `source_cnvscope` runs the
source CNVscope of every chromosome once, before the genes start; `call_depth`
then reads those outputs by path (`staged_run`).
"""

import os
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import vcflib

from .. import throttle, tools
from ..util import IntervalList

DEFAULT_BIN_SIZE = 1000
CONFIDENT_MAPQ = 30  # reads at or above this MAPQ carry "confident" depth
MIN_CONFIDENT_FRAC = 0.95  # bins used to set the diploid depth must be this mappable
NORM_BINS_PER_CHROM = 100
MIN_AUTOSOME_LEN = 10_000_000  # autosomes shorter than this are not used to set the diploid depth


def _inputs(bam: Any) -> List[str]:
    # a lifted Bam may read several BAMs at once (cnv_inputs); an input Bam reads bam.bam
    return list(getattr(bam, "cnv_inputs", None) or [bam.bam])


def _caller(bam: Any) -> str:
    return getattr(bam, "cnv_caller", "depth") or "depth"


def cnvscope_template(bam: Any, param: Dict[str, Any]) -> str:
    """Key identifying a depth run: two runs with the same key produce the same output."""
    return "|".join(
        [
            "native_depth" if _caller(bam) == "depth" else f"{_caller(bam)}_depth",
            ",".join(_inputs(bam)),
            str(param.get("region", "")),
            str(getattr(bam, "cnv_bin_size", DEFAULT_BIN_SIZE)),
        ]
    )


def _run(bam: Any, cmd: str) -> str:
    bam.logger.debug(f"Running command: {cmd}")
    with throttle.slot():
        result = subprocess.run(
            cmd, capture_output=True, text=True, shell=True, executable="/bin/bash"
        )
    if result.returncode:
        raise Exception(result.stderr)
    return result.stdout


def _fai_lengths(ref: str) -> Dict[str, int]:
    with open(ref + ".fai") as fh:
        return {f[0]: int(f[1]) for f in (line.split("\t") for line in fh)}


def _bins(region: str, lengths: Dict[str, int], size: int) -> List[Tuple[str, int, int]]:
    """Fixed-grid bins (0-based half-open) covering ``region`` ("chr" or "chr:s-e", comma separated)."""
    out = []
    for part in [p for p in region.replace(" ", ",").split(",") if p]:
        if ":" in part:
            chrom, span = part.rsplit(":", 1)
            s, e = (int(x) for x in span.split("-"))
            s -= 1
        else:
            chrom, s, e = part, 0, lengths[part]
        e = min(e, lengths.get(chrom, e))
        pos = (s // size) * size
        while pos < e:
            out.append((chrom, max(pos, s), min(pos + size, e)))
            pos += size
    return out


def _depth_table(
    bam: Any, bins: List[Tuple[str, int, int]], bed: str, inputs: List[str]
) -> Tuple[np.ndarray, np.ndarray]:
    """Mean depth per bin: (all reads, reads with MAPQ >= CONFIDENT_MAPQ), summed over ``inputs``."""
    with open(bed, "w") as fh:
        for c, s, e in bins:
            fh.write(f"{c}\t{s}\t{e}\n")

    def one(q: int) -> np.ndarray:
        # one bedcov per input so each is queried by its own index; the sums are added
        total = np.zeros(len(bins))
        for p in inputs:
            out = _run(bam, tools.bedcov_cmd(bam.samtools, bed, p, bam.ref, q))
            total += np.array([r[3] for r in tools.parse_bedcov(out)])
        return total

    with ThreadPoolExecutor(max_workers=2) as pool:
        f_all, f_conf = pool.submit(one, 0), pool.submit(one, CONFIDENT_MAPQ)
        return f_all.result(), f_conf.result()


_DIPLOID_CACHE: Dict[Tuple[str, str], float] = {}


def diploid_depth(bam: Any) -> float:
    """Typical per-base depth of a diploid unique region of this sample.

    Median confident depth over bins spread across the middle of the autosomes, taken from the
    sample's original BAM (also for a lifted Bam, whose merge holds the gene regions only).
    Cached per (input, reference).
    """
    src = _inputs(bam)[0]
    key = (src, bam.ref)
    if key in _DIPLOID_CACHE:
        return _DIPLOID_CACHE[key]
    lengths = _fai_lengths(bam.ref)
    size = getattr(bam, "cnv_bin_size", DEFAULT_BIN_SIZE)
    bins = []
    for chrom, ln in lengths.items():
        name = chrom[3:] if chrom.startswith("chr") else chrom
        if not name.isdigit() or not 1 <= int(name) <= 22 or ln < MIN_AUTOSOME_LEN:
            continue
        lo, hi = int(ln * 0.1), int(ln * 0.9)
        step = max(size, (hi - lo) // NORM_BINS_PER_CHROM)
        bins += [(chrom, p, p + size) for p in range(lo, hi - size, step)]
    if not bins:
        raise Exception(
            "Cannot estimate the diploid depth: no usable autosome in the reference index"
        )
    fd, bed = tempfile.mkstemp(suffix=".bed", dir=bam.tmpdir)
    os.close(fd)
    try:
        d_all, d_conf = _depth_table(bam, bins, bed, [src])
    finally:
        os.remove(bed)
    ok = (d_all > 0) & (d_conf >= MIN_CONFIDENT_FRAC * d_all)
    if ok.sum() < 20:
        raise Exception(
            f"Too few well-mapped bins ({int(ok.sum())}) in {src} to estimate the diploid depth"
        )
    value = float(np.median(d_conf[ok]))
    if value <= 0:
        raise Exception(f"Estimated diploid depth of {src} is zero")
    _DIPLOID_CACHE[key] = value
    return value


def _segments(rows: pd.DataFrame, cn: np.ndarray) -> List[Tuple[int, int]]:
    """Index ranges [i, j) of adjacent bins sharing a copy number and a mapping-quality class."""
    out, i, n = [], 0, len(rows)
    contig, start, end = rows["contig"].values, rows["start"].values, rows["end"].values
    hq = rows["MQ"].values > 50
    while i < n:
        j = i + 1
        while (
            j < n
            and contig[j] == contig[i]
            and start[j] == end[j - 1]
            and cn[j] == cn[i]
            and hq[j] == hq[i]
        ):
            j += 1
        out.append((i, j))
        i = j
    return out


def _write_vcf(bam: Any, rows: pd.DataFrame, cn: np.ndarray, out_vcf: str) -> None:
    lengths = _fai_lengths(bam.ref)
    chroms = list(dict.fromkeys(rows["contig"]))
    lines = ["##fileformat=VCFv4.2"]
    lines += [f"##contig=<ID={c},length={lengths[c]}>" for c in chroms]
    lines += [
        '##INFO=<ID=END,Number=1,Type=Integer,Description="Segment end">',
        '##INFO=<ID=CNT,Number=1,Type=Integer,Description="Copy number (2 x normalized depth, rounded)">',
        '##FORMAT=<ID=MQS,Number=1,Type=Float,Description="Mean MAPQ proxy of the segment">',
        '##FORMAT=<ID=DPS,Number=2,Type=Float,Description="Mean and SD of normalized depth">',
        '##FORMAT=<ID=DP0S,Number=2,Type=Float,Description="Mean and SD of normalized depth, all reads">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
    ]
    for i, j in _segments(rows, cn):
        r = rows.iloc[i:j]
        lines.append(
            "\t".join(
                [
                    str(r["contig"].iloc[0]),
                    str(int(r["start"].iloc[0]) + 1),
                    ".",
                    "N",
                    "<CNV>",
                    ".",
                    ".",
                    f"END={int(r['end'].iloc[-1])};CNT={int(cn[i])}",
                    "MQS:DPS:DP0S",
                    f"{r['MQ'].mean():.2f}:{r['DP'].mean():.4f},{r['DP'].std(ddof=0):.4f}"
                    f":{r['DP0'].mean():.4f},{r['DP0'].std(ddof=0):.4f}",
                ]
            )
        )
    raw = out_vcf + ".tmp.vcf"
    with open(raw, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    try:
        _run(bam, tools.compress_index_vcf_cmd(bam.bcftools, out_vcf, raw))
    finally:
        os.remove(raw)


CNVKIT_NULL_LOG2 = -20.0  # CNVkit's log2 for a bin with no coverage
CNVKIT_SEGMENT_METHOD = "haar"  # pure Python; cbs needs R, and hmm failed on integer weights in 0.9.14
_CNN_COLS = ["chromosome", "start", "end", "gene", "depth", "log2"]


def _cnvkit(bam: Any, args: str) -> str:
    return _run(bam, f"{tools._q(getattr(bam, 'cnvkit', 'cnvkit.py'))} {args}")


def _cnvkit_correct(bam: Any, rows: pd.DataFrame, dp_raw: np.ndarray, output: str) -> np.ndarray:
    """Replace ``rows["DP"]`` by CNVkit's bias-corrected, segmented depth; return per-bin copy number.

    The grid's confident depth (``dp_raw``, diploid = 1) goes through ``cnvkit fix`` against a flat
    reference (GC, edge and repeat-mask bias correction) and ``cnvkit segment``. CNVkit median-centres
    each run, which would hide a deletion or gain covering most of a short run, so the corrected
    depth is re-anchored to the sample's global diploid depth with the median per-bin ratio.
    Bins CNVkit drops (no coverage) keep their raw depth.
    """
    prefix = output[: output.index(".vcf")] if ".vcf" in output else output
    wd = tempfile.mkdtemp(prefix="cnvkit.", dir=bam.tmpdir)
    try:
        contig = rows["contig"].values
        start, end = rows["start"].values, rows["end"].values
        depth = dp_raw  # relative depth; CNVkit only needs the ratios
        pos = depth[depth > 0]
        ref_level = float(np.median(pos)) if len(pos) else 1.0
        log2 = np.where(depth > 0, np.log2(np.where(depth > 0, depth, 1.0) / ref_level), CNVKIT_NULL_LOG2)
        cnn = pd.DataFrame(
            {"chromosome": contig, "start": start, "end": end, "gene": "-", "depth": depth, "log2": log2}
        )
        t_cnn, a_cnn, r_cnn, cnr, cns, grid = (
            os.path.join(wd, f) for f in ("t.cnn", "a.cnn", "flat.cnn", "s.cnr", "s.cns", "grid.bed")
        )
        cnn.to_csv(t_cnn, sep="\t", index=False)
        with open(a_cnn, "w") as fh:
            fh.write("\t".join(_CNN_COLS) + "\n")
        cnn[["chromosome", "start", "end", "gene"]].to_csv(grid, sep="\t", index=False, header=False)
        _cnvkit(bam, f"reference -f {tools._q(bam.ref)} -t {tools._q(grid)} -o {tools._q(r_cnn)}")
        _cnvkit(bam, f"fix {tools._q(t_cnn)} {tools._q(a_cnn)} {tools._q(r_cnn)} -o {tools._q(cnr)}")
        _cnvkit(bam, f"segment {tools._q(cnr)} -m {CNVKIT_SEGMENT_METHOD} -o {tools._q(cns)}")
        fixed = pd.read_csv(cnr, sep="\t")
        segs = pd.read_csv(cns, sep="\t")
    finally:
        shutil.rmtree(wd, ignore_errors=True)

    # join the corrected bins back onto the grid (dropped bins stay NaN)
    key = pd.MultiIndex.from_arrays([contig, start])
    corr = pd.Series(
        np.power(2.0, fixed["log2"].values),
        index=pd.MultiIndex.from_arrays([fixed["chromosome"].astype(str).values, fixed["start"].values]),
    )
    corr = corr[~corr.index.duplicated()].reindex(key).values
    seg_log2 = np.full(len(rows), np.nan)
    order = {c: i for i, c in enumerate(dict.fromkeys(contig))}
    for c, g in segs.groupby(segs["chromosome"].astype(str)):
        if c not in order:
            continue
        idx = np.flatnonzero(contig == c)
        s_arr, e_arr = start[idx], end[idx]
        for s, e, v in zip(g["start"].values, g["end"].values, g["log2"].values):
            lo, hi = np.searchsorted(s_arr, s, "left"), np.searchsorted(e_arr, e, "right")
            seg_log2[idx[lo:hi]] = v
    ok = np.isfinite(corr) & (corr > 0) & (dp_raw > 0)
    factor = float(np.median(dp_raw[ok] / corr[ok])) if ok.sum() else 1.0
    new_dp = np.where(np.isfinite(corr), corr * factor, dp_raw)
    seg_dp = np.where(np.isfinite(seg_log2), np.power(2.0, seg_log2) * factor, new_dp)
    seg_dp = np.where(np.isfinite(corr), seg_dp, dp_raw)
    ratio = np.divide(new_dp, dp_raw, out=np.ones_like(new_dp), where=dp_raw > 0)
    rows["DP"] = new_dp
    rows["DP0"] = rows["DP0"].values * ratio  # keep the all-read depth on the corrected scale
    return np.rint(2 * seg_dp).astype(int)


def call_cnvscope(bam: Any, output: str, param: Dict[str, Any]) -> None:
    """Depth-based copy-number segmentation over ``param["region"]`` (replaces Sentieon CNVscope).

    With ``bam.cnv_caller == "cnvkit"`` the depth is GC/edge/repeat-corrected and segmented by
    CNVkit (`_cnvkit_correct`) instead of the plain rounding below; the outputs are identical in form.

    Writes the per-bin table ``param["output_norm"]`` (contig,start,end,MQ,DP,DP0) and a segment
    VCF ``output`` (INFO CNT; FORMAT MQS, DPS, DP0S) that ``call_depth`` reads. Depth is normalized
    so a diploid unique region is 1.0. MQ is a proxy: 60 x the fraction of depth carried by reads
    with MAPQ >= CONFIDENT_MAPQ. The default caller applies no GC or mappability correction.
    """
    if bam.use_existing and os.path.exists(output) and os.path.getsize(output):
        return
    if "region" not in param:
        raise Exception("A region is required for the depth caller.")
    norm_path = param["output_norm"]
    size = getattr(bam, "cnv_bin_size", DEFAULT_BIN_SIZE)
    bins = _bins(param["region"], _fai_lengths(bam.ref), size)
    if not bins:
        raise Exception(f"No bins for region {param['region']}")
    scale = diploid_depth(bam)
    fd, bed = tempfile.mkstemp(suffix=".bed", dir=bam.tmpdir)
    os.close(fd)
    try:
        d_all, d_conf = _depth_table(bam, bins, bed, _inputs(bam))
    finally:
        os.remove(bed)
    frac = np.where(d_all > 0, d_conf / np.where(d_all > 0, d_all, 1.0), 0.0)
    rows = pd.DataFrame(
        {
            "contig": [b[0] for b in bins],
            "start": [b[1] for b in bins],
            "end": [b[2] for b in bins],
            "MQ": 60.0 * np.minimum(frac, 1.0),
            "DP": d_conf / scale,
            "DP0": d_all / scale,
        }
    )
    if _caller(bam) == "cnvkit":
        cn = _cnvkit_correct(bam, rows, d_conf / scale, output)
    else:
        cn = np.rint(2 * rows["DP"].values).astype(int)
    rows.to_csv(norm_path, sep="\t", index=False)
    _write_vcf(bam, rows, cn, output)


def _segment_param(bam: Any, chrom: str, output_norm: str = "") -> Dict[str, Any]:
    # what a gene's segmentation passes call_cnvscope for one chromosome
    return {"region": chrom, "output_norm": output_norm}


def source_cnvscope(
    bams: Dict[str, Any], chroms: List[str], out_dir: str, max_parallel: int
) -> Dict[Tuple[str, str], Tuple[str, str]]:
    """Each input's source CNVscope per chromosome: {(input tag, chromosome): (VCF, TSV)}.

    The one place that decides how the stage runs CNVscope; callers see only these outputs. It runs
    the genes' own command per input and chromosome, max_parallel at a time: one call over several
    chromosomes normalises them together and changes the depths. A failed run is left out, and
    its genes run their own.
    """

    def run(job):
        tag, bam, chrom = job
        vcf, norm = (os.path.join(out_dir, f"{tag}.{chrom}{x}") for x in (".vcf.gz", ".norm.tsv"))
        try:
            call_cnvscope(bam, vcf, _segment_param(bam, chrom, norm))
        except Exception as e:
            bam.logger.error(f"Source CNVscope on {chrom} failed; its genes run their own: {e}")
            return None
        return (tag, chrom), (vcf, norm)

    jobs = [(tag, bam, chrom) for tag, bam in bams.items() for chrom in chroms]
    with ThreadPoolExecutor(max_workers=max(1, min(max_parallel, len(jobs)))) as pool:
        return dict(r for r in pool.map(run, jobs) if r)


def staged_table(
    bams: Dict[str, Any], outputs: Dict[Tuple[str, str], Tuple[str, str]]
) -> Dict[str, Tuple[str, str]]:
    """{command template: (VCF, TSV)}: each staged output under the gene command it stands for."""
    return {
        cnvscope_template(bams[tag], _segment_param(bams[tag], chrom)): paths
        for (tag, chrom), paths in outputs.items()
    }


def flank_factor(
    bam: Any,
    flank_regions: List[str],
    tag: str = "short_read",
    min_bins: int = 50,
    bounds: Tuple[float, float] = (0.8, 1.25),
) -> Tuple[float, int]:
    """Local depth baseline from single-copy flanks: (median depth of the flank bins, bin count).

    The depths are normalised to the genome-wide diploid depth, but a locus can sit above or below
    it (GC, mappability, library bias), so its copy numbers read high or low. The median MAPQ > 50
    depth of unique flanking bins around the gene is that locus's own diploid level; dividing the
    gene's depths by it rescales them to the local baseline. Returns (1.0, n) when the staged depth
    table is missing, fewer than min_bins flank bins are found, or the factor falls outside bounds.
    """
    table = getattr(bam, "source_cnv", None) or {}
    if isinstance(flank_regions, str):  # accept "a:1-2,b:3-4" as well as a list
        flank_regions = [r.strip() for r in flank_regions.split(",") if r.strip()]
    itv = IntervalList(region=",".join(flank_regions))
    chroms = set(itv.regions)
    values: List[float] = []
    for _, tsv in table.values():
        base = os.path.basename(tsv)
        chrom = next((c for c in chroms if base == f"{tag}.{c}.norm.tsv"), None)
        if chrom is None:
            continue
        df = pd.read_csv(tsv, sep="\t")
        df["contig"] = df["contig"].astype(str)
        df = df[(df["contig"] == chrom) & (df["MQ"] > 50)]
        mid = (df["start"] + df["end"]) // 2
        for region in flank_regions:
            c, rng = region.split(":")
            if c != chrom:
                continue
            a, b = (int(x) for x in rng.replace(",", "").split("-"))
            values.extend(df.loc[(mid >= a) & (mid <= b), "DP"].tolist())
    n = len(values)
    if n < min_bins:
        return 1.0, n
    factor = float(np.median(values))
    if not bounds[0] <= factor <= bounds[1]:
        return 1.0, n
    return factor, n


def staged_run(bam: Any, param: Dict[str, Any]) -> Optional[Tuple[str, str]]:
    """(VCF, TSV) the up-front stage wrote for this exact command, or None to run it here."""
    table = getattr(bam, "source_cnv", None)
    if table is None:  # no stage for this Bam (a lifted Bam, or a run without one)
        return None
    staged = table.get(cnvscope_template(bam, param))
    bam.logger.info(
        f"Source CNVscope for {param.get('region')}: "
        + (f"staged output {staged[0]}" if staged else "no staged output, running its own")
    )
    return staged


def _reference_order(chroms: Any, ref: str) -> List[str]:
    """chroms in the reference's .fai order; a set's own order follows the hash seed."""
    if len(chroms) < 2:
        return list(chroms)
    try:
        with open(ref + ".fai") as fh:
            rank = {line.split("\t", 1)[0]: i for i, line in enumerate(fh)}
    except OSError:
        rank = {}
    return sorted(chroms, key=lambda c: (rank.get(c, len(rank)), c))


def call_depth(bam: Any, output: str, param: Dict[str, Any]) -> None:
    # Import locally to avoid a circular import with bam_process at module load.
    from ..bam_process import Segment, DepthSegment

    orig_region = None
    if "region" in param:
        orig_region = param["region"]
        chrs = set()
        for parts in param["region"].split(","):
            chrs.add(parts.split(":")[0])
        param["region"] = ",".join(_reference_order(chrs, bam.ref))
    prefix = output[: output.index(".vcf")] if ".vcf" in output else output
    output_norm = param.get("output_norm", f"{prefix}.norm.tsv")
    if orig_region:
        tmp_output = prefix + ".tmp" + output[output.index(".vcf") :]
        tmp_output_norm = prefix + ".tmp.norm.tsv"
    else:
        tmp_output = output
        tmp_output_norm = output_norm
    param["output_norm"] = tmp_output_norm
    staged = staged_run(bam, param)
    if staged:  # the up-front stage ran this exact command
        tmp_output, tmp_output_norm = staged
    else:
        call_cnvscope(bam, tmp_output, param)
    if orig_region:
        param["region"] = orig_region
    vcf = vcflib.VCF(tmp_output)
    depth_fields = ["contig", "start", "end", "MQ", "DP", "DP0"]
    for v in vcf:
        dp = v.samples[0]["DPS"]
        dp0 = v.samples[0]["DP0S"]
        bam.segments.setdefault(v.chrom, []).append(
            Segment(
                v.chrom,
                v.pos,
                v.end,
                v.info["CNT"],
                v.samples[0]["MQS"],
                dp[0] / bam.dp_norm,
                dp[1] / bam.dp_norm,
                dp0[0] / bam.dp_norm,
                dp0[1] / bam.dp_norm,
            )
        )
    depth_df = pd.read_csv(tmp_output_norm, sep="\t")[depth_fields].copy()
    depth_df["contig"] = depth_df["contig"].astype(str)
    depth_df["DP"] /= bam.dp_norm
    depth_df["DP0"] /= bam.dp_norm
    high_mq_depth_df = depth_df[depth_df["MQ"] > 50]
    bam.overall_depth_std = high_mq_depth_df["DP"].std()
    # Optimize DepthSegment creation using vectorized operations
    bam.depths = {}
    bam._depth_coords_cache = {}
    for chr, group in depth_df.groupby("contig"):
        # Use vectorized numpy operations instead of itertuples
        values = group.values
        chr_depths = [DepthSegment(*row) for row in values]
        bam.depths[chr] = chr_depths
        # Cache coordinates directly from numpy arrays for better performance
        bam._depth_coords_cache[chr] = (
            values[:, 1].tolist(),  # start column
            values[:, 2].tolist(),  # end column
        )
    if orig_region:
        orig_region_itv = IntervalList(region=orig_region)

        # Vectorized filtering using the same logic as the earlier optimization
        mask = pd.Series(False, index=depth_df.index)

        for chr in orig_region_itv.regions:
            chr_mask = depth_df["contig"] == chr
            if not chr_mask.any():
                continue

            chr_regions = orig_region_itv.regions[chr]
            chr_rows = depth_df[chr_mask]

            i = 0
            while i < len(chr_regions):
                start, end = chr_regions[i], chr_regions[i + 1]

                # Vectorized overlap check
                overlap_mask = (chr_rows["start"] < end) & (chr_rows["end"] > start)
                mask[chr_rows[overlap_mask].index] = True

                i += 2

        trim_depth_df = depth_df[mask]
        trim_depth_df.to_csv(output_norm, sep="\t", index=None)
        trim_to_region_vcf(bam, orig_region, tmp_output, output)


def trim_to_region_vcf(bam: Any, region: str, in_vcf: str, out_vcf: str) -> None:
    if region.endswith(".bed"):
        region_opt = f"-T {region}"
    else:
        region_opt = f"-r {region}"
    cmd = (
        f"{bam.bcftools} view {region_opt} {in_vcf} "
        f"| {bam.bcftools} view -Oz -o {out_vcf} - && {bam.bcftools} index -f -t {out_vcf}"
    )
    result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
    result.check_returncode()

"""Command-line builders for the open-source tools the pipeline shells out to.

Every function returns a shell command string and does no I/O, so the commands
can be unit-tested without the tools installed. Callers run them with
``subprocess.run(cmd, shell=True)`` (under ``throttle.slot()`` for the heavy
ones). Paths are shell-quoted here; do not quote them again at the call site.

Tool map (replaces the former Sentieon commands):

* region extraction / merge  -> samtools view / merge / sort / index
* read re-alignment          -> minimap2 | samtools sort
* small-variant calling      -> gatk HaplotypeCaller
* VCF compress + index       -> bcftools view -Oz / bcftools index -t
* interval depth             -> samtools bedcov
"""

import os
import shlex
from typing import Iterable, List, Optional

DEFAULT_JAVA_OPTIONS = "-Xmx4g -XX:ParallelGCThreads=1"
DEFAULT_HC_MIN_MAPQ = 10

_q = shlex.quote


def _split_regions(regions: str) -> List[str]:
    """Split a "chr:s-e,chr:s-e" (or space separated) region string."""
    return [r for r in regions.replace(",", " ").split() if r]


def _alignment_format(path: str) -> str:
    return "cram" if path.endswith(".cram") else "bam"


def write_region_bed(regions: str, bed: str, padding: int = 0) -> str:
    """Write a merged, padded BED for ``regions`` and return its path."""
    from .util import IntervalList

    il = IntervalList(regions)
    if padding:
        il.pad(padding)
    with open(bed, "w") as fh:
        for chrom, spans in il.regions.items():
            for i in range(0, len(spans), 2):
                fh.write(f"{chrom}\t{spans[i]}\t{spans[i + 1]}\n")
    return bed


def extract_regions(
    samtools: str,
    ref: str,
    in_bam: str,
    bed: str,
    out: str,
    tmpdir: str,
    compression: Optional[int] = None,
) -> str:
    """Reads of ``in_bam`` overlapping ``bed`` -> coordinate-sorted, indexed ``out``.

    The multi-region iterator (-M) avoids emitting a read twice when regions
    overlap; the sort guarantees order across contigs.
    """
    fmt = _alignment_format(out)
    lvl = f" -l {compression}" if compression is not None else ""
    prefix = os.path.join(tmpdir, os.path.basename(out) + ".sorttmp")
    return (
        f"set -o pipefail; "
        f"{_q(samtools)} view -u -h -M -L {_q(bed)} --reference {_q(ref)} {_q(in_bam)} | "
        f"{_q(samtools)} sort{lvl} -O {fmt} --reference {_q(ref)} -T {_q(prefix)} -o {_q(out)} - && "
        f"{_q(samtools)} index {_q(out)}"
    )


def merge_bams_cmd(
    samtools: str,
    ref: str,
    in_bams: Iterable[str],
    out: str,
    compression: Optional[int] = None,
) -> str:
    """Merge coordinate-sorted inputs (same header lineage) into indexed ``out``."""
    fmt = _alignment_format(out)
    lvl = f" -l {compression}" if compression is not None else ""
    return (
        f"{_q(samtools)} merge -f{lvl} -O {fmt} --reference {_q(ref)} {_q(out)} "
        f"{' '.join(_q(b) for b in in_bams)} && {_q(samtools)} index {_q(out)}"
    )


def realign_cmd(
    samtools: str,
    minimap2: str,
    ref: str,
    in_bam: str,
    regions: str,
    target_fa: str,
    read_group: str,
    preset: str,
    out_bam: str,
    exclude_dups: bool = True,
    threads: int = 1,
) -> str:
    """Re-map reads from ``in_bam`` (restricted to ``regions``) onto ``target_fa``.

    ``read_group`` is the literal "@RG\\tID:..." header line. Dup/secondary/
    supplementary reads are dropped because ``samtools fastq`` discards FLAG bits.
    ``minimap2 -y`` carries the FASTQ comment (the ``-t`` tags) into the SAM.
    """
    dup = "-F 0x400 " if exclude_dups else ""
    return (
        f"set -o pipefail; "
        f"{_q(samtools)} view -F 0x100 -F 0x200 -F 0x800 {dup}-h --reference {_q(ref)} "
        f"{_q(in_bam)} {regions} | "
        f"{_q(samtools)} fastq -t | "
        f'{_q(minimap2)} -y -a -t {threads} -x {preset} -R "{read_group}" {_q(target_fa)} /dev/stdin | '
        f"{_q(samtools)} sort -O bam -o {_q(out_bam)} -"
    )


def haplotypecaller_cmd(
    gatk: str,
    ref: str,
    bam: str,
    output: str,
    regions: Optional[str] = None,
    ploidy: int = 2,
    given: Optional[str] = None,
    dbsnp: Optional[str] = None,
    tmpdir: Optional[str] = None,
    java_options: str = DEFAULT_JAVA_OPTIONS,
    min_mapq: int = DEFAULT_HC_MIN_MAPQ,
    extra: str = "",
) -> str:
    """GATK HaplotypeCaller. ``output`` must end in .vcf.gz (GATK then writes the .tbi).

    With ``given`` the call is genotype-given-alleles: only those alleles are emitted.
    """
    parts = [
        _q(gatk),
        f'--java-options "{java_options}"',
        "HaplotypeCaller",
        f"-R {_q(ref)}",
        f"-I {_q(bam)}",
        f"-O {_q(output)}",
        f"--sample-ploidy {int(ploidy)}",
        f"--minimum-mapping-quality {int(min_mapq)}",
        "--native-pair-hmm-threads 1",
        "--create-output-variant-index true",
    ]
    for r in _split_regions(regions or ""):
        parts.append(f"-L {_q(r)}")
    if given:
        parts.append(f"--alleles {_q(given)} --genotyping-mode GENOTYPE_GIVEN_ALLELES")
    if dbsnp:
        parts.append(f"--dbsnp {_q(dbsnp)}")
    if tmpdir:
        parts.append(f"--tmp-dir {_q(tmpdir)}")
    if extra:
        parts.append(extra)
    return " ".join(parts)


def compress_index_vcf_cmd(bcftools: str, out_vcf: str, src: str = "-") -> str:
    """``bcftools view`` of ``src`` (default stdin) -> bgzipped ``out_vcf`` + .tbi.

    Meant as the tail of a pipeline, e.g. ``f"bcftools view ... | {compress_index_vcf_cmd(...)}"``
    when ``src`` is "-", the leading pipe is the caller's.
    """
    return (
        f"{_q(bcftools)} view -Oz -o {_q(out_vcf)} {_q(src)} "
        f"&& {_q(bcftools)} index -f -t {_q(out_vcf)}"
    )


def index_vcf_cmd(bcftools: str, vcf: str) -> str:
    return f"{_q(bcftools)} index -f -t {_q(vcf)}"


def bedcov_cmd(samtools: str, bed: str, bam: str, ref: Optional[str] = None, min_mapq: int = 0) -> str:
    """``samtools bedcov`` -> "chrom start end sum_of_per_base_depth" lines on stdout."""
    r = f"--reference {_q(ref)} " if ref else ""
    return f"{_q(samtools)} bedcov -Q {int(min_mapq)} {r}{_q(bed)} {_q(bam)}"


def parse_bedcov(text: str) -> List[tuple]:
    """Parse bedcov output into (chrom, start, end, mean_depth) tuples."""
    out = []
    for line in text.splitlines():
        f = line.rstrip("\n").split("\t")
        if len(f) < 4:
            continue
        chrom, s, e, total = f[0], int(f[1]), int(f[2]), float(f[-1])
        length = e - s
        out.append((chrom, s, e, total / length if length > 0 else 0.0))
    return out

import random
import shutil
import subprocess

import pysam
import pytest

from genecaller.bam_process import Bam

need_tools = pytest.mark.skipif(
    not all(shutil.which(t) for t in ("samtools", "bcftools")),
    reason="samtools and bcftools are required",
)

CHROM_LEN = 200_000
READ_LEN = 100
BASE_DEPTH = 20


def _write_reads(bam_path, ref, depth_fn, mapq_fn):
    """Tile 100 bp reads over chr1/chr2 so mean depth follows depth_fn(chrom, pos) (a float)."""
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": c, "LN": CHROM_LEN} for c in ("chr1", "chr2")],
        "RG": [{"ID": "rg1", "SM": "S1"}],
    }
    rng = random.Random(1)
    with pysam.AlignmentFile(bam_path, "wb", header=header) as out:
        for tid, chrom in enumerate(("chr1", "chr2")):
            seq = ref[chrom]
            n = 0
            for pos in range(0, CHROM_LEN - READ_LEN, 5):
                # reads start every 5 bp -> 20x at weight 1.0
                w = depth_fn(chrom, pos) / BASE_DEPTH
                if rng.random() >= w:
                    continue
                a = pysam.AlignedSegment(out.header)
                a.query_name = f"{chrom}_{n}"
                n += 1
                a.query_sequence = seq[pos : pos + READ_LEN]
                a.flag = 0
                a.reference_id = tid
                a.reference_start = pos
                a.mapping_quality = mapq_fn(chrom, pos)
                a.cigartuples = [(0, READ_LEN)]
                a.query_qualities = pysam.qualitystring_to_array("I" * READ_LEN)
                a.set_tag("RG", "rg1")
                out.write(a)
    pysam.index(bam_path)


@pytest.fixture(scope="session")
def dataset(tmp_path_factory):
    d = tmp_path_factory.mktemp("data")
    rng = random.Random(7)
    ref = {c: "".join(rng.choice("ACGT") for _ in range(CHROM_LEN)) for c in ("chr1", "chr2")}
    fa = d / "ref.fa"
    with open(fa, "w") as fh:
        for c, s in ref.items():
            fh.write(f">{c}\n{s}\n")
    subprocess.run(["samtools", "faidx", str(fa)], check=True)

    def depth(chrom, pos):
        if chrom == "chr1" and 80_000 <= pos < 100_000:
            return BASE_DEPTH / 2  # heterozygous deletion
        return BASE_DEPTH

    def mapq(chrom, pos):
        return 0 if (chrom == "chr2" and 40_000 <= pos < 60_000) else 60

    bam = d / "s.bam"
    _write_reads(str(bam), ref, depth, mapq)
    return {"dir": d, "ref": str(fa), "bam": str(bam), "seq": ref}


@pytest.fixture
def make_bam(dataset, tmp_path):
    def _make(bam=None):
        tools = {
            t: shutil.which(t) or t
            for t in ("samtools", "bcftools", "minimap2", "gatk")
        }
        return Bam(
            bam or dataset["bam"],
            dataset["ref"],
            {"tools": tools, "tmpdir": str(tmp_path)},
        )

    return _make

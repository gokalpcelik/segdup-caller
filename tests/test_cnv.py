import os

import pandas as pd
import pytest
import vcflib

from genecaller.calling import cnv

from conftest import need_tools


@pytest.fixture(autouse=True)
def small_genome(monkeypatch):
    monkeypatch.setattr(cnv, "MIN_AUTOSOME_LEN", 100_000)
    monkeypatch.setattr(cnv, "NORM_BINS_PER_CHROM", 40)
    cnv._DIPLOID_CACHE.clear()


@need_tools
def test_diploid_depth_close_to_truth(make_bam):
    assert cnv.diploid_depth(make_bam()) == pytest.approx(20, rel=0.15)


@need_tools
def test_depth_caller_finds_het_deletion_and_low_mapq(make_bam, tmp_path):
    bam = make_bam()
    out = str(tmp_path / "x.vcf.gz")
    norm = str(tmp_path / "x.norm.tsv")
    cnv.call_cnvscope(bam, out, {"region": "chr1,chr2", "output_norm": norm})
    assert os.path.exists(out + ".tbi")
    df = pd.read_csv(norm, sep="\t")
    assert list(df.columns) == ["contig", "start", "end", "MQ", "DP", "DP0"]

    def mean(c, s, e, col="DP"):
        m = (df.contig == c) & (df.start >= s) & (df.end <= e)
        return df[m][col].mean()

    assert mean("chr1", 10_000, 70_000) == pytest.approx(1.0, rel=0.15)
    assert mean("chr1", 82_000, 98_000) == pytest.approx(0.5, rel=0.2)
    # MAPQ-0 block: confident depth ~0, all-read depth ~1, MQ proxy ~0
    assert mean("chr2", 42_000, 58_000) < 0.1
    assert mean("chr2", 42_000, 58_000, "DP0") == pytest.approx(1.0, rel=0.15)
    assert mean("chr2", 42_000, 58_000, "MQ") < 10

    v = [(r.chrom, r.pos, r.end, r.info["CNT"]) for r in vcflib.VCF(out)]
    assert any(c == "chr1" and p < 98_000 and e > 82_000 and cn == 1 for c, p, e, cn in v)
    assert any(c == "chr1" and cn == 2 for c, _, _, cn in v)


@need_tools
def test_call_depth_populates_segments_and_depths(make_bam, tmp_path):
    bam = make_bam()
    bam.prefix = str(tmp_path / "p")
    cnv.call_depth(bam, str(tmp_path / "c.vcf.gz"), {"region": "chr1:70000-110000"})
    assert bam.depths["chr1"] and bam.segments["chr1"]
    depths = bam.get_depth("chr1:85000-95000")
    assert depths and all(0.3 < d.mean < 0.7 for d in depths)
    assert bam.overall_depth_std >= 0

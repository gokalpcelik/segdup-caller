"""calc_logprob: all-read depth for liftover BAMs and haploid (male chrX) depth scaling."""
import logging

import pytest

from genecaller.bam_process import Bam, DepthSegment


def _bam(mean, mean0, n=20):
    """A Bam holding n 1 kb bins of constant depth on chr1 (depth: diploid autosome = 1)."""
    b = Bam.__new__(Bam)
    b.logger = logging.getLogger("test")
    b.overall_depth_std = 0.1
    segs = [DepthSegment("chr1", i * 1000, (i + 1) * 1000, 60.0, mean, mean0) for i in range(n)]
    b.depths = {"chr1": segs}
    b._depth_coords_cache = {"chr1": ([s.start for s in segs], [s.end for s in segs])}
    return b


def _best_cn(bam, baseline_cn=2, **kw):
    scores = {cn: bam.calc_logprob("chr1:1-20000", cn, baseline_cn=baseline_cn, **kw) for cn in range(0, 7)}
    return max(scores, key=scores.get)


def test_all_reads_counts_mapq0_reads():
    # SMN-like liftover: 4 copies' worth of reads, but only half pass MAPQ >= 30
    b = _bam(mean=1.0, mean0=2.0)
    assert _best_cn(b) == 2
    assert _best_cn(b, all_reads=True) == 4


def test_default_uses_confident_depth():
    b = _bam(mean=1.0, mean0=2.0)
    assert b.calc_logprob("chr1:1-20000", 2) > b.calc_logprob("chr1:1-20000", 4)


@pytest.mark.parametrize("all_reads", [False, True])
def test_haploid_baseline_rescales_depth(all_reads):
    # male non-PAR chrX: one copy has depth 0.5 on the autosomal-diploid scale
    b = _bam(mean=0.5, mean0=0.5)
    assert _best_cn(b, baseline_cn=1, all_reads=all_reads) == 1
    # IKBKG + IKBKGP1 lifted together: two haploid copies -> depth 1.0
    b = _bam(mean=1.0, mean0=1.0)
    assert _best_cn(b, baseline_cn=1, all_reads=all_reads) == 2


def test_diploid_baseline_unchanged():
    b = _bam(mean=1.0, mean0=1.0)
    assert _best_cn(b, baseline_cn=2) == 2
    assert _best_cn(b, baseline_cn=2, all_reads=True) == 2


def _ads(n=30, depth=60):
    return {i: [depth // 2, depth // 2] for i in range(n)}  # 50:50 allele depths (2+2)


def _best_cn_ad(bam, **kw):
    ratios = {3: 2 / 3, 4: 0.5}
    scores = {
        cn: bam.calc_logprob("chr1:1-20000", cn, ads=_ads(), ratio=ratios[cn], all_reads=True, **kw)
        for cn in ratios
    }
    return max(scores, key=scores.get)


def test_min_cv_keeps_ad_evidence_from_being_overridden():
    # lifted depth reads 3.1 copies (reads lost) while the allele ratios say 2+2
    b = _bam(mean=1.55, mean0=1.55)
    assert _best_cn_ad(b) == 3
    assert _best_cn_ad(b, min_cv=0.1) == 4


def test_min_cv_keeps_depth_discrimination():
    # a real 3-copy depth with a 2/3 allele ratio still wins
    b = _bam(mean=1.5, mean0=1.5)
    ads = {i: [40, 20] for i in range(30)}
    s3 = b.calc_logprob("chr1:1-20000", 3, ads=ads, ratio=2 / 3, all_reads=True, min_cv=0.1)
    s4 = b.calc_logprob("chr1:1-20000", 4, ads=ads, ratio=0.5, all_reads=True, min_cv=0.1)
    assert s3 > s4

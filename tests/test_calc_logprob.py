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

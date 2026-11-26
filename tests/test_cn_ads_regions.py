import gzip

import pytest
import yaml

from genecaller.util import IntervalList

DATA = "genecaller/data"

with open(f"{DATA}/hg38/genes.yaml") as fh:
    SMN = yaml.safe_load(fh)["SMN1"]


def _positions(vcf):
    with gzip.open(f"{DATA}/{vcf.removeprefix('data/')}", "rt") as fh:
        return [int(line.split("\t")[1]) - 1 for line in fh if line[0] != "#"]


# c.840 (exon 7) in each paralog's coordinates, 0-based
@pytest.mark.parametrize("paralog, exon7", [(0, 70951945), (1, 70076525)])
def test_smn_cn_ads_regions_keep_each_paralogs_sites(paralog, exon7):
    itv = IntervalList(region=",".join(SMN["config"]["cn_ads_regions"]))
    positions = _positions(SMN["diff_vcf"][paralog])
    kept = [p for p in positions if itv.get("chr5", p, p + 1)]
    assert exon7 in kept
    # the downstream flank sites are the ones the filter is there to drop
    assert 0 < len(kept) < len(positions)


def test_smn_cn_ads_regions_keep_matching_site_counts():
    itv = IntervalList(region=",".join(SMN["config"]["cn_ads_regions"]))
    counts = [
        sum(bool(itv.get("chr5", p, p + 1)) for p in _positions(vcf))
        for vcf in SMN["diff_vcf"]
    ]
    assert abs(counts[0] - counts[1]) <= 1

import glob
import re

import pandas as pd
import pytest
import yaml

from genecaller.calling.cnv import flank_factor
from genecaller.util import IntervalList

REGION = re.compile(r"^[\w.]+:\d+-\d+$")
GENES_YAMLS = sorted(glob.glob("genecaller/data/*/genes.yaml"))


def _flank_entries():
    for path in GENES_YAMLS:
        with open(path) as fh:
            cfg = yaml.safe_load(fh)
        for gene, sect in cfg.items():
            conf = sect.get("config") if isinstance(sect, dict) else None
            if isinstance(conf, dict) and "flank_regions" in conf:
                yield pytest.param(conf["flank_regions"], id=f"{path.split('/')[-2]}-{gene}")


@pytest.mark.parametrize("flank_regions", list(_flank_entries()))
def test_flank_regions_in_genes_yaml_parse(flank_regions):
    assert isinstance(flank_regions, list) and flank_regions
    assert all(REGION.match(r) for r in flank_regions)
    assert IntervalList(region=",".join(flank_regions)).regions


def test_genes_yaml_found():
    assert GENES_YAMLS and list(_flank_entries())


class FakeBam:
    def __init__(self, tsv):
        self.source_cnv = {"k": ("x.vcf.gz", str(tsv))}


@pytest.mark.parametrize(
    "regions",
    ["chr1:100-300,chr1:500-700", ["chr1:100-300", "chr1:500-700"]],
)
def test_flank_factor_accepts_string_or_list(tmp_path, regions):
    tsv = tmp_path / "short_read.chr1.norm.tsv"
    starts = list(range(100, 700, 2))
    pd.DataFrame(
        {"contig": "chr1", "start": starts, "end": [s + 2 for s in starts], "MQ": 60, "DP": 1.1}
    ).to_csv(tsv, sep="\t", index=False)
    factor, n = flank_factor(FakeBam(tsv), regions)
    assert n >= 50
    assert factor == pytest.approx(1.1)

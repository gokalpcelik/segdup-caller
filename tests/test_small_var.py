import os
import stat

from genecaller.calling.small_var import call_variant

from conftest import need_tools


def _fake_gatk(path):
    """A stand-in gatk that records its arguments and writes the -O file."""
    path.write_text(
        "#!/bin/bash\n"
        'echo "$@" > "$0.args"\n'
        'while [ $# -gt 0 ]; do [ "$1" = "-O" ] && out="$2"; shift; done\n'
        'echo "##fileformat=VCFv4.2" | gzip > "$out"\n'
    )
    path.chmod(path.stat().st_mode | stat.S_IEXEC)
    return str(path)


@need_tools
def test_call_variant_invokes_haplotypecaller(make_bam, tmp_path):
    bam = make_bam()
    bam.gatk = _fake_gatk(tmp_path / "gatk")
    out = str(tmp_path / "o.vcf.gz")
    call_variant(bam, out, {"region": "chr1:100-500,chr1:900-950", "ploidy": 3, "given": "/nonexistent.vcf.gz"})
    args = open(bam.gatk + ".args").read()
    assert "HaplotypeCaller" in args and "--sample-ploidy 3" in args
    assert "-L chr1:100-500" in args and "-L chr1:900-950" in args
    assert "--alleles" not in args  # missing given file is ignored
    assert os.path.exists(out)


@need_tools
def test_call_variant_given_alleles(make_bam, tmp_path):
    bam = make_bam()
    bam.gatk = _fake_gatk(tmp_path / "gatk")
    given = tmp_path / "g.vcf.gz"
    given.write_text("x")
    call_variant(bam, str(tmp_path / "o.vcf.gz"), {"region": "chr1:1-100", "given": str(given)})
    assert "GENOTYPE_GIVEN_ALLELES" in open(bam.gatk + ".args").read()

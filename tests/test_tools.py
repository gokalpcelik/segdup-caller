import os
import subprocess

import pysam

from genecaller import tools

from conftest import need_tools


def test_haplotypecaller_cmd_basic():
    cmd = tools.haplotypecaller_cmd(
        "gatk", "ref.fa", "in.bam", "out.vcf.gz", regions="chr1:10-20,chr2:5-9", ploidy=3
    )
    assert cmd.startswith("gatk --java-options")
    assert "HaplotypeCaller" in cmd
    assert "--sample-ploidy 3" in cmd
    assert "-L chr1:10-20 -L chr2:5-9" in cmd
    assert "--alleles" not in cmd and "--dbsnp" not in cmd


def test_haplotypecaller_cmd_given_and_dbsnp():
    cmd = tools.haplotypecaller_cmd(
        "gatk", "r.fa", "a.bam", "o.vcf.gz", given="g.vcf.gz", dbsnp="d.vcf.gz"
    )
    assert "--alleles g.vcf.gz --genotyping-mode GENOTYPE_GIVEN_ALLELES" in cmd
    assert "--dbsnp d.vcf.gz" in cmd


def test_paths_are_quoted():
    cmd = tools.haplotypecaller_cmd("gatk", "my ref.fa", "a b.bam", "o.vcf.gz")
    assert "'my ref.fa'" in cmd and "'a b.bam'" in cmd


def test_realign_cmd_has_no_sentieon():
    cmd = tools.realign_cmd(
        "samtools", "minimap2", "r.fa", "in.bam", "chr1:1-100", "t.fa",
        "@RG\\tID:x", "sr", "out.bam",
    )
    assert "sentieon" not in cmd
    assert "-F 0x400" in cmd and "minimap2" in cmd and "samtools sort" in cmd


def test_parse_bedcov():
    rows = tools.parse_bedcov("chr1\t0\t100\t2000\nchr1\t100\t200\t0\n")
    assert rows == [("chr1", 0, 100, 20.0), ("chr1", 100, 200, 0.0)]


def test_write_region_bed_merges_and_pads(tmp_path):
    bed = tools.write_region_bed("chr1:101-200,chr1:211-300", str(tmp_path / "r.bed"), padding=10)
    assert open(bed).read() == "chr1\t90\t310\n"


@need_tools
def test_extract_and_merge(dataset, tmp_path):
    bed = tools.write_region_bed("chr1:1000-2000", str(tmp_path / "r.bed"), padding=500)
    out = str(tmp_path / "clip.bam")
    cmd = tools.extract_regions("samtools", dataset["ref"], dataset["bam"], bed, out, str(tmp_path))
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
    assert os.path.exists(out + ".bai")
    with pysam.AlignmentFile(out) as f:
        reads = list(f.fetch(until_eof=True))
    assert reads and all(r.reference_name == "chr1" and r.reference_end > 499 and r.reference_start < 2500 for r in reads)
    merged = str(tmp_path / "merged.bam")
    cmd = tools.merge_bams_cmd("samtools", dataset["ref"], [out, out], merged)
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
    with pysam.AlignmentFile(merged) as f:
        assert len(list(f.fetch(until_eof=True))) == 2 * len(reads)
    assert os.path.exists(merged + ".bai")


@need_tools
def test_realign_runs_with_minimap2(dataset, tmp_path):
    import shutil

    if not shutil.which("minimap2"):
        import pytest

        pytest.skip("minimap2 not installed")
    target = str(tmp_path / "t.fa")
    subprocess.run(
        f"samtools faidx {dataset['ref']} chr1:5001-8000 > {target}", shell=True, check=True
    )
    out = str(tmp_path / "re.bam")
    cmd = tools.realign_cmd(
        "samtools", "minimap2", dataset["ref"], dataset["bam"], "chr1:5001-8000",
        target, "@RG\\tID:rg1\\tSM:S1", "sr", out,
    )
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
    with pysam.AlignmentFile(out) as f:
        reads = list(f.fetch(until_eof=True))
    assert len(reads) > 100
    assert sum(not r.is_unmapped for r in reads) / len(reads) > 0.9

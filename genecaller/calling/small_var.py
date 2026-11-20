"""Small-variant calling with GATK HaplotypeCaller.

Functions take a Bam instance as a state-holder for tool paths, reference, and
BAM file paths. They do not mutate the Bam's variant-call-irrelevant state.

Calls are made per region and per copy number (``--sample-ploidy``). With a
``given`` VCF the call genotypes only those alleles (GENOTYPE_GIVEN_ALLELES);
a ``dbsnp`` VCF only annotates IDs.
"""

import os
import subprocess
from typing import Any, Dict

from .. import throttle, tools
from ..logging import get_logger

logger = get_logger(__name__)


def _run(bam: Any, cmd: str) -> None:
    bam.logger.debug(f"Running command: {cmd}")
    with throttle.slot():
        result = subprocess.run(
            cmd, capture_output=True, text=True, shell=True, executable="/bin/bash"
        )
    if result.returncode:
        raise Exception(result.stderr)


def subset_input_vcf(
    bam: Any,
    input_vcf: str,
    region: str,
    output: str,
) -> None:
    cmd = (
        f"set -o pipefail; {bam.bcftools} view -r {region} -f PASS,. {input_vcf} "
        f"| {tools.compress_index_vcf_cmd(bam.bcftools, output)}"
    )
    _run(bam, cmd)


def call_variant(bam: Any, output: str, param: Dict[str, Any]) -> None:
    if bam.use_existing and os.path.exists(output) and os.path.getsize(output):
        return
    if (
        param.get("input_vcf")
        and param.get("ploidy") == 2
        and param.get("bam_type") == "orig"
        and param.get("seq_key") == "short_read"
        and not param.get("given")
    ):
        bam.logger.info(
            f"Using --input_vcf for {param['region']} (skipping the diploid call)"
        )
        subset_input_vcf(bam, param["input_vcf"], param["region"], output)
        return
    call_haplotypecaller(bam, output, param)


def call_haplotypecaller(bam: Any, output: str, param: Dict[str, Any]) -> None:
    if bam.use_existing and os.path.exists(output) and os.path.getsize(output):
        return
    if not output.endswith(".vcf.gz"):
        raise ValueError(f"HaplotypeCaller output must end in .vcf.gz: {output}")
    bam_path = bam.clipped_bam if bam.clipped_bam else bam.bam
    given = param.get("given", "")
    if given and not os.path.exists(given):
        given = ""
    region = param.get("region")
    if region and getattr(bam, "liftover_covered", None) is not None:  # regions-only merged BAM
        from ..util import IntervalList

        for chrom, spans in IntervalList.parse_region(region).items():
            for s, e in spans:
                bam.check_covered(chrom, s, e)
    cmd = tools.haplotypecaller_cmd(
        bam.gatk,
        bam.ref,
        bam_path,
        output,
        regions=region,
        ploidy=param.get("ploidy", 2),
        given=given or None,
        dbsnp=param.get("dbsnp"),
        tmpdir=bam.tmpdir,
        java_options=getattr(bam, "gatk_java_options", tools.DEFAULT_JAVA_OPTIONS),
        min_mapq=getattr(bam, "gatk_min_mapq", tools.DEFAULT_HC_MIN_MAPQ),
        extra=param.get("algo_opt", ""),
    )
    _run(bam, cmd)


# kept so existing call sites and out-of-tree configs keep working
call_dnascope = call_haplotypecaller

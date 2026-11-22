"""Conversion detector: a whole-region ratio shift is a copy-number error, not a fusion."""
import logging
from types import SimpleNamespace

from genecaller.conversion_detector import GeneConversionDetector


def _detector(segment, n_signals, copy_numbers):
    d = GeneConversionDetector.__new__(GeneConversionDetector)
    d.logger = logging.getLogger("test")
    d.copy_numbers = copy_numbers
    d.total_copy_number = sum(copy_numbers)
    d.signals = [object()] * n_signals
    d.min_num_signals = 2
    d.three_prime_terminus = None
    d.five_prime_terminus = None
    d.gene_strand = None
    d.fusion_detection = False
    d.segment_method = "HMMSEG"
    d.segmenter = SimpleNamespace(
        segment=lambda: None, breakpoints=[], segments=[segment]
    )
    return d


def _seg(start, end, n, alt, total, mean):
    return SimpleNamespace(
        start=start, end=end, num_sites=n, allele_counts=[total - alt, alt], mean_signal=mean
    )


def test_whole_region_shift_is_not_a_conversion():
    # called 2+1, but all 35 sites show ~0.53 GBAP1-like reads (2+2)
    d = _detector(_seg(100, 900, 35, 2, 3, 0.534), 35, [2, 1])
    assert d.detect_conversions() == []


def test_localized_conversion_is_kept():
    # 10 of 35 sites carry the extra allele; the segment is not the whole region
    d = _detector(_seg(100, 400, 10, 2, 3, 0.66), 35, [2, 1])
    assert len(d.detect_conversions()) == 1


def test_whole_region_matching_called_ratio_is_kept():
    # whole-region ratio fits the called 2/3 allele fraction
    d = _detector(_seg(100, 900, 35, 2, 3, 0.66), 35, [2, 1])
    assert len(d.detect_conversions()) == 1

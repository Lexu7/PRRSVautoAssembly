from scripts.fill_zero_coverage_with_ref import (
    fill_zero_with_reference_coordinates,
    read_primary_paf,
)


def test_refill_uses_alignment_coordinates_after_consensus_insertion():
    # The query has an insertion after position 2. Its fourth base maps to
    # reference position 3, not reference position 4.
    ref = "ACGT"
    consensus = "ACTGT"
    depths = [10, 10, 10, 0, 10]
    cigar = "2M1I2M"

    sequence, stats = fill_zero_with_reference_coordinates(
        consensus, ref, depths, query_start=0, target_start=0, cigar=cigar
    )

    assert sequence == "ACTGT"
    assert stats["filled_zero_with_ref_aligned"] == 1


def test_refill_masks_zero_depth_insertions_without_reference_coordinate():
    ref = "ACGT"
    consensus = "ACTGT"
    depths = [10, 10, 0, 10, 10]
    cigar = "2M1I2M"

    sequence, stats = fill_zero_with_reference_coordinates(
        consensus, ref, depths, query_start=0, target_start=0, cigar=cigar
    )

    assert sequence == "ACNGT"
    assert stats["zero_without_reference_coordinate"] == 1


def test_primary_paf_segments_are_merged_by_query_coordinate(tmp_path):
    paf = tmp_path / "split.paf"
    paf.write_text(
        "query\t6\t0\t3\t+\tref\t6\t0\t3\t3\t3\t60\ttp:A:P\tcg:Z:3M\n"
        "query\t6\t3\t6\t+\tref\t6\t3\t6\t3\t3\t60\ttp:A:P\tcg:Z:3M\n"
    )

    assert read_primary_paf(paf, 6) == [0, 1, 2, 3, 4, 5]

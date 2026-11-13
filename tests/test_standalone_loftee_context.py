import unittest

from standalone_loftee.context import (
    cds_coordinate,
    deletion_shift_maps_to_exon,
    feature_interval,
)
from standalone_loftee.transcripts import Interval, Transcript


class ContextCoordinateTest(unittest.TestCase):
    def test_insertion_uses_cds_start_before_hgvs_offset(self):
        self.assertEqual(cds_coordinate("276-277", insertion=True), 276)

    def test_deletion_uses_cds_end(self):
        self.assertEqual(cds_coordinate("160-161", insertion=False), 161)

    def test_vcf_anchor_is_removed_from_deletion_interval(self):
        self.assertEqual(feature_interval(27_983_022, "GCTGA", "G"), (27_983_023, 27_983_026))

    def test_deletion_shift_into_mapper_gap_is_rejected(self):
        transcript = Transcript(
            "ENST1", "ENSG1", "22", 1, "protein_coding",
            (Interval(100, 171, 1), Interval(174, 250, 2)),
            (Interval(100, 171, 1), Interval(174, 250, 2)),
            (), (),
        )
        self.assertFalse(deletion_shift_maps_to_exon(transcript, 165, 166, 6))
        self.assertTrue(deletion_shift_maps_to_exon(transcript, 160, 161, 6))


if __name__ == "__main__":
    unittest.main()

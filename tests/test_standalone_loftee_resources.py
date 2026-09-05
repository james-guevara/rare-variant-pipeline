import unittest

from standalone_loftee.resources import LofteeResources
from standalone_loftee.transcripts import Interval, Transcript


class FakeBigWig:
    def stats(self, chrom, start, end, type="mean", exact=False):
        assert chrom == "22"
        assert type == "mean"
        assert exact is False
        return [1.0]


def transcript(strand=1):
    exons = (
        Interval(100, 199, 1),
        Interval(300, 399, 2),
        Interval(500, 599, 3),
    )
    cds = (
        Interval(150, 199, 1),
        Interval(300, 399, 2),
        Interval(500, 550, 3),
    )
    if strand == -1:
        exons = tuple(
            Interval(item.start, item.end, 4 - item.rank) for item in reversed(exons)
        )
        cds = tuple(
            Interval(item.start, item.end, 4 - item.rank) for item in reversed(cds)
        )
    return Transcript("ENST1", "ENSG1", "22", strand, "protein_coding", exons, cds, (), ())


class ResourceTest(unittest.TestCase):
    def setUp(self):
        self.resources = object.__new__(LofteeResources)
        self.resources.gerp = FakeBigWig()

    def test_interval_gerp_uses_one_based_inclusive_coordinates(self):
        self.assertEqual(self.resources.interval_gerp("chr22", 10, 12), 3.0)

    def test_plus_strand_distance_follows_exons_to_stop(self):
        weighted, distance = self.resources.gerp_weighted_distance(transcript(), 180)
        self.assertEqual(weighted, 171.0)
        self.assertEqual(distance, 168)

    def test_minus_strand_distance_follows_reverse_transcript_order(self):
        weighted, distance = self.resources.gerp_weighted_distance(transcript(-1), 520)
        self.assertGreater(weighted, 0)
        self.assertGreater(distance, 0)


if __name__ == "__main__":
    unittest.main()

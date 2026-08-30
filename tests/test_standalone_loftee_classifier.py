import unittest

from standalone_loftee import TranscriptContext, classify


class ClassifierTest(unittest.TestCase):
    def test_end_trunc_matches_pinned_gerp_rule(self):
        result = classify(
            TranscriptContext(
                consequences=("stop_gained",),
                exon_number=(10, 10),
                percentile=0.98,
                gerp_weighted_distance=-60,
                distance_to_stop=140,
                last_exon_coding_length=100,
                conservation_too_short=True,
            )
        )
        self.assertEqual(result.lof, "LC")
        self.assertEqual(result.lof_filter, "END_TRUNC")
        self.assertEqual(
            result.lof_info,
            "PERCENTILE:0.98,GERP_DIST:-60,BP_DIST:140,"
            "DIST_FROM_LAST_EXON:40,50_BP_RULE:FAIL,PHYLOCSF_TOO_SHORT",
        )

    def test_single_exon_is_a_flag_not_a_filter(self):
        result = classify(
            TranscriptContext(
                consequences=("frameshift_variant",),
                exon_number=(1, 1),
                gerp_weighted_distance=0,
                distance_to_stop=500,
                last_exon_coding_length=100,
                conservation_too_short=True,
            )
        )
        self.assertEqual(result.lof, "HC")
        self.assertEqual(result.lof_flags, "SINGLE_EXON")

    def test_splice_filters_follow_perl_order(self):
        result = classify(
            TranscriptContext(
                consequences=("splice_donor_variant",),
                intron_number=(1, 2),
                intron_size=12,
                gc_to_gt_donor=True,
                noncanonical_intron=True,
                five_prime_utr=True,
                ancestral_allele=True,
            )
        )
        self.assertEqual(result.lof, "LC")
        self.assertEqual(
            result.lof_filter,
            "SMALL_INTRON,GC_TO_GT_DONOR,5UTR_SPLICE,ANC_ALLELE",
        )
        self.assertEqual(result.lof_flags, "NON_CAN_SPLICE")

    def test_explicit_utr_splice_consequence_is_low_confidence(self):
        result = classify(
            TranscriptContext(
                consequences=("splice_donor_variant", "5_prime_UTR_variant"),
                five_prime_utr=False,
            )
        )
        self.assertEqual(result.lof, "LC")
        self.assertEqual(result.lof_filter, "5UTR_SPLICE")

    def test_non_candidate_is_unannotated(self):
        result = classify(TranscriptContext(consequences=("missense_variant",)))
        self.assertEqual(result.lof, "")


if __name__ == "__main__":
    unittest.main()

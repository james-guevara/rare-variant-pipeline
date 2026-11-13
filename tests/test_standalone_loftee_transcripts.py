import tempfile
import unittest
from pathlib import Path

from standalone_loftee.transcripts import TranscriptStore, build_transcript_db


GTF = """\
##gtf-version 3
22\ttest\ttranscript\t100\t500\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1"; transcript_biotype "protein_coding";
22\ttest\texon\t100\t199\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1"; exon_number "1";
22\ttest\tCDS\t150\t199\t.\t+\t0\tgene_id "ENSG1"; transcript_id "ENST1"; exon_number "1";
22\ttest\texon\t300\t399\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1"; exon_number "2";
22\ttest\tCDS\t300\t350\t.\t+\t2\tgene_id "ENSG1"; transcript_id "ENST1"; exon_number "2";
22\ttest\ttranscript\t600\t900\t.\t-\t.\tgene_id "ENSG2"; transcript_id "ENST2"; transcript_biotype "protein_coding";
22\ttest\texon\t800\t900\t.\t-\t.\tgene_id "ENSG2"; transcript_id "ENST2"; exon_number "1";
22\ttest\tCDS\t820\t880\t.\t-\t0\tgene_id "ENSG2"; transcript_id "ENST2"; exon_number "1";
22\ttest\texon\t600\t699\t.\t-\t.\tgene_id "ENSG2"; transcript_id "ENST2"; exon_number "2";
22\ttest\tCDS\t620\t699\t.\t-\t2\tgene_id "ENSG2"; transcript_id "ENST2"; exon_number "2";
"""


class TranscriptStoreTest(unittest.TestCase):
    def test_builds_strand_aware_exons_introns_and_cds(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            gtf, database = tmp / "test.gtf", tmp / "transcripts.sqlite"
            gtf.write_text(GTF)
            build_transcript_db(gtf, database, seqname="22")
            with TranscriptStore(database) as store:
                plus = store.get("ENST1.9")
                minus = store.get("ENST2")
            self.assertIsNotNone(plus)
            self.assertIsNotNone(minus)
            assert plus is not None and minus is not None
            self.assertEqual(plus.cds_length, 101)
            self.assertEqual([(x.start, x.end) for x in plus.introns], [(200, 299)])
            self.assertEqual(plus.last_exon_coding_length(), 50)
            self.assertEqual([x.rank for x in minus.exons], [1, 2])
            self.assertEqual(minus.last_exon_coding_length(), 79)


if __name__ == "__main__":
    unittest.main()

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "build_target_bed.py"


class BuildTargetBedTest(unittest.TestCase):
    def test_selects_by_stable_gene_id_and_merges_padded_exons(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            genebayes = tmp / "genebayes.tsv"
            with genebayes.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["ensg", "post_mean"], delimiter="\t")
                writer.writeheader()
                writer.writerows([
                    {"ensg": "ENSG000001.9", "post_mean": "0.18"},
                    {"ensg": "ENSG000002", "post_mean": "0.029"},
                ])
            gtf = tmp / "genes.gtf"
            gtf.write_text(
                'chr22\ttest\texon\t101\t110\t.\t+\t.\tgene_id "ENSG000001.1";\n'
                'chr22\ttest\texon\t115\t120\t.\t+\t.\tgene_id "ENSG000001.1";\n'
                'chr22\ttest\texon\t500\t510\t.\t+\t.\tgene_id "ENSG000002";\n'
                'chr21\ttest\texon\t10\t20\t.\t+\t.\tgene_id "ENSG000001";\n'
            )
            output = tmp / "target.bed"
            completed = subprocess.run(
                [sys.executable, str(SCRIPT), "--genebayes", str(genebayes),
                 "--gtf", str(gtf), "--output", str(output), "--padding", "5",
                 "--chroms", "chr22"],
                check=True, capture_output=True, text=True,
            )
            self.assertEqual(output.read_text(), "chr22\t95\t125\n")
            self.assertIn("selected_genes=1", completed.stderr)
            self.assertIn("matched_genes=1", completed.stderr)

    def test_can_add_chr_prefix_for_chr_named_vcf(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            genebayes = tmp / "genebayes.tsv"
            genebayes.write_text("ensg\tpost_mean\nENSG000001\t0.03\n")
            gtf = tmp / "genes.gtf"
            gtf.write_text(
                '22\ttest\texon\t101\t110\t.\t+\t.\tgene_id "ENSG000001";\n'
            )
            output = tmp / "target.bed"
            subprocess.run(
                [sys.executable, str(SCRIPT), "--genebayes", str(genebayes),
                 "--gtf", str(gtf), "--output", str(output), "--padding", "0",
                 "--chroms", "22", "--add-chr-prefix"],
                check=True, capture_output=True, text=True,
            )
            self.assertEqual(output.read_text(), "chr22\t100\t110\n")

    def test_optional_gene_allowlist_limits_selected_genes(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            genebayes = tmp / "genebayes.tsv"
            genebayes.write_text(
                "ensg\tpost_mean\nENSG000001\t0.20\nENSG000002\t0.10\n"
            )
            allowlist = tmp / "genes.txt"
            allowlist.write_text("ENSG000002\n")
            gtf = tmp / "genes.gtf"
            gtf.write_text(
                'chr22\ttest\texon\t101\t110\t.\t+\t.\tgene_id "ENSG000001";\n'
                'chr22\ttest\texon\t201\t210\t.\t+\t.\tgene_id "ENSG000002";\n'
            )
            output = tmp / "target.bed"
            subprocess.run(
                [sys.executable, str(SCRIPT), "--genebayes", str(genebayes),
                 "--gtf", str(gtf), "--output", str(output), "--padding", "0",
                 "--gene-ids", str(allowlist)],
                check=True, capture_output=True, text=True,
            )
            self.assertEqual(output.read_text(), "chr22\t200\t210\n")


if __name__ == "__main__":
    unittest.main()

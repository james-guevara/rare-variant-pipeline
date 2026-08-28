import csv, hashlib, json, subprocess, sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts/gather_rare_burdens.py"

def package(root, chrom, lof, miss=None):
    root.mkdir(); files = {}
    for logical, name, text in (("plof_counts", "12.plof-per-sample-counts.tsv", lof), ("missense_counts", "12.missense-per-sample-counts.tsv", miss)):
        if text is None: continue
        path = root / name; path.write_text(text)
        files[logical] = {"file": name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "size_bytes": path.stat().st_size}
    (root / "targeted-output-manifest.json").write_text(json.dumps({"schema_version": 1, "cohort": "c", "chromosome": chrom, "status": "SUCCEEDED", "files": files}))

def test_gathers_and_zero_fills_only_complete_tiers(tmp_path):
    samples = tmp_path / "samples.tsv"; samples.write_text("FID\tIID\tSEX\nF1\tS1\tF\nF2\tS2\tM\n")
    p1, px = tmp_path / "chr1", tmp_path / "chrX"
    package(p1, "chr1", "SAMPLE\tlof_t1\tlof_t2\nS1\t1\t0\n", "SAMPLE\tmiss_t1\tmiss_t2\tmiss_t3\tmiss_t4\nS1\t0\t1\t0\t2\n")
    package(px, "chrX", "SAMPLE\tlof_t1\tlof_t2\nS2\t0\t1\n")
    out, strata = tmp_path / "rare.tsv", tmp_path / "strata.tsv"
    subprocess.run([sys.executable, str(SCRIPT), "--sample-manifest", str(samples), "--package", str(p1), "--package", str(px), "--expected-chromosomes", "chr1,chrX", "--output", str(out), "--strata-output", str(strata)], check=True)
    rows = list(csv.DictReader(out.open(), delimiter="\t"))
    assert rows[0]["lof_t1"] == "1" and rows[1]["lof_t2"] == "1"
    assert rows[1]["lof_t1"] == "0"
    assert rows[0]["miss_t1"] == ""

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

def test_gathers_completed_run_directly_with_sex_qc_roster(tmp_path):
    qc = tmp_path / "sex-qc.tsv"
    qc.write_text("sample_id\tinferred_karyotype\nS1\tXX-like\nS2\tambiguous\n")
    run = tmp_path / "run"; chrom = run / "chr1"; chrom.mkdir(parents=True)
    (chrom / "_SUCCESS").touch()
    (chrom / "12.plof-per-sample-counts.tsv").write_text(
        "SAMPLE\tlof_t1\tlof_t2\nS1\t1\t0\n"
    )
    (chrom / "12.missense-per-sample-counts.tsv").write_text(
        "SAMPLE\tmiss_t1\tmiss_t2\tmiss_t3\tmiss_t4\nS1\t0\t1\t0\t2\n"
    )
    out, strata = tmp_path / "rare.tsv", tmp_path / "strata.tsv"
    subprocess.run([
        sys.executable, str(SCRIPT), "--sex-qc", str(qc), "--run-base", str(run),
        "--expected-chromosomes", "chr1", "--output", str(out),
        "--strata-output", str(strata),
    ], check=True)
    rows = list(csv.DictReader(out.open(), delimiter="\t"))
    assert rows[0]["SEX"] == "F" and rows[0]["lof_t1"] == "1"
    assert rows[1]["SEX"] == "" and rows[1]["miss_t4"] == "0"

def test_gathers_nextflow_staged_package_root(tmp_path):
    samples = tmp_path / "samples.tsv"
    samples.write_text("FID\tIID\tSEX\nF1\tS1\tF\n")
    staged = tmp_path / "packages"
    staged.mkdir()
    package(
        staged / "package01",
        "chr22",
        "SAMPLE\tlof_t1\tlof_t2\nS1\t2\t1\n",
    )
    out, strata = tmp_path / "rare.tsv", tmp_path / "strata.tsv"
    subprocess.run([
        sys.executable, str(SCRIPT), "--sample-manifest", str(samples),
        "--package-root", str(staged), "--expected-chromosomes", "chr22",
        "--output", str(out), "--strata-output", str(strata),
    ], check=True)
    rows = list(csv.DictReader(out.open(), delimiter="\t"))
    assert rows[0]["lof_t1"] == "2" and rows[0]["lof_t2"] == "1"

def test_gathers_symlinked_nextflow_package_root(tmp_path):
    samples = tmp_path / "samples.tsv"
    samples.write_text("FID\tIID\tSEX\nF1\tS1\tF\n")
    source = tmp_path / "chr22"
    package(source, "chr22", "SAMPLE\tlof_t1\tlof_t2\nS1\t2\t1\n")
    staged = tmp_path / "packages"
    staged.mkdir()
    (staged / "package01").symlink_to(source, target_is_directory=True)
    out, strata = tmp_path / "rare.tsv", tmp_path / "strata.tsv"
    subprocess.run([
        sys.executable, str(SCRIPT), "--sample-manifest", str(samples),
        "--package-root", str(staged), "--expected-chromosomes", "chr22",
        "--output", str(out), "--strata-output", str(strata),
    ], check=True)
    rows = list(csv.DictReader(out.open(), delimiter="\t"))
    assert rows[0]["lof_t1"] == "2" and rows[0]["lof_t2"] == "1"

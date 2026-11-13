import csv
import json
import subprocess
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "inspect_cohort_vcfs.py"


def fixture(tmp_path, samples=("S1", "S2")):
    vcf = tmp_path / "input.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=1,length=248956422>\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"GQ\">\n"
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"DP\">\n"
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"AD\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(samples) + "\n"
    )
    plan = tmp_path / "plan.tsv"
    plan.write_text(
        "chromosome\tcontig\tcontig_length\tinput_vcf\tpreparation_state\n"
        f"chr1\tauto\t248956422\t{vcf}\tPENDING_DERIVED_RESOURCES\n"
    )
    manifest = tmp_path / "samples.tsv"
    manifest.write_text("FID\tIID\tSEX\n" + "".join(f"\t{s}\t\n" for s in samples))
    return plan, manifest


def test_resolves_bare_contig_and_marks_ready_for_zarr(tmp_path):
    plan, samples = fixture(tmp_path)
    inspected, report = tmp_path / "inspected.tsv", tmp_path / "report.json"
    subprocess.run([
        sys.executable, str(SCRIPT), "--preparation-plan", str(plan),
        "--sample-manifest", str(samples), "--output-plan", str(inspected),
        "--report", str(report),
    ], check=True)
    row = next(csv.DictReader(inspected.open(), delimiter="\t"))
    assert row["contig"] == "1"
    assert row["preparation_state"] == "READY_FOR_ZARR"
    assert json.loads(report.read_text())["samples"] == 2


def test_rejects_vcf_psam_sample_mismatch(tmp_path):
    plan, samples = fixture(tmp_path, ("S1", "S2"))
    samples.write_text("FID\tIID\tSEX\n\tS1\t\n\tS3\t\n")
    result = subprocess.run([
        sys.executable, str(SCRIPT), "--preparation-plan", str(plan),
        "--sample-manifest", str(samples), "--output-plan", str(tmp_path / "out"),
        "--report", str(tmp_path / "report"),
    ], capture_output=True, text=True)
    assert result.returncode != 0
    assert "VCF/PSAM sample mismatch" in result.stderr


def test_accepts_localized_allele_depth_fields(tmp_path):
    plan, samples = fixture(tmp_path)
    vcf = Path(next(csv.DictReader(plan.open(), delimiter="\t"))["input_vcf"])
    text = vcf.read_text().replace(
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">',
        '##FORMAT=<ID=LAD,Number=.,Type=Integer,Description="LAD">\n'
        '##FORMAT=<ID=LAA,Number=.,Type=Integer,Description="LAA">',
    )
    vcf.write_text(text)
    inspected, report = tmp_path / "inspected.tsv", tmp_path / "report.json"
    subprocess.run([
        sys.executable, str(SCRIPT), "--preparation-plan", str(plan),
        "--sample-manifest", str(samples), "--output-plan", str(inspected),
        "--report", str(report),
    ], check=True)
    assert json.loads(report.read_text())["chromosomes"][0]["allele_depth_format"] == "LAD+LAA"

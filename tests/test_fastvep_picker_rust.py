import subprocess
from pathlib import Path


ROOT = Path(__file__).parents[1]
CRATE = ROOT / "rust" / "fastvep-picker" / "Cargo.toml"
FIXTURE = ROOT / "tests" / "fixtures" / "fastvep-picker"


def picker_command(executable, output):
    return [
        str(executable),
        "--fastvep",
        str(FIXTURE / "input.vcf"),
        "--transcript-priority",
        str(FIXTURE / "transcript-priority.tsv"),
        "--consequence-ranks",
        str(FIXTURE / "consequence-ranks.tsv"),
        "--output",
        str(output),
    ]


def test_rust_picker_is_byte_identical_to_python_oracle(tmp_path):
    subprocess.run(
        ["cargo", "build", "--quiet", "--release", "--manifest-path", str(CRATE)],
        check=True,
    )
    rust_output = tmp_path / "rust.tsv"
    python_output = tmp_path / "python.tsv"
    rust_picker = CRATE.parent / "target" / "release" / "fastvep-picker"

    subprocess.run(picker_command(rust_picker, rust_output), check=True)
    subprocess.run(
        [
            "python3",
            str(ROOT / "scripts" / "pick_fastvep_consequences.py"),
            *picker_command("unused", python_output)[1:],
        ],
        check=True,
    )

    assert rust_output.read_bytes() == python_output.read_bytes()


def test_rust_picker_accepts_streamed_vcf(tmp_path):
    subprocess.run(
        ["cargo", "build", "--quiet", "--release", "--manifest-path", str(CRATE)],
        check=True,
    )
    output = tmp_path / "streamed.tsv"
    rust_picker = CRATE.parent / "target" / "release" / "fastvep-picker"
    command = picker_command(rust_picker, output)
    command[2] = "-"

    subprocess.run(command, input=(FIXTURE / "input.vcf").read_bytes(), check=True)

    oracle = tmp_path / "oracle.tsv"
    subprocess.run(
        [
            "python3",
            str(ROOT / "scripts" / "pick_fastvep_consequences.py"),
            *picker_command("unused", oracle)[1:],
        ],
        check=True,
    )
    assert output.read_bytes() == oracle.read_bytes()

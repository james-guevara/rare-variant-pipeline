#!/usr/bin/env python3
from __future__ import annotations

"""
Verify Part 6 completion and gather genotype chunks into a single file per chromosome.

Verification steps:
1. Compute expected chunks from Part 4 BED (ceil(lines / regions_per_chunk))
2. Verify Part 5 scatter file count matches
3. Verify Part 6 trace COMPLETED/CACHED count matches
4. Verify genotype file count matches

If all pass, concatenate chunks into {chrom}.family_genotypes.tsv
"""

import argparse
import csv
import shutil
import sys
from pathlib import Path


def count_lines(filepath: Path) -> int:
    """Count lines in a file."""
    with open(filepath) as f:
        return sum(1 for _ in f)


def compute_expected_chunks(bed_file: Path, regions_per_chunk: int) -> tuple[int, int]:
    """Compute expected chunk count from BED file lines."""
    bed_lines = count_lines(bed_file)
    expected = (bed_lines + regions_per_chunk - 1) // regions_per_chunk
    return bed_lines, expected


def count_scatter_files(scatter_dir: Path, chrom: str) -> int:
    """Count Part 5 scatter BED files for a chromosome."""
    pattern = f"{chrom}.chunk_*.bed"
    return len(list(scatter_dir.glob(pattern)))


def count_trace_completions(trace_file: Path, chrom: str) -> int:
    """Count COMPLETED/CACHED tasks for a chromosome in trace file."""
    count = 0
    # Pattern to match: "FAMILY_QUERY (chr1_chr1.chunk_00)"
    search_pattern = f"({chrom}_{chrom}.chunk_"

    with open(trace_file, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            status = row.get('status', '')
            name = row.get('name', '')
            if status in ('COMPLETED', 'CACHED') and search_pattern in name:
                count += 1
    return count


def count_genotype_files(genotypes_dir: Path, chrom: str) -> int:
    """Count genotype TSV files for a chromosome."""
    pattern = f"{chrom}.chunk_*.genotypes.tsv"
    return len(list(genotypes_dir.glob(pattern)))


def get_sorted_genotype_files(genotypes_dir: Path, chrom: str) -> list[Path]:
    """Get genotype files sorted by chunk number."""
    pattern = f"{chrom}.chunk_*.genotypes.tsv"
    files = list(genotypes_dir.glob(pattern))

    def chunk_number(p: Path) -> int:
        # Extract chunk number from "chr1.chunk_05.genotypes.tsv"
        name = p.stem  # "chr1.chunk_05.genotypes"
        parts = name.split('.')
        for part in parts:
            if part.startswith('chunk_'):
                return int(part.replace('chunk_', ''))
        return 0

    return sorted(files, key=chunk_number)


def concatenate_genotypes(genotype_files: list[Path], output_file: Path) -> int:
    """Concatenate genotype files using bulk copy, keeping header from first file only."""
    with open(output_file, 'w') as out:
        for i, gf in enumerate(genotype_files):
            with open(gf) as f:
                header = f.readline()
                if i == 0:
                    # Write header from first file
                    out.write(header)
                # Bulk copy remaining content (skips header for all files)
                shutil.copyfileobj(f, out)

    # Count lines in output
    return count_lines(output_file)


def main():
    parser = argparse.ArgumentParser(description='Verify and gather genotype chunks')
    parser.add_argument('--chrom', required=True, help='Chromosome name (e.g., chr1)')
    parser.add_argument('--reformat-dir', required=True, type=Path, help='Part 4 reformat output directory')
    parser.add_argument('--scatter-dir', required=True, type=Path, help='Part 5 scatter output directory')
    parser.add_argument('--genotypes-dir', required=True, type=Path, help='Part 6 genotypes output directory')
    parser.add_argument('--trace-file', required=True, type=Path, help='Part 6 trace file')
    parser.add_argument('--output', required=True, type=Path, help='Output file path')
    parser.add_argument('--regions-per-chunk', type=int, default=1000, help='Regions per chunk (default: 1000)')

    args = parser.parse_args()

    print(f"=== Verifying {args.chrom} ===")
    errors = []

    # 1. Compute expected chunks from Part 4 BED
    bed_file = args.reformat_dir / f"{args.chrom}.consequential.bed"
    if not bed_file.exists():
        print(f"ERROR: Missing BED file: {bed_file}", file=sys.stderr)
        sys.exit(1)

    bed_lines, expected = compute_expected_chunks(bed_file, args.regions_per_chunk)
    print(f"Part 4 BED lines: {bed_lines}")
    print(f"Expected chunks (ceil({bed_lines} / {args.regions_per_chunk})): {expected}")

    # 2. Verify Part 5 scatter output
    scatter_count = count_scatter_files(args.scatter_dir, args.chrom)
    print(f"Part 5 scatter files: {scatter_count}")
    if expected != scatter_count:
        errors.append(f"Part 5 scatter mismatch! Expected {expected}, got {scatter_count}")

    # 3. Verify Part 6 trace completions
    completed = count_trace_completions(args.trace_file, args.chrom)
    print(f"Part 6 completed tasks: {completed}")
    if expected != completed:
        errors.append(f"Part 6 trace mismatch! Expected {expected}, got {completed}")

    # 4. Verify genotype files exist
    genotype_count = count_genotype_files(args.genotypes_dir, args.chrom)
    print(f"Genotype files: {genotype_count}")
    if expected != genotype_count:
        errors.append(f"Genotype file mismatch! Expected {expected}, got {genotype_count}")

    # Report errors and exit if any
    if errors:
        print("\n=== VERIFICATION FAILED ===", file=sys.stderr)
        for err in errors:
            print(f"ERROR: {err}", file=sys.stderr)
        sys.exit(1)

    # All checks passed - concatenate
    print("\n=== All checks passed, concatenating ===")
    genotype_files = get_sorted_genotype_files(args.genotypes_dir, args.chrom)
    total_lines = concatenate_genotypes(genotype_files, args.output)
    print(f"Created {args.output} with {total_lines} lines")


if __name__ == '__main__':
    main()

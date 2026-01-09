#!/usr/bin/env python3
"""
Extract VCF records for entire families when any family member carries a variant.

This script efficiently queries VCF files for family-based variant analysis. When any
family member has an alternate allele at a variant site, genotypes for all family
members are extracted for that variant.
"""
import sys
import argparse
from pathlib import Path
from cyvcf2 import VCF
import numpy as np
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(
        description="Extract VCF records for entire families if ANY member has an ALT genotype."
    )
    p.add_argument("--vcf", required=True, help="Input VCF (bgzipped + indexed)")
    p.add_argument("--ped", required=True, help="PED file (Col1=FamID, Col2=IID)")
    p.add_argument("--out", required=True, help="Output TSV")
    p.add_argument(
        "--region",
        required=True,
        help="BED file (TAB-delimited): chrom<TAB>start0<TAB>end (0-based, end-exclusive)",
    )
    return p.parse_args()

def load_pedigree(ped_path: str) -> dict:
    """
    Load pedigree file and create sample-to-family mapping.

    Args:
        ped_path: Path to PED file (tab-delimited, columns: FamID, IID, ...)

    Returns:
        Dictionary mapping sample IDs to family IDs

    Raises:
        SystemExit: If PED file not found
    """
    sample_to_fam = {}
    try:
        with open(ped_path, "r") as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                cols = line.strip().split()
                if len(cols) < 2:
                    continue
                fam_id, sample_id = cols[0], cols[1]

                # skip typical header row like: FAM IID ...
                if fam_id.lower() in {"fam", "famid", "family", "familyid"} and sample_id.lower() in {"iid", "sample", "sampleid"}:
                    continue

                sample_to_fam[sample_id] = fam_id
    except FileNotFoundError:
        sys.exit(f"ERROR: Could not find PED file at {ped_path}")
    return sample_to_fam


def normalize_chrom(chrom: str, vcf_seqnames: set) -> str | None:
    """
    Normalize chromosome name to match VCF contig names.

    Handles chr-prefix differences (e.g., 'chr1' vs '1').

    Args:
        chrom: Chromosome name from BED file
        vcf_seqnames: Set of sequence names from VCF header

    Returns:
        Normalized chromosome name if match found, None otherwise
    """
    if chrom in vcf_seqnames:
        return chrom
    if not chrom.startswith("chr") and f"chr{chrom}" in vcf_seqnames:
        return f"chr{chrom}"
    if chrom.startswith("chr"):
        stripped = chrom[3:]
        if stripped in vcf_seqnames:
            return stripped
    return None


def read_bed_queries(bed_path: str, vcf_seqnames: set) -> list[str]:
    """
    Parse BED file and convert to cyvcf2 query strings.

    BED format is 0-based, half-open [start, end). cyvcf2 uses 1-based, inclusive coords.

    Args:
        bed_path: Path to BED file (tab-delimited: chrom, start0, end)
        vcf_seqnames: Set of valid chromosome names from VCF

    Returns:
        List of query strings in format "chrom:start1-end_inclusive"

    Raises:
        SystemExit: If BED file invalid or no usable regions found
    """
    p = Path(bed_path)
    if not p.is_file():
        sys.exit(f"ERROR: --region must be a BED file path: {bed_path}")

    queries = []
    with p.open("r") as f:
        for ln, line in enumerate(f, start=1):
            if not line.strip(): continue
            if line.startswith("#") or line.startswith("track") or line.startswith("browser"): continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3: sys.exit(f"ERROR: BED must be TAB-delimited with >=3 cols at line {ln}: {line.strip()}")

            raw_chrom, start0_str, end_str = cols[0], cols[1], cols[2]
            norm_chrom = normalize_chrom(raw_chrom, vcf_seqnames)
            if not norm_chrom:
                print(f"WARNING: skipping BED line {ln} (chrom not in VCF contigs): {raw_chrom}", file=sys.stderr)
                continue

            try:
                start0 = int(start0_str)
                end_excl = int(end_str)
            except ValueError: sys.exit(f"ERROR: Non-integer start/end in BED at line {ln}: {line.strip()}")
            if end_excl <= start0: sys.exit(f"ERROR: Invalid BED interval (end <= start) at line {ln}: {line.strip()}")

            # cyvcf2 region uses 1-based inclusive coords: chrom:start-end
            start1 = start0 + 1
            end_inclusive = end_excl  # because BED end is exclusive; inclusive end = end_excl
            queries.append(f"{norm_chrom}:{start1}-{end_inclusive}")

    if not queries:
        sys.exit("ERROR: No usable regions found in BED (or none matched VCF contigs).")
    return queries


def main():
    """
    Main execution function.

    Workflow:
    1. Load pedigree and build family-sample mappings
    2. Open VCF and validate samples against pedigree
    3. Parse BED file to create genomic region queries
    4. For each region:
       a. Identify variants where any family member has ALT allele
       b. Extract genotypes for ALL family members at those sites
    5. Write results to TSV

    Output columns: #CHROM, POS0, END, POS, REF, ALT, SAMPLE, GT, GQ, DP, AD, FAMILY
    """
    args = parse_args()

    sample_to_fam = load_pedigree(args.ped)

    try: vcf = VCF(args.vcf, strict_gt=True)
    except Exception as e: sys.exit(f"ERROR: opening VCF failed: {e}")

    vcf_samples = list(vcf.samples)
    vcf_seqnames = set(vcf.seqnames)

    # HARD REQUIREMENT: every VCF sample must be in PED IID column
    missing_in_ped = sorted(set(vcf_samples) - set(sample_to_fam.keys()))
    if missing_in_ped:
        preview = ", ".join(missing_in_ped[:25])
        more = "" if len(missing_in_ped) <= 25 else f" ... (+{len(missing_in_ped)-25} more)"
        sys.exit(f"ERROR: {len(missing_in_ped)} VCF samples missing from PED IID: {preview}{more}")

    extra_in_ped = sorted(set(sample_to_fam.keys()) - set(vcf_samples))
    if extra_in_ped:
        preview = ", ".join(extra_in_ped[:25])
        more = "" if len(extra_in_ped) <= 25 else f" ... (+{len(extra_in_ped)-25} more)"
        print(f"WARNING: {len(extra_in_ped)} PED samples not present in VCF: {preview}{more}", file=sys.stderr)


    # For each VCF sample index, what family is it in? (family ID for each sample in the VCF’s sample order)
    family_id_by_sample_index = []
    for sample_id in vcf_samples: family_id_by_sample_index.append(sample_to_fam[sample_id])
    # For each family, which VCF sample indices belong to it?
    sample_indices_by_family_id = defaultdict(list)
    for sample_index, family_id in enumerate(family_id_by_sample_index):
        sample_indices_by_family_id[family_id].append(sample_index)


    queries = read_bed_queries(args.region, vcf_seqnames)

    with open(args.out, "w") as out:
        header = ["#CHROM", "POS0", "END", "POS", "REF", "ALT", "SAMPLE", "GT", "GQ", "DP", "AD", "FAMILY"]
        out.write("\t".join(header) + "\n")

        line_count = 0

        for region_str in queries:
            try: iterator = vcf(region_str)
            except Exception as e: sys.exit(f"ERROR: querying region {region_str} failed: {e}")

            for variant in iterator:
                # PASS 1: identify families with at least one ALT genotype (HET or HOM_ALT)
                # gt_types: 0=HOM_REF, 1=HET, 2=UNKNOWN, 3=HOM_ALT
                carrier_idx = np.nonzero((variant.gt_types == 1) | (variant.gt_types == 3))[0]
                if carrier_idx.size == 0: continue
                # Identify families with at least one carrier
                hit_fams = {family_id_by_sample_index[i] for i in carrier_idx}

                # Variant-level fields
                chrom = variant.CHROM
                pos1 = int(variant.POS)
                pos0 = pos1 - 1
                ref = variant.REF
                alt = ",".join(variant.ALT) if variant.ALT else "."
                end_excl = pos0 + len(ref)  # BED-style end (exclusive) for the REF span

                all_genotypes = variant.genotypes
                ad_list = variant.format("AD")
                dp_list = variant.format("DP")
                gq_list = variant.format("GQ")

                # PASS 2: output all members of families with at least one carrier
                for family_id in hit_fams:
                    for sample_index in sample_indices_by_family_id[family_id]:
                        sample_id = vcf_samples[sample_index]

                        # Parse genotype
                        gt_data = all_genotypes[sample_index]  # [a1, a2, phased]
                        if not gt_data:
                            gt_str = "."
                        else:
                            phased = bool(gt_data[-1])
                            alleles = gt_data[:-1]
                            sep = "|" if phased else "/"
                            allele_strs = [("." if a == -1 else str(a)) for a in alleles]
                            gt_str = allele_strs[0] if len(allele_strs) == 1 else sep.join(allele_strs)

                        # Parse AD (allelic depths)
                        if ad_list is not None and ad_list[sample_index] is not None:
                            try:
                                ad_str = ",".join(map(str, ad_list[sample_index]))
                            except TypeError:
                                ad_str = str(ad_list[sample_index])
                        else:
                            ad_str = "."

                        # Parse DP (total depth)
                        if dp_list is not None and dp_list[sample_index] is not None:
                            try:
                                dp_val = str(dp_list[sample_index][0])
                            except (TypeError, IndexError):
                                dp_val = str(dp_list[sample_index])
                        else:
                            dp_val = "."

                        # Parse GQ (genotype quality)
                        if gq_list is not None and gq_list[sample_index] is not None:
                            try:
                                gq_val = str(gq_list[sample_index][0])
                            except (TypeError, IndexError):
                                gq_val = str(gq_list[sample_index])
                        else:
                            gq_val = "."

                        out.write(
                            f"{chrom}\t{pos0}\t{end_excl}\t{pos1}\t{ref}\t{alt}\t"
                            f"{sample_id}\t{gt_str}\t{gq_val}\t{dp_val}\t{ad_str}\t{family_id}\n"
                        )
                        line_count += 1

    print(f"Success! Written {line_count} lines to {args.out}")


if __name__ == "__main__":
    main()


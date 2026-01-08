"""
Resolve multi-allelic genotypes from VCF query output.

For samples with non-reference genotypes (e.g., GT=1/2), this script:
- Explodes each allele into a separate row
- Flags compound heterozygotes
- Splits allelic depths (AD) into ref and alt components
"""
import polars as pl
import argparse
import sys

def resolve_genotypes(tsv_path, output_path):
    """
    Process genotype TSV to resolve multi-allelic sites.

    Args:
        tsv_path: Input TSV from bcftools query
        output_path: Output TSV path

    Output columns include:
        - ALT_resolved: Specific ALT allele for this row
        - is_compound_het: Flag for compound heterozygotes
        - AD_ref, AD_alt: Split allelic depths
    """
    try:
        # 1. Load Data
        df = pl.read_csv(tsv_path, separator='\t', infer_schema_length=0)

        # 2. Extract Allele Indices (Non-Reference only)
        # Regex finds '1', '2' etc. (skips '0')
        df = df.with_columns(
            pl.col("GT")
            .str.extract_all(r"([1-9][0-9]*)") 
            .alias("allele_indices_str")
        )

        # 3. Flag COMPOUND HETEROZYGOTES (Before Exploding)
        df = df.with_columns(
            (
                (pl.col("allele_indices_str").list.len() > 1) & 
                (pl.col("allele_indices_str").list.n_unique() > 1)
            ).alias("is_compound_het")
        )

        # 4. Filter & Explode
        df = df.drop_nulls("allele_indices_str")
        df = df.explode("allele_indices_str")
        
        # Calculate 0-based lookup index (VCF Allele 1 -> Index 0)
        df = df.with_columns(
            pl.col("allele_indices_str")
            .cast(pl.Int32)
            .sub(1)
            .alias("lookup_index")
        )

        # 5. Resolve ALT
        df = df.with_columns(
            pl.col("ALT")
            .str.split(",") 
            .list.get(pl.col("lookup_index"))
            .alias("ALT_resolved")
        )

        # 6. Resolve AD (Ref + Specific Alt)
        # Logic: AD list is [Ref, Alt1, Alt2...]
        # Ref is always index 0.
        # Specific Alt is index (lookup_index + 1).
        
        # We first split AD into a list
        df = df.with_columns(
            pl.col("AD")
            .str.split(",")
            .alias("AD_list")
        )
        
        # Extract Ref and Specific Alt depths
        df = df.with_columns(
            pl.col("AD_list").list.get(0).alias("AD_ref"),
            pl.col("AD_list").list.get(pl.col("lookup_index") + 1).alias("AD_alt")
        )
        
        # 7. Clean Up
        # Rename original columns to *_raw
        df = df.rename({
            "ALT": "ALT_raw", 
            "ALT_resolved": "ALT",
            "AD": "AD_raw"
        })
        
        # Drop helpers
        df = df.drop(["allele_indices_str", "lookup_index", "AD_list"])

        # 8. Reorder for readability
        # Put is_compound_het after GT
        # Put AD_ref, AD_alt after AD_raw
        cols = df.columns
        
        if "is_compound_het" in cols and "GT" in cols:
            cols.remove("is_compound_het")
            cols.insert(cols.index("GT") + 1, "is_compound_het")
            
        if "AD_ref" in cols and "AD_alt" in cols and "AD_raw" in cols:
            cols.remove("AD_ref")
            cols.remove("AD_alt")
            idx_ad = cols.index("AD_raw") + 1
            cols.insert(idx_ad, "AD_ref")
            cols.insert(idx_ad + 1, "AD_alt")

        df = df.select(cols)

        # 9. Save
        df.write_csv(output_path, separator='\t')
        print(f"Processed genotypes (Split AD). Saved to {output_path}", file=sys.stderr)

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Resolve VCF genotypes, flag compound hets, and split AD.")
    parser.add_argument("input", help="Path to bcftools query output")
    parser.add_argument("output", help="Path to output TSV")
    args = parser.parse_args()
    
    resolve_genotypes(args.input, args.output)

#!/usr/bin/env python3
"""Extract sparse carrier genotypes for normalized allele pointers from VCZ."""

import argparse
import csv
import re
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import zarr


POINTER = re.compile(r"^z([0-9]+)_a([0-9]+)$")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--zarr", required=True, type=Path)
    parser.add_argument("--alleles", required=True, type=Path)
    parser.add_argument("--carriers-output", required=True, type=Path)
    parser.add_argument("--summary-output", required=True, type=Path)
    args = parser.parse_args()

    if args.alleles.suffix == ".parquet":
        alleles = pq.read_table(args.alleles).to_pylist()
    else:
        with args.alleles.open(newline="") as handle:
            alleles = list(csv.DictReader(handle, delimiter="\t"))
    pointers = []
    for row in alleles:
        match = POINTER.match(row["record_id"])
        if match is None:
            raise ValueError("invalid normalized allele pointer: {}".format(row["record_id"]))
        pointers.append((int(match.group(1)), int(match.group(2))))
    if len(set(pointers)) != len(pointers):
        raise ValueError("normalized allele pointers are not unique")

    root = zarr.open_group(str(args.zarr), mode="r")
    genotypes = root["call_genotype"]
    masks = root["call_genotype_mask"]
    sample_ids = [str(value) for value in root["sample_id"][:]]
    gq_array = root.get("call_GQ")
    dp_array = root.get("call_DP")
    laa_array = root.get("call_LAA")
    lad_array = root.get("call_LAD")
    phased_array = root.get("call_genotype_phased")
    variant_filter = root.get("variant_filter")
    filter_ids = [str(value) for value in root["filter_id"][:]]
    chunk_size = genotypes.chunks[0]

    carrier_rows = []
    summaries = []
    rows_by_chunk = {}
    for row_number, (variant_index, alt_index) in enumerate(pointers):
        if variant_index >= genotypes.shape[0]:
            raise ValueError("variant index out of bounds: {}".format(variant_index))
        rows_by_chunk.setdefault(variant_index // chunk_size, []).append(
            (row_number, variant_index, alt_index)
        )

    for chunk_id in sorted(rows_by_chunk):
        chunk_start = chunk_id * chunk_size
        chunk_end = min(chunk_start + chunk_size, genotypes.shape[0])
        gt_chunk = genotypes[chunk_start:chunk_end, :, :]
        mask_chunk = masks[chunk_start:chunk_end, :, :]
        gq_chunk = gq_array[chunk_start:chunk_end, :] if gq_array is not None else None
        dp_chunk = dp_array[chunk_start:chunk_end, :] if dp_array is not None else None
        laa_chunk = laa_array[chunk_start:chunk_end, :, :] if laa_array is not None else None
        lad_chunk = lad_array[chunk_start:chunk_end, :, :] if lad_array is not None else None
        phased_chunk = (
            phased_array[chunk_start:chunk_end, :] if phased_array is not None else None
        )
        filter_chunk = (
            variant_filter[chunk_start:chunk_end, :] if variant_filter is not None else None
        )

        for row_number, variant_index, alt_index in rows_by_chunk[chunk_id]:
            local_index = variant_index - chunk_start
            gt = gt_chunk[local_index]
            called = (~mask_chunk[local_index]) & (gt >= 0)
            dosage = ((gt == alt_index) & called).sum(axis=1)
            carrier_indexes = np.flatnonzero(dosage > 0)
            allele = alleles[row_number]

            an = int(called.sum())
            ac = int(dosage.sum())
            summaries.append({
                "record_id": allele["record_id"],
                "Gene": allele.get("Gene", ""),
                "SYMBOL": allele.get("SYMBOL", ""),
                "lof_tier": allele.get("lof_tier", ""),
                "genotype_ac": ac,
                "genotype_an": an,
                "genotype_af": ac / an if an else None,
                "carrier_count": int(len(carrier_indexes)),
                "hom_alt_count": int(np.count_nonzero(dosage == 2)),
            })

            site_filter = "."
            if filter_chunk is not None:
                site_filter = ";".join(
                    filter_ids[index]
                    for index, present in enumerate(filter_chunk[local_index])
                    if present
                ) or "."

            for sample_index in carrier_indexes.tolist():
                genotype = gt[sample_index]
                genotype_called = called[sample_index]
                called_ploidy = int(genotype_called.sum())
                sample_alt_dosage = int(dosage[sample_index])
                separator = "|" if (
                    phased_chunk is not None and phased_chunk[local_index, sample_index]
                ) else "/"
                gt_text = separator.join(
                    str(int(value)) if is_called else "."
                    for value, is_called in zip(genotype, genotype_called)
                )
                ad_ref = None
                ad_alt = None
                if laa_chunk is not None and lad_chunk is not None:
                    local_alt_slots = np.flatnonzero(
                        laa_chunk[local_index, sample_index] == alt_index
                    )
                    if len(local_alt_slots) == 1:
                        ad_ref_value = int(lad_chunk[local_index, sample_index, 0])
                        ad_alt_value = int(
                            lad_chunk[local_index, sample_index, int(local_alt_slots[0]) + 1]
                        )
                        ad_ref = ad_ref_value if ad_ref_value >= 0 else None
                        ad_alt = ad_alt_value if ad_alt_value >= 0 else None
                ad_text = (
                    "{},{}".format(ad_ref, ad_alt)
                    if ad_ref is not None and ad_alt is not None else ""
                )
                carrier_rows.append({
                    "record_id": allele["record_id"],
                    "#CHROM": "chr" + allele.get("CHROM", "").removeprefix("chr"),
                    "CHROM": allele.get("CHROM", ""),
                    "POS": int(allele["POS"]),
                    "REF": allele["REF"],
                    "ALT": allele["ALT"],
                    "Gene": allele.get("Gene", ""),
                    "SYMBOL": allele.get("SYMBOL", ""),
                    "Consequence": allele.get("Consequence", ""),
                    "LoF": allele.get("LoF", ""),
                    "lof_tier": allele.get("lof_tier", ""),
                    "miss_tier": allele.get("miss_tier", ""),
                    "miss_n_flag": (
                        int(allele["miss_n_flag"])
                        if allele.get("miss_n_flag") not in (None, "") else None
                    ),
                    "candidate_genes": allele.get("candidate_genes", ""),
                    "ClinPred_rankscore": allele.get("ClinPred_rankscore"),
                    "AlphaMissense_rankscore": allele.get("AlphaMissense_rankscore"),
                    "popEVE_converted_rankscore": allele.get("popEVE_converted_rankscore"),
                    "MPC_rankscore": allele.get("MPC_rankscore"),
                    "genebayes_post_mean": (
                        float(allele["genebayes_post_mean"])
                        if allele.get("genebayes_post_mean") not in (None, "") else None
                    ),
                    "variant_index": variant_index,
                    "alt_index": alt_index,
                    "sample_index": sample_index,
                    "sample_id": sample_ids[sample_index],
                    "FILTER": site_filter,
                    "GT": gt_text,
                    "AD": ad_text,
                    "genotype": gt_text,
                    # QC must use these mask-derived values, never the spelling
                    # of GT. Fixed-width stores may render a haploid call as
                    # `1/.`, while another source may render the same call as
                    # `1` or `./1`.
                    "called_ploidy": called_ploidy,
                    "alt_dosage": sample_alt_dosage,
                    "is_heterozygous": (
                        called_ploidy == 2 and sample_alt_dosage == 1
                    ),
                    "is_homozygous_alt": (
                        called_ploidy == 2 and sample_alt_dosage == 2
                    ),
                    "is_hemizygous_alt": (
                        called_ploidy == 1 and sample_alt_dosage == 1
                    ),
                    "GQ": int(gq_chunk[local_index, sample_index]) if gq_chunk is not None else None,
                    "DP": int(dp_chunk[local_index, sample_index]) if dp_chunk is not None else None,
                })

    args.carriers_output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pylist(carrier_rows), args.carriers_output, compression="zstd")
    pq.write_table(pa.Table.from_pylist(summaries), args.summary_output, compression="zstd")
    print("alleles={:,}".format(len(alleles)))
    print("carriers={:,}".format(len(carrier_rows)))
    print("carrier_samples={:,}".format(len({row["sample_id"] for row in carrier_rows})))


if __name__ == "__main__":
    main()

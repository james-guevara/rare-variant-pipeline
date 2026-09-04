#!/usr/bin/env bash
# Run one manifest-resolved targeted chromosome workflow in any environment.
# Every stage writes a named checkpoint and is skipped when that checkpoint is
# already nonempty. Override paths with the corresponding environment variable.
set -euo pipefail
export PATH="${HOME}/.local/bin:${PATH}"

repo=${RVP_REPO:-/opt/rvp}
: "${RUN_ROOT:?RUN_ROOT is required; use run_targeted_manifest.py}"
: "${ZARR_STORE:?ZARR_STORE is required; use run_targeted_manifest.py}"
: "${ANNOTATION_ROOT:?ANNOTATION_ROOT is required; use run_targeted_manifest.py}"
: "${LOFTEE_ROOT:?LOFTEE_ROOT is required; use run_targeted_manifest.py}"
: "${GENEBAYES:?GENEBAYES is required; use run_targeted_manifest.py}"
run_root=$RUN_ROOT
zarr_store=$ZARR_STORE
target_bed=${TARGET_BED:-}
if test "${ALL_OBSERVED:-0}" = 1; then target_bed=""; fi
fastvep_root=${FASTVEP_ROOT:-/opt/fastvep}
annotation_root=$ANNOTATION_ROOT
loftee_root=$LOFTEE_ROOT
genebayes=$GENEBAYES
missense_candidates=${MISSENSE_CANDIDATES:-}
missense_dbnsfp=${MISSENSE_DBNSFP:-}
postprocess_config=${POSTPROCESS_CONFIG:-}
sample_sex_qc=${SAMPLE_SEX_QC:-}
sex_chromosome_regions=${SEX_CHROMOSOME_REGIONS:-}
population_af_max=${POPULATION_AF_MAX:-0.01}
cohort_af_max=${COHORT_AF_MAX:-0.01}
lof_t1_min=${LOF_T1_MIN_GENEBAYES_POST_MEAN:-0.18}
lof_t2_min=${LOF_T2_MIN_GENEBAYES_POST_MEAN:-0.03}
missense_clinpred_min=${MISSENSE_CLINPRED_RANKSCORE_MIN:-0.4298}
missense_alphamissense_min=${MISSENSE_ALPHAMISSENSE_RANKSCORE_MIN:-0.9603}
missense_popeve_min=${MISSENSE_POPEVE_CONVERTED_RANKSCORE_MIN:-0.9209}
missense_mpc_min=${MISSENSE_MPC_RANKSCORE_MIN:-0.8947}
miss_t1_pass_count=${MISS_T1_PASS_COUNT:-4}
miss_t2_pass_count=${MISS_T2_PASS_COUNT:-3}
miss_t3_pass_count=${MISS_T3_PASS_COUNT:-2}
miss_t4_pass_count=${MISS_T4_PASS_COUNT:-1}
synonymous_controls=${SYNONYMOUS_TIERED_CONTROLS:-0}
family_genotypes=${FAMILY_GENOTYPES:-0}
sample_manifest=${SAMPLE_MANIFEST:-}
cohort=${COHORT:-g2mh}
chromosome=${CHROMOSOME:-chr22}
contig=${CONTIG:-${chromosome#chr}}
contig_length=${CONTIG_LENGTH:-50818468}

reference=$annotation_root/Homo_sapiens.GRCh38.dna.chromosome.$contig.fa
ancestor=$loftee_root/loftee-grch38/human_ancestor.fa.gz
gerp=$loftee_root/loftee-grch38/gerp_conservation_scores.homo_sapiens.GRCh38.bw
conservation=$loftee_root/loftee-grch38/loftee.sql
transcripts=$loftee_root/ensembl-115/$chromosome.transcripts.sqlite
gff3=$annotation_root/Homo_sapiens.GRCh38.115.$chromosome.gff3
transcript_cache=$gff3.fastvep.cache
transcript_priority=$annotation_root/vep115.$chromosome.transcript-priority.tsv
consequence_ranks=$annotation_root/vep115.consequence-ranks.tsv
fastvep=$fastvep_root/fastvep
fastvep_picker=${FASTVEP_PICKER:-fastvep-picker}
sex_chromosome=0
if test "$chromosome" = chrX || test "$chromosome" = chrY; then
  sex_chromosome=1
  test -r "$sample_sex_qc" || {
    echo "ERROR: sex chromosome requires SAMPLE_SEX_QC" >&2
    exit 1
  }
  test -r "$sex_chromosome_regions" || {
    echo "ERROR: sex chromosome requires SEX_CHROMOSOME_REGIONS" >&2
    exit 1
  }
fi

mkdir -p "$run_root"
if ! test -w "$run_root"; then
  echo "ERROR: run root is not writable: $run_root" >&2
  exit 1
fi
write_probe="$run_root/.write-test.$$"
if ! touch "$write_probe" || ! rm -f "$write_probe"; then
  echo "ERROR: run root failed write/delete preflight: $run_root" >&2
  exit 1
fi
if test -f "$run_root/_SUCCESS"; then
  printf 'workflow\tALREADY_SUCCEEDED\trun_root=%s\n' "$run_root"
  exit 0
fi
exec > >(tee -a "$run_root/workflow.log") 2>&1

required=(
  "$zarr_store/zarr.json" "$fastvep" "$reference" "$reference.fai"
  "$ancestor" "$ancestor.fai" "$ancestor.gzi" "$gerp" "$conservation"
  "$transcripts" "$gff3" "$transcript_cache" "$transcript_priority"
  "$consequence_ranks" "$genebayes"
)
if test -n "$target_bed"; then required+=("$target_bed"); fi
if test -n "$missense_dbnsfp"; then required+=("$missense_dbnsfp"); fi
for path in "${required[@]}"; do
  test -r "$path" || { echo "ERROR: missing input: $path" >&2; exit 1; }
done
command -v bcftools >/dev/null

cd "$repo"
python=${RVP_PYTHON:-$repo/.venv/bin/python}
if ! test -x "$python"; then
  command -v uv >/dev/null
  uv sync --frozen
fi
test -x "$python"

run_stage() {
  local output=$1
  shift
  if test -s "$output"; then
    printf 'checkpoint\tSKIP\t%s\n' "$output"
    return
  fi
  local start=$SECONDS
  rm -f "$output"
  if ! "$@"; then
    rm -f "$output"
    return 1
  fi
  if ! test -s "$output"; then
    rm -f "$output"
    echo "ERROR: stage completed without a nonempty output: $output" >&2
    return 1
  fi
  printf 'checkpoint\tDONE\t%s\telapsed_seconds=%s\n' "$output" "$((SECONDS-start))"
}

alleles=$run_root/01.target-alleles.parquet
raw_vcf=$run_root/02.target-sites.raw.vcf
normalized_vcf=$run_root/03.target-sites.normalized.vcf
picked=$run_root/04.fastvep-picked.tsv
loftee=$run_root/05.loftee.tsv
plof_all=$run_root/06.plof-genebayes.tsv
plof_tiered=$run_root/06.plof-tiered.tsv
carriers=$run_root/07.plof-tiered.carriers.parquet
carriers_raw=$run_root/07.plof-tiered.carriers.raw.parquet
summary=$run_root/07.plof-tiered.genotype-summary.parquet
plof_regions=$run_root/08.plof-region-filtered.parquet
plof_qc=$run_root/09.plof-genotype-qc.parquet
plof_pop_annotated=$run_root/10.plof-population-af-annotated.parquet
plof_pop_eligible=$run_root/10.plof-population-af-eligible.parquet
plof_cohort_annotated=$run_root/11.plof-cohort-af-annotated.parquet
plof_burden_eligible=$run_root/11.plof-burden-eligible.parquet
plof_counts=$run_root/12.plof-per-sample-counts.tsv
plof_totals=$run_root/12.plof-tier-totals.tsv
plof_sample_gene=$run_root/12.plof-primary-sample-gene.tsv
plof_sample_burden=$run_root/12.plof-primary-sample-burden.tsv
plof_sensitivity_gene=$run_root/12.plof-sensitivity-sample-gene.tsv
plof_sensitivity_burden=$run_root/12.plof-sensitivity-sample-burden.tsv
missense_tiered=$run_root/06.missense-tiered.parquet
generated_missense_candidates=$run_root/01.observed-missense-candidates.parquet
missense_carriers=$run_root/07.missense-tiered.carriers.parquet
missense_carriers_raw=$run_root/07.missense-tiered.carriers.raw.parquet
missense_summary=$run_root/07.missense-tiered.genotype-summary.parquet
missense_family=$run_root/07.missense-tiered.family-genotypes.parquet
missense_regions=$run_root/08.missense-region-filtered.parquet
missense_qc=$run_root/09.missense-genotype-qc.parquet
missense_pop_annotated=$run_root/10.missense-population-af-annotated.parquet
missense_pop_eligible=$run_root/10.missense-population-af-eligible.parquet
missense_cohort_annotated=$run_root/11.missense-cohort-af-annotated.parquet
missense_burden_eligible=$run_root/11.missense-burden-eligible.parquet
missense_counts=$run_root/12.missense-per-sample-counts.tsv
missense_totals=$run_root/12.missense-tier-totals.tsv
synonymous_tiered=$run_root/06.synonymous-tiered.tsv
synonymous_carriers=$run_root/07.synonymous-tiered.carriers.parquet
synonymous_carriers_raw=$run_root/07.synonymous-tiered.carriers.raw.parquet
synonymous_summary=$run_root/07.synonymous-tiered.genotype-summary.parquet
synonymous_family=$run_root/07.synonymous-tiered.family-genotypes.parquet
synonymous_regions=$run_root/08.synonymous-region-filtered.parquet
synonymous_qc=$run_root/09.synonymous-genotype-qc.parquet
synonymous_pop_annotated=$run_root/10.synonymous-population-af-annotated.parquet
synonymous_pop_eligible=$run_root/10.synonymous-population-af-eligible.parquet
synonymous_cohort_annotated=$run_root/11.synonymous-cohort-af-annotated.parquet
synonymous_burden_eligible=$run_root/11.synonymous-burden-eligible.parquet
synonymous_counts=$run_root/12.synonymous-per-sample-counts.tsv
synonymous_totals=$run_root/12.synonymous-tier-totals.tsv
synonymous_sample_gene=$run_root/12.synonymous-sample-gene.tsv
synonymous_sample_burden=$run_root/12.synonymous-sample-burden.tsv
plof_family=$run_root/07.plof-tiered.family-genotypes.parquet

if test "$family_genotypes" = 1; then
  test -r "$sample_manifest" || {
    echo "ERROR: FAMILY_GENOTYPES requires readable SAMPLE_MANIFEST" >&2
    exit 1
  }
fi

target_args=()
if test -n "$target_bed"; then target_args=(--bed "$target_bed"); fi
run_stage "$alleles" "$python" scripts/extract_zarr_target_alleles.py \
  --zarr "$zarr_store" "${target_args[@]}" --chrom "$chromosome" --output "$alleles"
run_stage "$raw_vcf" "$python" scripts/allele_parquet_to_sites_vcf.py \
  --input "$alleles" --output "$raw_vcf" --output-contig "$contig" \
  --contig-length "$contig_length"
run_stage "$normalized_vcf" bcftools norm -f "$reference" -m -any -Ov \
  -o "$normalized_vcf" "$raw_vcf"

if test -z "$missense_candidates" && test -n "$missense_dbnsfp"; then
  if ! test -s "$generated_missense_candidates"; then
    rm -f "$generated_missense_candidates"
    if ! "$python" scripts/build_missense_candidate_alleles.py \
        --dbnsfp "$missense_dbnsfp" --all-genes --chrom "$chromosome" \
        --clinpred-rankscore-min "$missense_clinpred_min" \
        --alphamissense-rankscore-min "$missense_alphamissense_min" \
        --popeve-converted-rankscore-min "$missense_popeve_min" \
        --mpc-rankscore-min "$missense_mpc_min" \
        --miss-t1-pass-count "$miss_t1_pass_count" \
        --miss-t2-pass-count "$miss_t2_pass_count" \
        --miss-t3-pass-count "$miss_t3_pass_count" \
        --miss-t4-pass-count "$miss_t4_pass_count" \
        --observed-vcf "$normalized_vcf" \
        --observed-output "$generated_missense_candidates"; then
      rm -f "$generated_missense_candidates"
      exit 1
    fi
  fi
  test -s "$generated_missense_candidates"
  missense_candidates=$generated_missense_candidates
fi

if ! test -s "$picked"; then
  start=$SECONDS
  rm -f "$picked"
  if ! "$fastvep" annotate --input "$normalized_vcf" --gff3 "$gff3" \
      --fasta "$reference" --transcript-cache "$transcript_cache" --hgvs \
      --symbol --canonical --output-format vcf --output - \
    | "$fastvep_picker" --fastvep - \
        --transcript-priority "$transcript_priority" \
        --consequence-ranks "$consequence_ranks" --output "$picked"; then
    rm -f "$picked"
    exit 1
  fi
  if ! test -s "$picked"; then
    rm -f "$picked"
    echo "ERROR: FastVEP/picker completed without a nonempty output: $picked" >&2
    exit 1
  fi
  printf 'checkpoint\tDONE\t%s\telapsed_seconds=%s\n' "$picked" "$((SECONDS-start))"
else
  printf 'checkpoint\tSKIP\t%s\n' "$picked"
fi

run_stage "$loftee" "$python" scripts/run_standalone_loftee.py \
  --input "$picked" --transcripts "$transcripts" --reference "$reference" \
  --ancestor "$ancestor" --gerp "$gerp" --conservation "$conservation" \
  --output "$loftee"

if ! test -s "$plof_all" || ! test -s "$plof_tiered"; then
  "$python" scripts/join_genebayes_lof_tiers.py --input "$loftee" \
    --genebayes "$genebayes" --output "$plof_all" \
    --qualifying-output "$plof_tiered" \
    --lof-t1-min "$lof_t1_min" --lof-t2-min "$lof_t2_min"
fi
test -s "$plof_all"
test -s "$plof_tiered"

extract_carriers=$carriers
if test "$sex_chromosome" = 1; then extract_carriers=$carriers_raw; fi
if ! test -s "$extract_carriers" || ! test -s "$summary" || \
   { test "$family_genotypes" = 1 && ! test -s "$plof_family"; }; then
  sex_extraction_args=()
  if test "$sex_chromosome" = 1; then
    sex_extraction_args=(--sample-sex-qc "$sample_sex_qc" \
      --sex-chromosome-regions "$sex_chromosome_regions")
  fi
  family_args=()
  if test "$family_genotypes" = 1; then
    family_args=(--sample-manifest "$sample_manifest" --family-output "$plof_family")
  fi
  "$python" scripts/extract_zarr_allele_genotypes.py --zarr "$zarr_store" \
    --alleles "$plof_tiered" --carriers-output "$extract_carriers" \
    --summary-output "$summary" "${sex_extraction_args[@]}" \
    "${family_args[@]}"
fi
if test "$sex_chromosome" = 1; then
  run_stage "$carriers" "$python" scripts/annotate_sex_chromosome_carriers.py \
    --input "$carriers_raw" --sample-qc "$sample_sex_qc" \
    --regions "$sex_chromosome_regions" --output "$carriers"
fi
test -s "$carriers"
test -s "$summary"

if test -n "$postprocess_config"; then
  test -r "$postprocess_config" || {
    echo "ERROR: missing postprocess config: $postprocess_config" >&2
    exit 1
  }
  run_stage "$plof_regions" "$python" scripts/postprocess/filter_regions.py \
    --cohort "$cohort" --chrom "$chromosome" --resources "$postprocess_config" \
    --input "$carriers" --output "$plof_regions"
  run_stage "$plof_qc" "$python" scripts/postprocess/qc_genotype.py \
    --cohort "$cohort" --chrom "$chromosome" --resources "$postprocess_config" \
    --input "$plof_regions" --output "$plof_qc"
  lof_af_enabled=$(
    "$python" -c \
      'import json,sys; print(int(bool(json.load(open(sys.argv[1])).get("dbnsfp_af_dir"))))' \
      "$postprocess_config"
  )
  if test "$lof_af_enabled" = 1; then
    run_stage "$plof_pop_annotated" "$python" scripts/postprocess/join_pop_af.py \
      --cohort "$cohort" --chrom "$chromosome" --resources "$postprocess_config" \
      --input "$plof_qc" --output "$plof_pop_annotated"
    run_stage "$plof_pop_eligible" "$python" scripts/apply_population_af_filter.py \
      --input "$plof_pop_annotated" --output "$plof_pop_eligible" \
      --column gnomAD4.1_joint_AF --max-af "$population_af_max"
    if ! test -s "$plof_cohort_annotated" || ! test -s "$plof_burden_eligible"; then
      summary_prefix=genotype
      if test "$sex_chromosome" = 1; then summary_prefix=primary_genotype; fi
      "$python" scripts/apply_cohort_af_filter.py \
        --input "$plof_pop_eligible" --allele-summary "$summary" \
        --max-af "$cohort_af_max" --summary-prefix "$summary_prefix" \
        --annotated-output "$plof_cohort_annotated" \
        --eligible-output "$plof_burden_eligible"
    fi
    test -s "$plof_cohort_annotated"
    test -s "$plof_burden_eligible"
    if ! test -s "$plof_counts" || ! test -s "$plof_totals"; then
      count_eligibility=()
      if test "$sex_chromosome" = 1; then
        count_eligibility=(--eligibility-col primary_analysis_eligible)
      fi
      "$python" scripts/postprocess/count_carriers.py \
        --input "$plof_burden_eligible" --sample-col sample_id \
        --group-col lof_tier --out-counts "$plof_counts" \
        --out-totals "$plof_totals" "${count_eligibility[@]}"
    fi
    test -s "$plof_counts"
    test -s "$plof_totals"
    if ! test -s "$plof_sample_gene" || ! test -s "$plof_sample_burden"; then
      collapse_args=()
      if test "$sex_chromosome" = 1; then
        collapse_args=(
          --eligibility-col primary_analysis_eligible
          --burden-available-col burden_count_available
          --sensitivity-sample-gene-output "$plof_sensitivity_gene"
          --sensitivity-sample-burden-output "$plof_sensitivity_burden"
        )
      fi
      "$python" scripts/collapse_lof_carriers.py \
        --input "$plof_burden_eligible" \
        --sample-gene-output "$plof_sample_gene" \
        --sample-burden-output "$plof_sample_burden" \
        "${collapse_args[@]}"
    fi
    test -s "$plof_sample_gene"
    test -s "$plof_sample_burden"
    if test "$sex_chromosome" = 1; then
      test -s "$plof_sensitivity_gene"
      test -s "$plof_sensitivity_burden"
    fi
  fi
fi

if test "$synonymous_controls" = 1; then
  run_stage "$synonymous_tiered" "$python" scripts/select_synonymous_tiered.py \
    --picked "$picked" --genebayes "$genebayes" --output "$synonymous_tiered"
  if awk 'NR > 1 { found=1; exit } END { exit !found }' "$synonymous_tiered"; then
  extract_synonymous_carriers=$synonymous_carriers
  if test "$sex_chromosome" = 1; then extract_synonymous_carriers=$synonymous_carriers_raw; fi
  if ! test -s "$extract_synonymous_carriers" || ! test -s "$synonymous_summary" || \
     { test "$family_genotypes" = 1 && ! test -s "$synonymous_family"; }; then
    sex_extraction_args=()
    if test "$sex_chromosome" = 1; then
      sex_extraction_args=(--sample-sex-qc "$sample_sex_qc" \
        --sex-chromosome-regions "$sex_chromosome_regions")
    fi
    family_args=()
    if test "$family_genotypes" = 1; then
      family_args=(--sample-manifest "$sample_manifest" --family-output "$synonymous_family")
    fi
    "$python" scripts/extract_zarr_allele_genotypes.py --zarr "$zarr_store" \
      --alleles "$synonymous_tiered" --carriers-output "$extract_synonymous_carriers" \
      --summary-output "$synonymous_summary" "${sex_extraction_args[@]}" \
      "${family_args[@]}"
  fi
  if test "$sex_chromosome" = 1; then
    run_stage "$synonymous_carriers" "$python" scripts/annotate_sex_chromosome_carriers.py \
      --input "$synonymous_carriers_raw" --sample-qc "$sample_sex_qc" \
      --regions "$sex_chromosome_regions" --output "$synonymous_carriers"
  fi
  test -s "$synonymous_carriers"
  test -s "$synonymous_summary"

  if test -n "$postprocess_config"; then
    run_stage "$synonymous_regions" "$python" scripts/postprocess/filter_regions.py \
      --cohort "$cohort" --chrom "$chromosome" --resources "$postprocess_config" \
      --input "$synonymous_carriers" --output "$synonymous_regions"
    run_stage "$synonymous_qc" "$python" scripts/postprocess/qc_genotype.py \
      --cohort "$cohort" --chrom "$chromosome" --resources "$postprocess_config" \
      --input "$synonymous_regions" --output "$synonymous_qc"
    run_stage "$synonymous_pop_annotated" "$python" scripts/postprocess/join_pop_af.py \
      --cohort "$cohort" --chrom "$chromosome" --resources "$postprocess_config" \
      --input "$synonymous_qc" --output "$synonymous_pop_annotated"
    run_stage "$synonymous_pop_eligible" "$python" scripts/apply_population_af_filter.py \
      --input "$synonymous_pop_annotated" --output "$synonymous_pop_eligible" \
      --column gnomAD4.1_joint_AF --max-af "$population_af_max"
    if ! test -s "$synonymous_cohort_annotated" || ! test -s "$synonymous_burden_eligible"; then
      summary_prefix=genotype
      if test "$sex_chromosome" = 1; then summary_prefix=primary_genotype; fi
      "$python" scripts/apply_cohort_af_filter.py \
        --input "$synonymous_pop_eligible" --allele-summary "$synonymous_summary" \
        --max-af "$cohort_af_max" --summary-prefix "$summary_prefix" \
        --annotated-output "$synonymous_cohort_annotated" \
        --eligible-output "$synonymous_burden_eligible"
    fi
    count_eligibility=()
    collapse_args=()
    if test "$sex_chromosome" = 1; then
      count_eligibility=(--eligibility-col primary_analysis_eligible)
      collapse_args=(--eligibility-col primary_analysis_eligible \
        --burden-available-col burden_count_available)
    fi
    "$python" scripts/postprocess/count_carriers.py \
      --input "$synonymous_burden_eligible" --sample-col sample_id \
      --group-col lof_tier --out-counts "$synonymous_counts" \
      --out-totals "$synonymous_totals" "${count_eligibility[@]}"
    "$python" scripts/collapse_lof_carriers.py \
      --input "$synonymous_burden_eligible" \
      --sample-gene-output "$synonymous_sample_gene" \
      --sample-burden-output "$synonymous_sample_burden" "${collapse_args[@]}"
  fi
  else
    printf 'synonymous\tNO_TIERED_ALLELES\tchromosome=%s\n' "$chromosome"
  fi
fi

if test -n "$missense_candidates"; then
  test -r "$missense_candidates" || {
    echo "ERROR: missing missense candidates: $missense_candidates" >&2
    exit 1
  }
  run_stage "$missense_tiered" "$python" scripts/select_missense_candidates.py \
    --picked "$picked" --candidates "$missense_candidates" \
    --output "$missense_tiered"
  extract_missense_carriers=$missense_carriers
  if test "$sex_chromosome" = 1; then
    extract_missense_carriers=$missense_carriers_raw
  fi
  if ! test -s "$extract_missense_carriers" || ! test -s "$missense_summary" || \
     { test "$family_genotypes" = 1 && ! test -s "$missense_family"; }; then
    sex_extraction_args=()
    if test "$sex_chromosome" = 1; then
      sex_extraction_args=(--sample-sex-qc "$sample_sex_qc" \
        --sex-chromosome-regions "$sex_chromosome_regions")
    fi
    family_args=()
    if test "$family_genotypes" = 1; then
      family_args=(--sample-manifest "$sample_manifest" --family-output "$missense_family")
    fi
    "$python" scripts/extract_zarr_allele_genotypes.py --zarr "$zarr_store" \
      --alleles "$missense_tiered" --carriers-output "$extract_missense_carriers" \
      --summary-output "$missense_summary" "${sex_extraction_args[@]}" \
      "${family_args[@]}"
  fi
  if test "$sex_chromosome" = 1; then
    run_stage "$missense_carriers" "$python" scripts/annotate_sex_chromosome_carriers.py \
      --input "$missense_carriers_raw" --sample-qc "$sample_sex_qc" \
      --regions "$sex_chromosome_regions" --output "$missense_carriers"
  fi
  test -s "$missense_carriers"
  test -s "$missense_summary"

  if test -n "$postprocess_config"; then
    test -r "$postprocess_config" || {
      echo "ERROR: missing postprocess config: $postprocess_config" >&2
      exit 1
    }
    run_stage "$missense_regions" "$python" scripts/postprocess/filter_regions.py \
      --cohort "$cohort" --chrom "$chromosome" --resources "$postprocess_config" \
      --input "$missense_carriers" --output "$missense_regions"
    run_stage "$missense_qc" "$python" scripts/postprocess/qc_genotype.py \
      --cohort "$cohort" --chrom "$chromosome" --resources "$postprocess_config" \
      --input "$missense_regions" --output "$missense_qc"
    run_stage "$missense_pop_annotated" "$python" scripts/postprocess/join_pop_af.py \
      --cohort "$cohort" --chrom "$chromosome" --resources "$postprocess_config" \
      --input "$missense_qc" --output "$missense_pop_annotated"
    run_stage "$missense_pop_eligible" "$python" scripts/apply_population_af_filter.py \
      --input "$missense_pop_annotated" --output "$missense_pop_eligible" \
      --column gnomAD4.1_joint_AF --max-af "$population_af_max"
    if ! test -s "$missense_cohort_annotated" || ! test -s "$missense_burden_eligible"; then
      summary_prefix=genotype
      if test "$sex_chromosome" = 1; then summary_prefix=primary_genotype; fi
      "$python" scripts/apply_cohort_af_filter.py \
        --input "$missense_pop_eligible" --allele-summary "$missense_summary" \
        --max-af "$cohort_af_max" --summary-prefix "$summary_prefix" \
        --annotated-output "$missense_cohort_annotated" \
        --eligible-output "$missense_burden_eligible"
    fi
    test -s "$missense_cohort_annotated"
    test -s "$missense_burden_eligible"
    if ! test -s "$missense_counts" || ! test -s "$missense_totals"; then
      count_eligibility=()
      if test "$sex_chromosome" = 1; then
        count_eligibility=(--eligibility-col primary_analysis_eligible)
      fi
      "$python" scripts/postprocess/count_carriers.py \
        --input "$missense_burden_eligible" --sample-col sample_id \
        --group-col miss_tier --out-counts "$missense_counts" \
        --out-totals "$missense_totals" "${count_eligibility[@]}"
    fi
    test -s "$missense_counts"
    test -s "$missense_totals"
  fi
fi

touch "$run_root/_SUCCESS"
printf 'workflow\tSUCCEEDED\trun_root=%s\n' "$run_root"

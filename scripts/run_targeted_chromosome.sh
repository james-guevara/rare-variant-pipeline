#!/usr/bin/env bash
# Run one manifest-resolved targeted chromosome workflow in any environment.
# Every stage writes a named checkpoint and is skipped when that checkpoint is
# already nonempty. Override paths with the corresponding environment variable.
set -euo pipefail
export PATH="${HOME}/.local/bin:${PATH}"

repo=${RVP_REPO:-/opt/rvp}
: "${RUN_ROOT:?RUN_ROOT is required; use run_targeted_manifest.py}"
: "${ZARR_STORE:?ZARR_STORE is required; use run_targeted_manifest.py}"
: "${TARGET_BED:?TARGET_BED is required; use run_targeted_manifest.py}"
: "${ANNOTATION_ROOT:?ANNOTATION_ROOT is required; use run_targeted_manifest.py}"
: "${LOFTEE_ROOT:?LOFTEE_ROOT is required; use run_targeted_manifest.py}"
: "${GENEBAYES:?GENEBAYES is required; use run_targeted_manifest.py}"
run_root=$RUN_ROOT
zarr_store=$ZARR_STORE
target_bed=$TARGET_BED
fastvep_root=${FASTVEP_ROOT:-/opt/fastvep}
annotation_root=$ANNOTATION_ROOT
loftee_root=$LOFTEE_ROOT
genebayes=$GENEBAYES
missense_candidates=${MISSENSE_CANDIDATES:-}
postprocess_config=${POSTPROCESS_CONFIG:-}
sample_sex_qc=${SAMPLE_SEX_QC:-}
sex_chromosome_regions=${SEX_CHROMOSOME_REGIONS:-}
population_af_max=${POPULATION_AF_MAX:-0.01}
cohort_af_max=${COHORT_AF_MAX:-0.01}
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
  "$zarr_store/zarr.json" "$target_bed" "$fastvep" "$reference" "$reference.fai"
  "$ancestor" "$ancestor.fai" "$ancestor.gzi" "$gerp" "$conservation"
  "$transcripts" "$gff3" "$transcript_cache" "$transcript_priority"
  "$consequence_ranks" "$genebayes"
)
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
  "$@"
  test -s "$output"
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
missense_tiered=$run_root/06.missense-tiered.parquet
missense_carriers=$run_root/07.missense-tiered.carriers.parquet
missense_carriers_raw=$run_root/07.missense-tiered.carriers.raw.parquet
missense_summary=$run_root/07.missense-tiered.genotype-summary.parquet
missense_regions=$run_root/08.missense-region-filtered.parquet
missense_qc=$run_root/09.missense-genotype-qc.parquet
missense_pop_annotated=$run_root/10.missense-population-af-annotated.parquet
missense_pop_eligible=$run_root/10.missense-population-af-eligible.parquet
missense_cohort_annotated=$run_root/11.missense-cohort-af-annotated.parquet
missense_burden_eligible=$run_root/11.missense-burden-eligible.parquet
missense_counts=$run_root/12.missense-per-sample-counts.tsv
missense_totals=$run_root/12.missense-tier-totals.tsv

run_stage "$alleles" "$python" scripts/extract_zarr_target_alleles.py \
  --zarr "$zarr_store" --bed "$target_bed" --chrom "$chromosome" --output "$alleles"
run_stage "$raw_vcf" "$python" scripts/allele_parquet_to_sites_vcf.py \
  --input "$alleles" --output "$raw_vcf" --output-contig "$contig" \
  --contig-length "$contig_length"
run_stage "$normalized_vcf" bcftools norm -f "$reference" -m -any -Ov \
  -o "$normalized_vcf" "$raw_vcf"

if ! test -s "$picked"; then
  start=$SECONDS
  "$fastvep" annotate --input "$normalized_vcf" --gff3 "$gff3" \
    --fasta "$reference" --transcript-cache "$transcript_cache" --hgvs \
    --symbol --canonical --output-format vcf --output - \
  | "$python" scripts/pick_fastvep_consequences.py --fastvep - \
      --transcript-priority "$transcript_priority" \
      --consequence-ranks "$consequence_ranks" --output "$picked"
  test -s "$picked"
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
    --qualifying-output "$plof_tiered"
fi
test -s "$plof_all"
test -s "$plof_tiered"

extract_carriers=$carriers
if test "$sex_chromosome" = 1; then extract_carriers=$carriers_raw; fi
if ! test -s "$extract_carriers" || ! test -s "$summary"; then
  "$python" scripts/extract_zarr_allele_genotypes.py --zarr "$zarr_store" \
    --alleles "$plof_tiered" --carriers-output "$extract_carriers" \
    --summary-output "$summary"
fi
if test "$sex_chromosome" = 1; then
  run_stage "$carriers" "$python" scripts/annotate_sex_chromosome_carriers.py \
    --input "$carriers_raw" --sample-qc "$sample_sex_qc" \
    --regions "$sex_chromosome_regions" --output "$carriers"
fi
test -s "$carriers"
test -s "$summary"

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
  if ! test -s "$extract_missense_carriers" || ! test -s "$missense_summary"; then
    "$python" scripts/extract_zarr_allele_genotypes.py --zarr "$zarr_store" \
      --alleles "$missense_tiered" --carriers-output "$extract_missense_carriers" \
      --summary-output "$missense_summary"
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
      "$python" scripts/apply_cohort_af_filter.py \
        --input "$missense_pop_eligible" --allele-summary "$missense_summary" \
        --max-af "$cohort_af_max" \
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

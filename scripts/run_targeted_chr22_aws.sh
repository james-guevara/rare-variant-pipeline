#!/usr/bin/env bash
# Run the integrated chr22 targeted LoF workflow on the persistent AWS host.
# Every stage writes a named checkpoint and is skipped when that checkpoint is
# already nonempty. Override paths with the corresponding environment variable.
set -euo pipefail
export PATH="${HOME}/.local/bin:${PATH}"

repo=${RVP_REPO:-/home/ubuntu/work/rare-variant-pipeline-targeted-aws}
run_root=${RUN_ROOT:-/fsx/rare-variant-pilot/targeted-workflows/g2mh/chr22-lof-v1}
zarr_store=${ZARR_STORE:-/fsx/rare-variant-pilot/g2mh-vcz-v3/v1/chr22.sharded-v3.zarr}
target_bed=${TARGET_BED:-/home/ubuntu/work/rare-variant-missense-scope/targets/lof-plus-missense-candidates.chr22.bed}
fastvep_root=${FASTVEP_ROOT:-/home/ubuntu/work/fastvep-runtime}
annotation_root=${ANNOTATION_ROOT:-$fastvep_root/ensembl-115}
loftee_root=${LOFTEE_ROOT:-/fsx/loftee-parity/resources}
genebayes=${GENEBAYES:-/home/ubuntu/work/standalone-loftee-test/port/GeneBayes.Supplementary_Table_1.tsv}
missense_candidates=${MISSENSE_CANDIDATES:-}
postprocess_config=${POSTPROCESS_CONFIG:-}
population_af_max=${POPULATION_AF_MAX:-0.01}
cohort_af_max=${COHORT_AF_MAX:-0.01}

reference=$annotation_root/Homo_sapiens.GRCh38.dna.chromosome.22.fa
ancestor=$loftee_root/loftee-grch38/human_ancestor.fa.gz
gerp=$loftee_root/loftee-grch38/gerp_conservation_scores.homo_sapiens.GRCh38.bw
conservation=$loftee_root/loftee-grch38/loftee.sql
transcripts=$loftee_root/ensembl-115/chr22.transcripts.sqlite
gff3=$annotation_root/Homo_sapiens.GRCh38.115.chr22.gff3
transcript_cache=$gff3.fastvep.cache
transcript_priority=$annotation_root/vep115.chr22.transcript-priority.tsv
consequence_ranks=$annotation_root/vep115.consequence-ranks.tsv
fastvep=$fastvep_root/fastvep

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
summary=$run_root/07.plof-tiered.genotype-summary.parquet
missense_tiered=$run_root/06.missense-tiered.parquet
missense_carriers=$run_root/07.missense-tiered.carriers.parquet
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
  --zarr "$zarr_store" --bed "$target_bed" --chrom chr22 --output "$alleles"
run_stage "$raw_vcf" "$python" scripts/allele_parquet_to_sites_vcf.py \
  --input "$alleles" --output "$raw_vcf" --output-contig 22 --contig-length 50818468
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

if ! test -s "$carriers" || ! test -s "$summary"; then
  "$python" scripts/extract_zarr_allele_genotypes.py --zarr "$zarr_store" \
    --alleles "$plof_tiered" --carriers-output "$carriers" \
    --summary-output "$summary"
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
  if ! test -s "$missense_carriers" || ! test -s "$missense_summary"; then
    "$python" scripts/extract_zarr_allele_genotypes.py --zarr "$zarr_store" \
      --alleles "$missense_tiered" --carriers-output "$missense_carriers" \
      --summary-output "$missense_summary"
  fi
  test -s "$missense_carriers"
  test -s "$missense_summary"

  if test -n "$postprocess_config"; then
    test -r "$postprocess_config" || {
      echo "ERROR: missing postprocess config: $postprocess_config" >&2
      exit 1
    }
    run_stage "$missense_regions" "$python" scripts/postprocess/filter_regions.py \
      --cohort g2mh --chrom chr22 --resources "$postprocess_config" \
      --input "$missense_carriers" --output "$missense_regions"
    run_stage "$missense_qc" "$python" scripts/postprocess/qc_genotype.py \
      --cohort g2mh --chrom chr22 --resources "$postprocess_config" \
      --input "$missense_regions" --output "$missense_qc"
    run_stage "$missense_pop_annotated" "$python" scripts/postprocess/join_pop_af.py \
      --cohort g2mh --chrom chr22 --resources "$postprocess_config" \
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
      "$python" scripts/postprocess/count_carriers.py \
        --input "$missense_burden_eligible" --sample-col sample_id \
        --group-col miss_tier --out-counts "$missense_counts" \
        --out-totals "$missense_totals"
    fi
    test -s "$missense_counts"
    test -s "$missense_totals"
  fi
fi

touch "$run_root/_SUCCESS"
printf 'workflow\tSUCCEEDED\trun_root=%s\n' "$run_root"

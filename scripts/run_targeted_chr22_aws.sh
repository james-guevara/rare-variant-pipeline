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
command -v uv >/dev/null
command -v bcftools >/dev/null

cd "$repo"
uv sync --frozen
python=$repo/.venv/bin/python

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

printf 'workflow\tSUCCEEDED\trun_root=%s\n' "$run_root"

#!/usr/bin/env bash
# Build one chromosome's Ensembl 115 FastVEP and standalone-LOFTEE resources.
set -euo pipefail

chromosome=${1:?usage: prepare_ensembl_chromosome_resources.sh CHROMOSOME CONTIG}
contig=${2:?usage: prepare_ensembl_chromosome_resources.sh CHROMOSOME CONTIG}
source_root=${ENSEMBL_SOURCE_ROOT:?ENSEMBL_SOURCE_ROOT is required}
annotation_root=${ANNOTATION_ROOT:?ANNOTATION_ROOT is required}
loftee_root=${LOFTEE_ROOT:?LOFTEE_ROOT is required}
gtf=${ENSEMBL_GTF:?ENSEMBL_GTF is required}
fastvep=${FASTVEP:?FASTVEP is required}
python=${RVP_PYTHON:?RVP_PYTHON is required}
repo=${RVP_REPO:?RVP_REPO is required}

release=115
gff3_source=$source_root/Homo_sapiens.GRCh38.$release.gff3.gz
fasta_gz=$source_root/Homo_sapiens.GRCh38.dna.chromosome.$contig.fa.gz
gff3=$annotation_root/Homo_sapiens.GRCh38.$release.$chromosome.gff3
fasta=$annotation_root/Homo_sapiens.GRCh38.dna.chromosome.$contig.fa
cache=$gff3.fastvep.cache
transcripts=$loftee_root/ensembl-$release/$chromosome.transcripts.sqlite

mkdir -p "$source_root" "$annotation_root" "$(dirname "$transcripts")"
if ! test -s "$gff3_source"; then
  curl --fail --location --retry 5 \
    --output "$gff3_source" \
    "https://ftp.ensembl.org/pub/release-$release/gff3/homo_sapiens/$(basename "$gff3_source")"
fi
if ! test -s "$fasta_gz"; then
  curl --fail --location --retry 5 \
    --output "$fasta_gz" \
    "https://ftp.ensembl.org/pub/release-$release/fasta/homo_sapiens/dna/$(basename "$fasta_gz")"
fi

if ! test -s "$gff3"; then
  gzip -dc "$gff3_source" \
    | awk -F '\t' -v contig="$contig" 'substr($0,1,1) == "#" || $1 == contig' \
    > "$gff3"
fi
if ! test -s "$fasta"; then
  gzip -dc "$fasta_gz" > "$fasta"
fi
if ! test -s "$fasta.fai"; then
  "$python" -c 'import pysam,sys; pysam.faidx(sys.argv[1])' "$fasta"
fi
if ! test -s "$cache"; then
  "$fastvep" cache --gff3 "$gff3" --fasta "$fasta" --output "$cache"
fi
if ! test -s "$transcripts"; then
  cd "$repo"
  "$python" scripts/build_loftee_transcript_db.py \
    --gtf "$gtf" --seqname "$contig" --output "$transcripts"
fi

printf 'resource_ready\t%s\tgff3=%s\tfasta=%s\tcache=%s\ttranscripts=%s\n' \
  "$chromosome" "$(stat -c %s "$gff3")" "$(stat -c %s "$fasta")" \
  "$(stat -c %s "$cache")" "$(stat -c %s "$transcripts")"

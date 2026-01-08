micromamba activate python3.12_env_default
python family_query.py --vcf vcfs/chr22_jointcall_VQSR_combined.vcf.gz --ped wgs.psam --out out/chr22.family_genotypes.tsv --region out/chr22.bed

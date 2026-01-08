module load singularitypro
BIND_FOLDER="/expanse/projects/sebat1/"
BCFTOOLS_CONTAINER="containers/bcftools:1.22--h3a4d415_1"

mkdir -p out

# Loop through 1..22, X, and Y
for i in {1..22} X Y; do
    CHR="chr${i}"
    echo "Processing ${CHR}..."
    singularity exec --bind "${BIND_FOLDER}" "${BCFTOOLS_CONTAINER}" \
		bcftools +split-vep -p CSQ -HH -A '\t' -d -s :missense+ -f '%CHROM\t%POS0\t%END\t%POS\t%REF\t%ALT\t%ID\t%QUAL\t%INFO\t%CSQ\n' out_vep/${CHR}.vep.vcf.gz | \
  		sed -E '2,$s/;?CSQ[^=]*=[^;\t]*//g' | \
  		sed -E '1s/CSQ//g; 1s/\[[0-9]+\]//g; 1s/\(null\)/INFO/g' | \
      	bgzip > out/${CHR}.variants.tsv.gz
done




#for i in {1..22} X Y; do
#    CHR="chr${i}"
#    echo "Converting ${CHR} variants file to bed file..."
#	bedtools merge -i out/${CHR}.variants.tsv.gz > out/${CHR}.bed 
#done
#
#
#for i in {1..22} X Y; do
#    CHR="chr${i}"
#    echo "Creating ${CHR} genotypes file..."
#	singularity exec --bind "${BIND_FOLDER}" "${BCFTOOLS_CONTAINER}" \
#		bcftools query -R out/${CHR}.bed -HH -i 'GT="alt"' -f '[%CHROM\t%POS0\t%END\t%POS\t%REF\t%ALT\t%SAMPLE\t%GT\t%GQ\t%DP\t%AD\n]' vcfs/${CHR}_jointcall_VQSR_combined.vcf.gz | bgzip > out/${CHR}.genotypes.tsv.gz 
#done
#
#

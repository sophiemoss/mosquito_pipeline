## To do Principal Components Analysis and make PCA plots, one way is to make an allelic matrix from your VCF and then use this in the subsequent R script.

## We're going to use plink here to make the matrix, and for that we have to rename the X and Mt chromosomes, otherwise plink removes these variants as it thinks that they are germline

## Example of renaming mitochondira within the vcf from Mt to anop_mito, and X to anop_X, otherwise plink removes the variants
vim chrom_map.txt 
Mt  anop_mito
X   anop_X

bcftools annotate --rename-chrs chrom_map.txt FMISSING_MAF_AC0_DP5_GQ20_gatk_filtered_miss_40_mac_bi_snps_melas_2019_plusglobal.2023_07_25.genotyped.ann.vcf.gz -Oz -o pca_filteredvcf_renamedchr_melas2019plusglobal.vcf.gz
tabix -p vcf pca_filteredvcf_renamedchr_melas2019plusglobal.vcf.gz

## Check that this has properly renamed the chromosomes. Look at unique_chromosomes.txt using 'cat unique_chromosomes.txt'
zgrep -v '^##' pca_filteredvcf_renamedchr_melas2019plusglobal.vcf.gz | cut -f1 | sort | uniq > unique_chromosomes.txt

## Create a distance matrix for PCA with melas plus global samples
plink --vcf final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.vcf.gz --distance square --double-id --out wholgenome_melas_plusglobal --threads 6 --allow-extra-chr

## You might want to also make matrices for separate chromosomes
## Here is an example of how to do that for just the variants on the mitochondrial (anop_mito) genome.
## create mitochondria only vcf using bcftools, then use plink to make distance matrix
bcftools view -r anop_mito final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.vcf.gz -Oz -o mito_only_final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.vcf.gz
tabix -p vcf mito_only_final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.vcf.gz
plink --vcf mito_only_final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.vcf.gz --distance square --double-id --out mito_only_melas_plusglobal --threads 6 --allow-extra-chr


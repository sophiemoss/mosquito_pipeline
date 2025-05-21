## Admixture

# Identify which samples are in your vcf:
bcftools query -l yourfile.vcf.gz

# STEP 1 EXTRA FILTERING for MAF if you need to do it

# STEP 2 CONVERT CHR NAMES IN VCF TO INTEGERS

zgrep -v "^#" PBO_exposed_F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_filt_coluzzii.2025_02_28.genotyped.ann.vcf.gz | cut -f1 | sort | uniq

zcat PBO_exposed_F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_filt_coluzzii.2025_02_28.genotyped.ann.vcf.gz | awk 'BEGIN{OFS=FS="\t"} /^#/ {print; next} {gsub(/^NC_064604.1$/, "1", $1); gsub(/^NC_064669.1$/, "2", $1); gsub(/^NC_064670.1$/, "3", $1); gsub(/^NC_064671.1$/, "4", $1); print}' | bgzip > admixture.vcf.gz

tabix -p vcf admixture.vcf.gz
zgrep -v "^#" admixture.vcf.gz | cut -f1 | sort | uniq
bcftools query -l admixture.vcf.gz

# STEP 3 MAKE BED AND BIM FILES

plink --vcf admixture.vcf.gz --set-missing-var-ids @:# --keep-allele-order --const-fid --allow-extra-chr --make-bed --out guinea_coluzzii_pbo_exposed

# STEP 4 RUN ADMIXTURE
# Make a file called K_runs.txt with values 1 to 10 and then run the following command
# Run this in screen as it takes a while

cat K_runs.txt | xargs -I {} sh -c 'admixture --cv=10 -j20 -s 14062 guinea_coluzzii_pbo_exposed.bed {} | tee log{}.cv10.seed82763.out'

# Inspect the CV values for the inflection point to check what the best K value is for your dataset
# Inspect using grep -h CV *out

# Plot admixture using R script
# Change the directory and metadata below as needed

conda create -n radmix r-essentials r-base
install.packages(c("unikn", "countrycode", "optparse"))
setwd("/mnt/storage11/sophie/admixture")

Rscript generate_admix_barplot_colours.R \
-d /mnt/storage11/sophie/admixture \
--prefix melas_global_gambiaealigned \
--kval 4 \
-m metadata_melasplusglobal_admixture.tsv \
--filter_N 1 \
--label_id sample \
--label_region region \
--label_country country \
--label_site site \
--country_code TRUE


### ____ Troubleshooting ___ ###
bcftools query -f '[%CHROM\t%POS\t%GT\n]' F_MISS_MAF_admixture_modified.vcf.gz | \
awk '{print $1, $2, $3}' | sort | uniq -c | grep -w '1' > invariant_sites.txt

bcftools query -f '[%CHROM\t%POS\t%GT\n]' PBO_exposed_F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_filt_coluzzii.2025_02_28.genotyped.ann.vcf.gz | \
awk '{if ($3 == "./.") print $0}' | sort | uniq -c | \ awk '{if ($1 == <number_of_samples>) print $2, $3}' > all_missing_sites.txt


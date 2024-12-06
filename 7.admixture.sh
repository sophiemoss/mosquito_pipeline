## Admixture

# Identify which samples are in your vcf:
bcftools query -l yourfile.vcf.gz

# STEP 1 EXTRA FILTERING for MAF if you need to do it
# MAF > 0.01 filter has already been applied to my VCF so I am not doing any more filtering

# STEP 2 CONVERT CHR NAMES IN VCF TO INTEGERS

zgrep -v "^#" final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.vcf.gz | cut -f1 | sort | uniq

zcat final_filteredvcf.vcf.gz | awk 'BEGIN{OFS=FS="\t"} /^#/ {print; next} {gsub(/^2L$/, "1", $1); gsub(/^2R$/, "2", $1); gsub(/^3L$/, "3", $1); gsub(/^3R$/, "4", $1); gsub(/^anop_mito$/, "5", $1); gsub(/^anop_X$/, "6", $1); gsub(/^Y_unplaced$/, "7", $1); print}' | bgzip > admixture_modified_v2.vcf.gz

tabix -p vcf admixture_modified_v2.vcf.gz
zgrep -v "^#" admixture_modified_v2.vcf.gz | cut -f1 | sort | uniq
bcftools query -l admixture_modified_v2.vcf.gz

# STEP 3 MAKE BED AND BIM FILES

plink --vcf admixture_modified_v2.vcf.gz --set-missing-var-ids @:# --keep-allele-order --const-fid --allow-extra-chr --make-bed --out melas_global_gambiaealigned

# STEP 4 RUN ADMIXTURE
# Make a file called K_runs.txt with values 1 to 10 and then run the following command
# Run this in screen as it takes a while

cat K_runs.txt | xargs -I {} sh -c 'admixture --cv=10 -j20 -s 14062 melas_global_gambiaealigned.bed {} | tee log{}.cv10.seed14062.out'

# Inspect the CV values for the inflection point to check what the best K value is for your dataset

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
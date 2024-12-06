# mosquito_pipeline

Pipeline for whole genome sequence analysis of Anopheles mosquitoes.

Step 1: Acquire fastq reads of samples and use fastq2matrix to generate bam files and vcf files for each sample.

Step 2: Conduct basic statistics on samples to check quality.

Step 3: Make a genomics database vcf using the samples that you deem good enough quality to keep.

Step 4: Filter the genomics database vcf to retain only good quality SNPs.

Step 5: Conducting Principal Components Analysis with your filtered vcf

Step 6: Create PCA using this R script

Step 7: Calculating and plotting admixture

Step 8: Creating a maximum likelihood tree

Step 9: Calculating FST (using python and jupyter notebook)

Step 10: Calculating genetic diversity metrics, nucleotide diversity and tajimas D

Step 11: Selection, calculating Garud's H12 statistic (using python and jupyter notebook)

Step 12: Selection, calculating iHS (using python and jupyter notebook)

Step 13: Selection, calculating XPEHH (using python and jupyter notebook)

Step 14: Using DELLY to analyse structural variants

Other scripts:

calculate_n50.py can be used to calculate the n50 of a sequence, for example a reference genome.

check_sex.py can be used to check the sex of mosquito samples using the ratio of coverage between the X chromosome and an autosome.

chromo_coverage.py and create_coverage_plots_sambamba_generic.py are part of basic_statistics and can be used to assess the coverage of sequence data across the genome, including visualising this in plots.

generate_admix_barplot_colours.R is an R script used as part of 7.admixture.sh


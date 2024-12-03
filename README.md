# mosquito_pipeline

Pipeline for whole genome sequence analysis of Anopheles mosquitoes.

Step 1: Acquire fastq reads of samples and use fastq2matrix to generate bam files and vcf files for each sample.

Step 2: Conduct basic statistics on samples to check quality.

Step 3: Make a genomics database vcf using the samples that you deem good enough quality to keep.

Step 4: Filter the genomics database vcf to retain only good quality SNPs.


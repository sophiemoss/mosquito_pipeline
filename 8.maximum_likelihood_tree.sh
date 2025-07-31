## For a maximum likelihood tree we want to use fasta sequences, for example of the mitochondria

## Use the vcf2fasta.py script to make the fasta sequences. You have to input a multisample vcf.
## If you are just creating a mitochondrial tree, you need to first subset your vcf to just me the mitochondrial genome
## Then subset the reference to be just the mitochondrial genome.

python vcf2fasta.py \
        --vcf "filtered_vcf_file.vcf.gz" \
        --ref "path_to_reference_genome" \
        --threads 10 \
        --whole-genome

## Then you can use muscle to align the mitochondrial fasta files to each other
# muscle -align seqs.fa -output seqs.afa
# e.g.

muscle -align mito_only_F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.ann.fa -output 2022_gambiae_bijagos_mito.afa
muscle -align mito_only_FMISSING_MAF_AC0_DP5_GQ20_gatk_filtered_miss_40_mac_bi_snps_melas_2019_plusglobal.2023_07_25.genotyped.ann.fa -output 2019_melas_global_bijagos_mito.afa

## Alternatively, I usually download alignments to my local disk and view them in aliview, then trim the alignment
## Then reupload the trimmed alignments to server

## Once you have aligned the fasta files and trimmed if necessary, you can create your maximum likelihood tree!
## These are different to neighbour joining trees and take a bit longer to run (read about the differences)

# Run raxml-ng and plot
raxml-ng --all --msa Anopheles_mito_18022024_aligned_bijagos_with_global_melasv6.fa --model GTR --prefix Anopheles_mito_0823 --seed 706046 --bs-metric tbe --tree rand{10} --bs-trees 1000

# RaXml outputs the best fitting tree, and you can then upload this to iTOL to see it!

raxml-ng --all --msa cox1_bed_edit_combined_trimmed.afa --model GTR --prefix AnDarCOX1 --seed 183645 --bs-metric tbe --tree rand{10} --bs-trees 1000

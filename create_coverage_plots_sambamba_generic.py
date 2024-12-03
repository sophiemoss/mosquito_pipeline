import subprocess
import pandas as pd
import sys

# Check if chromosome names were provided as command-line arguments
if len(sys.argv) < 2:
    print("Usage: python script_name.py <chromosome1> <chromosome2> ...")
    sys.exit(1)

# Read chromosome names from command-line arguments
chromosome_names = sys.argv[1:]

# Log the chromosome names being used
print(f"Filtering for chromosomes: {', '.join(chromosome_names)}")

# Use the fastq2vcfsamples.txt file which was created using:
# ls *_1.fastq.gz | sed 's/_1.fastq.gz//' > fastq2vcfsamples.txt

# Execute the first bash command to calculate depth of coverage in windows across the genome
subprocess.run('for f in *mkdup.bam ; do sambamba depth window -w 50000 "$f" -o "$f.chromocov" ; done', shell=True)
subprocess.run('cat fastq2vcfsamples.txt | parallel -j 1 "mv {}.mkdup.bam.chromocov {}_chromocov.csv"', shell=True)

# Read the file containing sample names
with open('fastq2vcfsamples.txt', 'r') as file:
    sample_names = file.read().splitlines()

# Loop through each sample name to filter out contigs
for sample_name in sample_names:
    # Construct the file paths
    input_file = sample_name + '_chromocov.csv'
    output_file = sample_name + '_chromocov.filtered.csv'

    # Read the input file
    df = pd.read_csv(input_file, sep='\t')

    # Filter the data based on the input chromosome names
    sub = df[df['# chrom'].isin(chromosome_names)]

    # Export the filtered data to a new file
    sub.to_csv(output_file, sep='\t', index=None)

# Use the other Python script to create figures of coverage over the genome
subprocess.run('cat fastq2vcfsamples.txt | parallel -j 30 "python chromo_coverage.py --sample {} --csvfile {}_chromocov.filtered.csv" > chromocoverage.txt', shell=True)

# Single line outside Python
# cat fastq2vcfsamples.txt | parallel -j 30 "python chromo_coverage.py --sample {} --csvfile {}_chromocov.filtered.csv" > chromocoveragelog.txt

# %% 
import subprocess
import pandas as pd
working_directory = "/mnt/storage11/sophie/bijagos_mosq_wgs/2019_melas_fq2vcf_gambiae_aligned/genomics_database_melas2019plusglobal/genomics_database_melas2019plusglobal_vcf/melas_2019_plusglobal_filtering/ir_snps"
os.chdir(working_directory)

# %% Define the input file with positions (CHR and POS)
vcf_filename = "final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.ann.vcf.gz"  # Replace with your VCF file path

# %% Define the input file with gene positions
positions_file = "genes_of_interest.txt"

# %% Prepare a list to collect data
data = []

# %% Read positions from the input file
with open(positions_file, "r") as infile:
    next(infile)  # Skip header line
    for line in infile:
        columns = line.strip().split("\t")
        if len(columns) >= 4:
            gene, chr, start, stop = columns

            # Define the bcftools command
            cmd = f'bcftools query -r {chr}:{start}-{stop} -f"%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO\t[%SAMPLE:%GT;]\n" {vcf_filename}'

            # Execute the bcftools command
            try:
                result = subprocess.run(cmd, stdout=subprocess.PIPE, text=True, check=True, shell=True)
                bcftools_output = result.stdout.splitlines()

                for line in bcftools_output:
                    fields = line.split('\t')
                    Chrom = fields[0]
                    Pos = fields[1]
                    Ref = fields[2]
                    Alt = fields[3]
                    Filter = fields[4]
                    Info = fields[5]
                    sample_genotype_concatenated = fields[6]

                    # Split the concatenated string using the semicolon as a delimiter
                    sample_genotype_pairs = sample_genotype_concatenated.split(';')

                    # Process sample genotypes to find samples with specified genotypes
                    samples_with_genotypes = []
                    for sample_genotype in sample_genotype_pairs:
                        parts = sample_genotype.split(':')
                        if len(parts) == 2:
                            sample, genotype = parts
                            if genotype in ["0/0","0/1", "0|1", "1/1", "1|1"]:  # Include specified genotypes
                                samples_with_genotypes.append(f"{sample}:{genotype}")

                    # Add the data to the list
                    data.append([gene, Chrom, Pos, Ref, Alt, Filter, Info, '; '.join(samples_with_genotypes)])

            except subprocess.CalledProcessError as e:
                print(f"Error processing gene: {gene}, CHR: {chr}, range: {start}-{stop}. Error: {e}")

# %% Create DataFrame
df = pd.DataFrame(data, columns=["Gene", "Chrom", "Position", "Ref", "Alt", "Filter", "Info", "SamplesWithGenotypes"])

# %% Optionally, save the DataFrame to a file
output_file = "snp_info_from_filtered_vcf_v2.csv"
df.to_csv(output_file, index=False)
print("DataFrame created and saved.")

# %%

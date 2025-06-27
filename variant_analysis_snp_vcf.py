###______###
# This script is designed to identify statistically different proportions of SNPs between two different populations within a VCF
# For example, in this VCF, PBO exposed resistant and susceptible mosquitoes
# You may want to filter your VCF to only contain missense variants before running this script.

# %% 
from cyvcf2 import VCF #for parsing the vcf
import pandas as pd
import os
from scipy.stats import fisher_exact # for doing the statistical test to compare allele counts
from tqdm import tqdm
import subprocess
import gzip
from collections import defaultdict
import re

# %% Create function to count the variants in a vcf, for estimating the time this will take
def count_variants_fast(vcf_path):
    result = subprocess.run(
        f"bcftools view -H {vcf_path} | wc -l",
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    return int(result.stdout.strip())

# %% 
# Set working directory
os.chdir('/mnt/storage11/sophie/guinea_coluzzi/filtered_vcf/snps')
cwd = os.getcwd()
print(f"Starting script in {cwd}")


# %% 
# Define your sample groups by finding the names in the VCF header
pbo_samples = [f"pbo-{i}" for i in range(1, 11)]
sus_samples = [f"sus-{i}" for i in range(1, 11)]

# %% 
# Load VCF
vcf_path = "PBO_exposed_filtered_missense_variants.vcf.gz"  # Replace with your VCF file name
vcf = VCF(vcf_path)
print(f"Loaded vcf {vcf_path}")

# %% 
# Map the sample names to their indices
sample_to_idx = {sample: idx for idx, sample in enumerate(vcf.samples)}
pbo_idxs = [sample_to_idx[s] for s in pbo_samples if s in sample_to_idx]
sus_idxs = [sample_to_idx[s] for s in sus_samples if s in sample_to_idx]
print("Mapped sample names to indices")

# %% 
# Estimate time taken for running this
variant_count = count_variants_fast(vcf_path)

# Create a container for the results
results = []
print("Running analysis for each variant...")

# Process each variant
for variant in tqdm(vcf, total=variant_count, desc="Processing variants", dynamic_ncols=True): # iterates over every variant in the vcf
    chrom = variant.CHROM 
    pos = variant.POS
    ref = variant.REF
    alt = variant.ALT[0]  # First (and only) ALT allele, assuming biallelic variants
    gts = variant.genotypes  # List of tuples (allele1, allele2, phased)

    def count_alleles(idx_list):
        alt_alleles = 0
        called_genotypes = 0
        for i in idx_list:
            gt = gts[i][:2] # First two elements: allelle1 and allele2
            if -1 in gt:
                continue  # Skip this genotype as it is missing (-1)(./.)
            alt_alleles += gt.count(1) # Count how many 1s (ALT) there are 
            called_genotypes += 1
        return alt_alleles, called_genotypes

    # Count ALT alleles and called genotypes for each group
    pbo_alt, pbo_called = count_alleles(pbo_idxs)
    sus_alt, sus_called = count_alleles(sus_idxs)

    # Create 2x2 table for Fisher's exact test
    pbo_total = 2 * pbo_called
    sus_total = 2 * sus_called
    
    # Calculate frequencies
    pbo_freq = pbo_alt / (pbo_total) if pbo_called > 0 else None
    sus_freq = sus_alt / (sus_total) if sus_called > 0 else None

    pbo_ref = pbo_total - pbo_alt
    sus_ref = sus_total - sus_alt

    # Valid 2x2 table of ALT vs REF allele counts
    if pbo_called > 0 and sus_called > 0:
        contingency = [[pbo_alt, pbo_ref], [sus_alt, sus_ref]]
        _, pval = fisher_exact(contingency)
    else:
        pval = None


    # Add FILTER field
    filter_status = variant.FILTER or 'PASS'

    # Add INFO as string
    info_dict = dict(variant.INFO)
    info_str = ";".join([f"{k}={v}" for k, v in info_dict.items()])

    # Sample genotypes
    sample_gts = [
        "/".join(map(str, gts[i][:2])) if -1 not in gts[i][:2] else './.'
        for i in range(len(vcf.samples))
    ]
    sample_gt_str = ";".join([f"{sample}:{gt}" for sample, gt in zip(vcf.samples, sample_gts)])

    results.append({
        "CHROM": chrom,
        "POS": pos,
        "REF": ref,
        "ALT": alt,
        "FILTER": filter_status,
        "INFO": info_str,
        "SAMPLE_GT": sample_gt_str,
        "pbo_freq": pbo_freq,
        "sus_freq": sus_freq,
        "pbo_called": pbo_called,
        "sus_called": sus_called,
        "p_value": pval
    })

# Convert to DataFrame
df = pd.DataFrame(results)
df_sorted = df.sort_values(by="p_value", ascending=True)

# Annotate the dataframe with information from the GFF file

# %%
# Annotate variants with gene names from a GFF file

gff_path = "/mnt/storage11/sophie/reference_genomes/AcolN3_ncbi_dataset/data/GCF_943734685.1/acoln3_genomic.gff"  # <-- Update this path to your actual GFF file
print(f"Annotating variants using GFF file: {gff_path}")

# Build a list of gene regions per chromosome

gene_regions = defaultdict(list)

# Support gzipped or plain GFF
open_func = gzip.open if gff_path.endswith(".gz") else open

with open_func(gff_path, "rt") as gff_file:
    for line in gff_file:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) != 9 or parts[2].lower() != "gene":
            continue
        chrom, source, feature, start, end, score, strand, phase, attributes = parts
        start = int(start)
        end = int(end)

        # Store full attributes string (starting from ID=...)
        gene_regions[chrom].append((start, end, attributes))

# Function to find gene name by position
def find_gene(chrom, pos):
    if chrom not in gene_regions:
        return None
    for start, end, gene in gene_regions[chrom]:
        if start <= pos <= end:
            return gene
    return None

# Apply to DataFrame
df_sorted["gene_name"] = df_sorted.apply(lambda row: find_gene(row["CHROM"], row["POS"]), axis=1)
print("Annotation complete: added 'gene_name' column.")

# %% 
# Export to CSV
df_sorted.to_csv("missense_variant_frequencies_with_stats.csv", index=False)

print("Analysis complete. Output saved to 'missense_variant_frequencies_with_stats.csv'.")

# %%
# Filter to statistically significant variants (p < 0.05)
significant_df = df_sorted[df_sorted["p_value"] < 0.05]

if not significant_df.empty:
    significant_df.to_csv("significantly_different_missense_variants_with_stats.csv", index=False)
    print("Significant variants saved to 'significantly_different_missense_variants_with_stats.csv'.")
else:
    print("No significantly different variants between the populations chosen.")

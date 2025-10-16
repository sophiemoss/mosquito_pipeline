# %% IMPORTS
from cyvcf2 import VCF
import pandas as pd
import os
from scipy.stats import fisher_exact
from tqdm import tqdm
import subprocess
import gzip
from collections import defaultdict

# %% COUNT VARIANTS FAST
def count_variants_fast(vcf_path):
    result = subprocess.run(
        f"bcftools view -H {vcf_path} | wc -l",
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    return int(result.stdout.strip())

# %% SETUP
os.chdir('/mnt/storage11/sophie/guinea_coluzzi/structural_variants')
print(f"Current directory: {os.getcwd()}")

# Define sample groups
pbo_samples = [f"pbo-{i}" for i in range(1, 11)]
sus_samples = [f"sus-{i}" for i in range(1, 11)]

# %% LOAD VCF
vcf_path = "merged_genotyped_structural_variants_miss20.ann.vcf.gz"
vcf = VCF(vcf_path)
print(f"Loaded VCF: {vcf_path}")

# Sample index mapping
sample_to_idx = {sample: idx for idx, sample in enumerate(vcf.samples)}
pbo_idxs = [sample_to_idx[s] for s in pbo_samples if s in sample_to_idx]
sus_idxs = [sample_to_idx[s] for s in sus_samples if s in sample_to_idx]

variant_count = count_variants_fast(vcf_path)

# %% Define function to count alleles
def count_alleles(idx_list, gts):
    carriers = 0
    called = 0
    for i in idx_list:
        gt = gts[i][:2]
        if -1 in gt:
            continue
        if 1 in gt:
            carriers += 1
        called += 1
    return carriers, called

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

    # Calculate frequencies

    pbo_total = 2 * pbo_called
    sus_total = 2 * sus_called

    pbo_ref = pbo_total - pbo_alt
    sus_ref = sus_total - sus_alt

    pbo_freq = pbo_alt / (pbo_total) if pbo_called > 0 else None
    sus_freq = sus_alt / (sus_total) if sus_called > 0 else None

    # Fisher's Exact Test setup. Creates 2 x 2 contingency tables.
    if pbo_called > 0 and sus_called > 0:
        contingency = [[pbo_alt, pbo_ref], [sus_alt, sus_ref]]
        _, pval = fisher_exact(contingency)
    else:
        pval = None

    # Add FILTER field

    filter_status = variant.FILTER or "PASS"
    
    # Add info as string

    info_dict = dict(variant.INFO)
    info_str = ";".join([f"{k}={v}" for k, v in info_dict.items()])

    sv_type = info_dict.get("SVTYPE", alt.strip("<>"))  # fallback to ALT if missing


    # Collect genotypes
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
        "SVTYPE": sv_type,
        "FILTER": filter_status,
        "INFO": info_str,
        "SAMPLE_GT": sample_gt_str,
        "pbo_freq": pbo_freq,
        "sus_freq": sus_freq,
        "pbo_called": pbo_total,
        "sus_called": sus_total,
        "p_value": pval
    })

df = pd.DataFrame(results)
df_sorted = df.sort_values(by="p_value", ascending=True)

# %% GFF ANNOTATION
gff_path = "/mnt/storage11/sophie/reference_genomes/AcolN3_ncbi_dataset/data/GCF_943734685.1/acoln3_genomic.gff"
print(f"Annotating variants using GFF: {gff_path}")
gene_regions = defaultdict(list)
open_func = gzip.open if gff_path.endswith(".gz") else open

with open_func(gff_path, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) != 9 or parts[2].lower() != "gene":
            continue
        chrom, _, _, start, end, _, _, _, attr = parts
        gene_regions[chrom].append((int(start), int(end), attr))

def find_gene(chrom, pos):
    if chrom not in gene_regions:
        return None
    for start, end, attr in gene_regions[chrom]:
        if start <= pos <= end:
            return attr
    return ""

df_sorted["gene_name"] = df_sorted.apply(lambda row: find_gene(row["CHROM"], row["POS"]), axis=1)
print("Annotation complete.")

#Reorder dataframe

cols = [
    "CHROM", "POS", "REF", "ALT", "SVTYPE", "FILTER", "INFO", "SAMPLE_GT",
    "pbo_freq", "sus_freq", "pbo_called", "sus_called", "p_value", "gene_name"
]

df_sorted = df_sorted[cols]


# %% SAVE RESULTS
df_sorted.to_csv("structural_variant_frequencies_with_stats.csv", index=False)
print("Saved: structural_variant_frequencies_with_stats.csv")

sig_df = df_sorted[df_sorted["p_value"] < 0.05]
if not sig_df.empty:
    sig_df.to_csv("significantly_different_structural_variants.csv", index=False)
    print("Saved: significantly_different_structural_variants.csv")
else:
    print("No significant differences found.")

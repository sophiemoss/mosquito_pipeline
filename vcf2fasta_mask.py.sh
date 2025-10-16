# From a multi-sample gVCF, create per-sample masks for areas with poor coverage
# Make sure that low depth variants are set to missing in the VCF (this script doesn't mask on DP, it only masks on ./.)
import os
import subprocess
import argparse
from itertools import islice

def run_cmd(cmd):
    print(f"[RUN] {cmd}")
    subprocess.run(cmd, shell=True, check=True)

def batch_iterator(iterable, batch_size):
    """Yield successive batches from iterable"""
    iterable = iter(iterable)
    while True:
        batch = list(islice(iterable, batch_size))
        if not batch:
            break
        yield batch

def main():
    parser = argparse.ArgumentParser(description="Per-sample masked consensus with batch processing")
    parser.add_argument("--vcf", required=True, help="Input multisample VCF")
    parser.add_argument("--ref", required=True, help="Reference FASTA")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--batch_size", type=int, default=5, help="Number of samples per batch")
    parser.add_argument("--tree", action="store_true", help="Run IQTree on combined FASTA")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # --- Get sample names from VCF ---
    samples = subprocess.check_output(
        f"bcftools query -l {args.vcf}", shell=True
    ).decode().splitlines()
    samples = [s.strip() for s in samples]

    masked_fastas = []

    for batch in batch_iterator(samples, args.batch_size):
        print(f"[INFO] Processing batch: {batch}")

        for sample in batch:
            print(f"[INFO] Processing sample {sample}")

            consensus_fasta = os.path.join(args.outdir, f"{sample}.consensus.fa")
            low_depth_bed = os.path.join(args.outdir, f"{sample}.lowdepth.bed")
            het_bed = os.path.join(args.outdir, f"{sample}.het.bed")
            mask_bed = os.path.join(args.outdir, f"{sample}.mask.bed")
            masked_fasta = os.path.join(args.outdir, f"{sample}.masked.fa")

            # 1. Generate sample-specific consensus FASTA
            run_cmd(
                f"bcftools consensus -f {args.ref} -s {sample} {args.vcf} > {consensus_fasta}"
            )

            # 2. Low-depth positions (./.)
            run_cmd(
                f"bcftools query -s {sample} -f '%CHROM\\t%POS0\\t%POS\\t[%GT]\\n' {args.vcf} | "
                f"awk '$4==\"./.\" {{print $1\"\\t\"$2\"\\t\"$3}}' > {low_depth_bed}"
            )

            # 3. Heterozygous calls: exclude homozygous ref (0/0), homozygous alt (1/1), and missing (./.)
            run_cmd(
                f"bcftools query -s {sample} -f '%CHROM\\t%POS0\\t%POS\\t[%GT]\\n' {args.vcf} | "
                f"awk '$4!=\"0/0\" && $4!=\"1/1\" && $4!=\"./.\" {{print $1\"\\t\"$2\"\\t\"$3}}' > {het_bed}"
            )

            # 4. Merge mask regions
            run_cmd(
                f"cat {low_depth_bed} {het_bed} | sort -k1,1 -k2,2n | bedtools merge -i - > {mask_bed}"
            )

            # 5. Mask consensus FASTA
            run_cmd(
                f"bedtools maskfasta -fi {consensus_fasta} -bed {mask_bed} -fo {masked_fasta}"
            )

            # 6. Fix FASTA headers -> sample_chrX
            run_cmd(f"sed -i 's/^>/>{sample}_/' {masked_fasta}")

            masked_fastas.append(masked_fasta)

    # 7. Concatenate all masked FASTAs
    combined_fasta = os.path.join(args.outdir, "combined_masked_consensus.fa")
    run_cmd(f"cat {' '.join(masked_fastas)} > {combined_fasta}")

    # 8. Optional IQTree
    if args.tree:
        run_cmd(f"iqtree -s {combined_fasta} -m GTR+G+ASC -nt AUTO")

if __name__ == "__main__":
    main()

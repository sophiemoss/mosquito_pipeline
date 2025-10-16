import sys
import argparse
import fastq2matrix as fm
import os

def main(args):
    args.region_arg = ""

    # 1) Variant calling
    if args.variant_caller == "gatk":
        if args.bed:
            args.region_arg = "-L %s" % args.bed
        fm.run_cmd("gatk HaplotypeCaller -R %(ref)s %(region_arg)s -I %(bam)s -O %(out)s.vcf.gz" % vars(args))

    elif args.variant_caller == "bcftools":
        if args.bed:
            args.region_arg = "-R %s" % args.bed  # mpileup supports -R even on older builds
        fm.run_cmd("bcftools mpileup -f %(ref)s %(region_arg)s %(bam)s | bcftools call -mv -Oz -o %(out)s.vcf.gz" % vars(args))

    elif args.variant_caller == "freebayes":
        if args.bed:
            args.region_arg = "-t %s" % args.bed
        fm.run_cmd("freebayes -f %(ref)s %(region_arg)s %(bam)s | bgzip -c > %(out)s.vcf.gz" % vars(args))
    else:
        quit("Unknown variant caller! Exiting!")

    # 2) Index VCF
    fm.run_cmd("tabix -f %(out)s.vcf.gz" % vars(args))

    # 3) Build depth mask (positions with depth < cutoff -> mask to N)
    #    Clamp negatives and ensure start<end to avoid consensus crashes at start-of-interval
    if args.bed:
        fm.run_cmd(
            "bedtools coverage -a %(bed)s -b %(bam)s -d | "
            "awk -v c=%(depth_cutoff)s '{ "
            "  if($NF<c){ s=$2+$(NF-1)-2; e=$2+$(NF-1)-1; "
            "    if(s<0)s=0; if(e<=s)e=s+1; "
            "    print $1\"\\t\"s\"\\t\"e "
            "  } "
            "}' > %(out)s.depth_mask.bed" % vars(args)
        )
    else:
        fm.run_cmd(
            "bedtools genomecov -ibam %(bam)s -d | "
            "awk -v c=%(depth_cutoff)s '{ "
            "  if($NF<c){ s=$2-1; e=$2; if(s<0)s=0; if(e<=s)e=s+1; "
            "    print $1\"\\t\"s\"\\t\"e "
            "  } "
            "}' > %(out)s.depth_mask.bed" % vars(args)
        )

    # 4) Decide whether to pass a mask
    num_lines = 0
    for l in fm.cmd_out("wc -l %(out)s.depth_mask.bed" % vars(args)):
        num_lines = int(l.strip().split()[0])
    args.mask_arg = "-m %(out)s.depth_mask.bed -M N" % vars(args) if num_lines > 0 else ""

    # 5) Region names (BED col 4) and regions file for faidx (1-based inclusive)
    region_names = {}
    if args.bed:
        regions_file = args.out + ".regions.txt"
        with open(regions_file, "w") as O:
            for l in open(args.bed):
                if not l.strip() or l.startswith(("#", "track", "browser")):
                    continue
                row = l.strip().split()
                if len(row) < 3:
                    continue
                chrom, start0, end0 = row[0], int(row[1]), int(row[2])
                start1 = start0 + 1     # 1-based inclusive for samtools faidx
                end1   = end0
                r1 = f"{chrom}:{start1}-{end1}"
                O.write(r1 + "\n")
                # store both possible header styles to catch either case
                if len(row) > 3:
                    name = row[3]
                    region_names[f"{chrom}:{start0}-{end0}"] = name   # 0-based style
                    region_names[r1] = name                            # 1-based style

        # ---- FIXED CONSENSUS COMMAND (no '-r', uses xargs) ----
        consensus_cmd = (
            "xargs -a {regions} -I{{}} samtools faidx {ref} \"{{}}\" | "
            "bcftools consensus {vcf} {mask} -M N"
        ).format(
            regions=regions_file,
            ref=args.ref,
            vcf=f"{args.out}.vcf.gz",
            mask=args.mask_arg
        )
    else:
        # Whole-reference consensus: stream all contigs via faidx, no -R needed
        fai = args.ref + ".fai"
        if not os.path.exists(fai):
            fm.run_cmd("samtools faidx %(ref)s" % vars(args))
        consensus_cmd = (
            "cut -f1 {fai} | xargs -I{{}} samtools faidx {ref} \"{{}}\" | "
            "bcftools consensus {vcf} {mask} -M N"
        ).format(
            fai=fai,
            ref=args.ref,
            vcf=f"{args.out}.vcf.gz",
            mask=args.mask_arg
        )

    # 6) Write consensus FASTA, rewriting headers with friendly names
    with open(args.out + ".consensus.fa", "w") as O:
        for l in fm.cmd_out(consensus_cmd):
            if not l:
                continue
            if l[0] == ">":
                r = l.strip()[1:]
                O.write(">%s %s\n" % (args.out, region_names.get(r, r)))
            else:
                O.write(l if l.endswith("\n") else l + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='tbprofiler script', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--bam', type=str, help='Bam file', required=True)
    parser.add_argument('--ref', type=str, help='Reference file', required=True)
    parser.add_argument('--out', type=str, help='Output prefix', required=True)
    parser.add_argument('--bed', type=str, help='Bed file with regions to be extracted (optional 4th column to contain region name which will be present in consensus file)')
    parser.add_argument('--variant-caller', choices=["gatk","freebayes","bcftools"], default="freebayes", type=str, help='Variant caller to be used')
    parser.add_argument('--depth-cutoff', default=10, type=int, help='Sites with less than this cutoff will be marked as missing in the consensus file')
    parser.set_defaults(func=main)

    args = parser.parse_args()
    args.func(args)

# Calculate heterozygosity from a filtered VCF, per sample, in 1000 bp windows.
# Converts VCF into PLINK format, then run --het in 1000 bp windows in a loop.

#!/usr/bin/env bash
set -euo pipefail

# Usage: ./quick_window_het.sh <input.vcf[.gz]> <out_prefix> [window_bp=1000000]
VCF="${1:?Usage: $0 <input.vcf[.gz]> <out_prefix> [window_bp=1000000]}"
OUT="${2:?Usage: $0 <input.vcf[.gz]> <out_prefix> [window_bp=1000000]}"
WIN="${3:-1000000}"

for x in plink bcftools awk Rscript; do
  command -v "$x" >/dev/null 2>&1 || { echo "ERROR: $x not found in PATH"; exit 1; }
done

WORK="${OUT}.win${WIN}.work"
mkdir -p "$WORK"

# 1) Index VCF if needed
[ -f "${VCF}.csi" ] || [ -f "${VCF}.tbi" ] || bcftools index -c "$VCF"

# 2) Convert VCF → PLINK binary once
if [ ! -f "$WORK/data.bed" ]; then
  echo "[*] Converting VCF to PLINK binary..."
  plink \
    --vcf "$VCF" \
    --double-id \
    --snps-only just-acgt \
    --biallelic-only strict \
    --make-bed \
    --allow-extra-chr \
    --out "$WORK/data"
fi

# 3) Output header
OUTTSV="${OUT}.window_${WIN}bp.heterozygosity.tsv"
echo -e "CHR\tSTART\tEND\tFID\tIID\tN_SITES\tO_HOM\tE_HOM\tF\tHET_OBS" > "$OUTTSV"

echo "[*] Scanning genome in ${WIN} bp windows..."
bcftools index -s "$VCF" | awk '$1!~/^#/{print $1"\t"$2}' | while read -r CHR LEN; do
  START=0
  while [ "$START" -lt "$LEN" ]; do
    END=$(( START + WIN ))
    [ "$END" -gt "$LEN" ] && END="$LEN"

    plink --bfile "$WORK/data" \
          --chr "$CHR" \
          --from-bp $((START+1)) \
          --to-bp "$END" \
          --het --allow-extra-chr \
          --out "$WORK/het_${CHR}_${START}_${END}" >/dev/null 2>&1 || true

    HETFILE="$WORK/het_${CHR}_${START}_${END}.het"
    if [ -s "$HETFILE" ]; then
      awk -v C="$CHR" -v S="$START" -v E="$END" 'NR>1{
        FID=$1; IID=$2; O=$3; Eexp=$4; N=$5; F=$6;
        het = (N>0 ? 1-(O/N) : "NA");
        printf "%s\t%s\t%s\t%s\t%s\t%d\t%d\t%.6f\t%.6f\t%s\n",
               C,S,E,FID,IID,N,O,Eexp,F,het;
      }' "$HETFILE" >> "$OUTTSV"
      rm -f "$HETFILE" "$WORK/het_${CHR}_${START}_${END}.log"
    fi

    START=$(( END ))
  done
done

echo "[*] Heterozygosity table: $OUTTSV"

# 4) Plot with R using RColorBrewer, continuous genome + chromosome markers
PLOTPDF="${OUT}.heterozygosity.pdf"
echo "[*] Plotting per-sample heterozygosity with chromosome markers -> $PLOTPDF"

export OUTTSV
export PLOTPDF

Rscript - <<'EOF'
suppressPackageStartupMessages(library(RColorBrewer))

tbl <- read.table(Sys.getenv("OUTTSV"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
if (nrow(tbl) == 0L) {
  stop("No data rows found in heterozygosity table.")
}

# Midpoint of each window
tbl$WindowMid <- (tbl$START + tbl$END)/2

# Keep chromosome order as it appears in the table (which follows bcftools index -s)
chr_levels <- unique(tbl$CHR)
tbl$CHR <- factor(tbl$CHR, levels = chr_levels)

# Chromosome lengths (use max END per chr), offsets to build continuous genome axis
chr_lengths <- tapply(tbl$END, tbl$CHR, max)
chr_lengths <- as.numeric(chr_lengths)              # ensure numeric
names(chr_lengths) <- chr_levels

chr_offsets <- c(0, cumsum(head(chr_lengths, -1)))
names(chr_offsets) <- chr_levels

# Global genomic coordinate for plotting
tbl$GlobalMid <- tbl$WindowMid + chr_offsets[as.character(tbl$CHR)]

# Set up colours per sample
samples <- unique(tbl$IID)
n <- length(samples)
if (n <= 12) {
  cols <- brewer.pal(n, "Set3")
} else {
  cols <- colorRampPalette(brewer.pal(12, "Set3"))(n)
}
names(cols) <- samples

# Plot
pdf(Sys.getenv("PLOTPDF"), width=12, height=6)
plot(NULL,
     xlim = range(tbl$GlobalMid, na.rm=TRUE),
     ylim = c(0, 1),
     xlab = "Genomic coordinate (continuous)",
     ylab = "Observed heterozygosity",
     main = "Per-sample heterozygosity (windowed)",
     xaxt = "n")

# Draw chromosome boundary lines and x-axis labels at chr midpoints
boundary_positions <- chr_offsets[-1]  # exclude 0
abline(v = boundary_positions, lty = 3)

chr_midpoints <- chr_offsets + (chr_lengths / 2)
axis(1, at = chr_midpoints, labels = chr_levels, las = 2, cex.axis = 0.7)

# Draw lines per sample
for (s in samples) {
  d <- tbl[tbl$IID == s & !is.na(tbl$HET_OBS), ]
  # Order by global coordinate to avoid zig-zags
  d <- d[order(d$GlobalMid), ]
  lines(d$GlobalMid, d$HET_OBS, col = cols[s], type = "l")
}

legend("topright", legend = samples, col = cols, lty = 1, cex = 0.7, bty = "n", ncol = 1)
dev.off()
EOF

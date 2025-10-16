# Collect bamfile, index them
for bam in *.bam; do
    echo "[INFO] Indexing $bam"
    samtools index "$bam"
done

# Splice the bam files for the region we're interested in
for bam in *.bam; do
    base=$(basename "$bam" .bam)
    samtools view -hb "$bam" "AdarGF1_3RL:52164884-52166357" > "${base}.pf47rec.bam"
    samtools index "${base}.region.bam"
done

# Convert bam to consensus fasta files
# Running now for darlingi 
Region: NC_064875.1:52164884-52166357 

python /mnt/storage11/sophie/gitrepos/mosquito_pipeline/bam2consensus_v2.py \
  --bam /mnt/storage11/sophie/pf47/anopheles_pfs47Rec/darlingi_bam_files/AnDar_008.mkdup.bam \
  --ref AnoDarl_H01.genomic.fasta \
  --out AnDar_008 \
  --bed targets.bed \
  --variant-caller freebayes \
  --depth-cutoff 10

# Looped

for bam in /mnt/storage11/sophie/pf47/anopheles_pfs47Rec/darlingi_bam_files/*.mkdup.bam
do
    sample=$(basename "$bam" .mkdup.bam)
    echo "Processing $sample"

    python /mnt/storage11/sophie/gitrepos/mosquito_pipeline/bam2consensus_v2.py \
        --bam "$bam" \
        --ref AnoDarl_H01.genomic.fasta \
        --out "$sample" \
        --bed targets.bed \
        --variant-caller freebayes \
        --depth-cutoff 10
done

# Run for gambiae s.l.

Region: 2L:31084065-31086356

python /mnt/storage11/sophie/gitrepos/mosquito_pipeline/bam2consensus_v2.py \
  --bam /mnt/storage11/sophie/pf47/anopheles_pfs47Rec/darlingi_bam_files/AnDar_008.mkdup.bam \
  --ref AnoDarl_H01.genomic.fasta \
  --out AnDar_008 \
  --bed targets.bed \
  --variant-caller freebayes \
  --depth-cutoff 10

# Looped

for bam in /mnt/storage11/sophie/pf47/anopheles_pfs47Rec/gambiae_bam_files/*.mkdup.bam
do
    sample=$(basename "$bam" .mkdup.bam)
    echo "Processing $sample"

    python /mnt/storage11/sophie/gitrepos/mosquito_pipeline/bam2consensus_v2.py \
        --bam "$bam" \
        --ref Anopheles_gambiae.AgamP4.dna.toplevel.fa \
        --out "$sample" \
        --bed targets.bed \
        --variant-caller freebayes \
        --depth-cutoff 10
done

# create a directory called bam_coverage_output which mosdepth will put its results in

#!/bin/bash

# Define the BED file containing the regions
regions_bed="AgamP4_chr.bed"

# Read each BAM file from bamlist.txt and run mosdepth
while read bamfile; do
    # Extract the sample name from the BAM file path (assuming it's the file name without the extension)
    sample_name=$(basename "$bamfile" .bam)
    
    # Run mosdepth for the current sample
    mosdepth --by $regions_bed "bam_coverage_output/${sample_name}" "$bamfile"
done < bamlist.txt

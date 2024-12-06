import pandas as pd
import numpy as np
import os
import glob

# Correct path to the directory containing the mosdepth outputs
output_dir = "bam_coverage_output"

# List all .regions.bed.gz files produced by mosdepth
coverage_files = glob.glob(os.path.join(output_dir, "*.regions.bed.gz"))

# Initialize an empty DataFrame for aggregating coverage data
coverage_data = pd.DataFrame()

for file_path in coverage_files:
    # Extract the sample name from the file path
    sample_name = os.path.basename(file_path).split(".")[0]
    # Read the coverage data from the gzipped BED file
    temp_df = pd.read_csv(file_path, sep='\t', header=None, compression='gzip')
    # Correctly assign column names, including a placeholder for the region name/ID
    temp_df.columns = ['chrom', 'start', 'end', 'region', sample_name]
    # If it's the first file, initialize the coverage_data DataFrame with it
    if coverage_data.empty:
        coverage_data = temp_df
    else:
        # Merge the new data with the existing DataFrame on chrom, start, end, and region
        coverage_data = pd.merge(coverage_data, temp_df, on=['chrom', 'start', 'end', 'region'], how='outer')

# Calculate the average coverage across all samples for each region
coverage_data['Average_Coverage'] = coverage_data.drop(columns=['chrom', 'start', 'end', 'region']).mean(axis=1)

# Specify the output file path
output_file_path = "final_coverage_report.csv"

# Save the aggregated coverage data to a CSV file
coverage_data.to_csv(output_file_path, index=False, sep='\t')

print(f"Coverage report generated: {output_file_path}")

# %% 
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Read your dataset (assuming it's a CSV file)
df = pd.read_csv("27_05_2025_hmmIBD_fraction.csv")

# Create a pivot table for the heatmap
heatmap_data = df.pivot(index="sample1", columns="sample2", values="fraction")

# Set ggplot style
plt.style.use("ggplot")

# Create the heatmap with plasma colormap
plt.figure(figsize=(12, 10))
sns.heatmap(heatmap_data, cmap="plasma", annot=False, fmt=".2f", linewidths=0.5, linecolor='gray')

# Customize the plot
plt.title("IBD Fraction Heatmap")
plt.xlabel("Sample 2")
plt.ylabel("Sample 1")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()

# Show the plot
plt.show()

# %% Format differently
#### With full symmetry and only lower triangle

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Load data
df = pd.read_csv("27_05_2025_hmmIBD_fraction.csv")  # Update this to your file name

# Create a list of all unique samples
all_samples = pd.unique(df[['sample1', 'sample2']].values.ravel('K'))

# Initialize empty DataFrame
matrix = pd.DataFrame(index=all_samples, columns=all_samples, dtype=float)

# Fill in the matrix with fractions
for _, row in df.iterrows():
    s1, s2, val = row['sample1'], row['sample2'], row['fraction']
    matrix.loc[s1, s2] = val
    matrix.loc[s2, s1] = val  # Ensure symmetry by filling in for the opposite comparisons too

# Create a mask for the upper triangle so that only the lower triangle shows
mask = np.triu(np.ones_like(matrix, dtype=bool))

# Plot
plt.figure(figsize=(14, 12))
sns.heatmap(matrix, cmap='plasma', mask=mask, square=True, cbar_kws={"shrink": 0.8}, linewidths=0.5, linecolor='gray')

plt.title("IBD Fraction Heatmap")
plt.xlabel("Sample 2")
plt.ylabel("Sample 1")
plt.xticks(rotation=45, ha="right")
plt.yticks(rotation=0)
plt.tight_layout()
plt.show()

# %%

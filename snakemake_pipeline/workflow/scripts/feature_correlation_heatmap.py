import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

corr_matrix = pd.read_csv("microbe_gene_count_normalized_matrix.csv", index_col=0)

threshold = 0.8

# # Mask rows and columns where all values are <= 0.7
# rows_to_keep = (corr_matrix > threshold).any(axis=1)  # Keep rows where any value > 0.7
# cols_to_keep = (corr_matrix > threshold).any(axis=0)  # Keep columns where any value > 0.7

# # Step 2: Subset the correlation matrix to include only significant rows and columns
# filtered_corr_matrix = corr_matrix.loc[rows_to_keep, cols_to_keep]

corr_matrix_masked = corr_matrix.copy()
corr_matrix_masked[corr_matrix <= threshold] = np.nan  # Set values <= threshold to NaN

max_val = corr_matrix_masked.max().max()

# Step 3: Plot the heatmap
plt.figure(figsize=(15, 15))
sns.heatmap(corr_matrix_masked, 
            cmap='Blues', 
            annot=False, 
            fmt=".2f", 
            vmin=threshold, 
            vmax=max_val, 
            cbar_kws={'label': 'Correlation'},
            square=True)

plt.title("Heatmap of Normalized Microbe-Gene Correlations")
plt.show()

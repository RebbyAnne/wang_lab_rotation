import pandas as pd
import sys
import matplotlib.pyplot as plt
import seaborn as sns

# DISTANCE_MATRIX = sys.argv[1]  
# GENE_FEATURES_MATRIX_F = sys.argv[2]
# MICROBE_FEATURES_MATRIX_F = sys.argv[3]
DISTANCE_MATRIX = "workflow/out/node_distances/noIsolates_4Lw100Nw0.1p10q128k_featureDistances.csv"
GENE_FEATURES_MATRIX_F = "workflow/out/matrices/gene_features_matrix.csv"
MICROBE_FEATURES_MATRIX_F = "workflow/out/matrices/microbe_features_matrix.csv"

dist_threshold = 3.7

def main(): 
    distance_df = pd.read_csv(DISTANCE_MATRIX, index_col=0)

    gene_features_matrix = pd.read_csv(GENE_FEATURES_MATRIX_F, index_col=0)
    microbe_features_matrix = pd.read_csv(MICROBE_FEATURES_MATRIX_F, index_col=0)

    gene_features = gene_features_matrix.columns.tolist()
    microbe_features = microbe_features_matrix.columns.tolist()

    matching_genes = [gene for gene in gene_features if gene in distance_df.columns]
    matching_microbes = [microbe for microbe in microbe_features if microbe in distance_df.columns]

    df = distance_df.loc[matching_genes, matching_microbes]
    # Remove columns where all values are >= 5
    df_filtered_cols = df.loc[:, (df < dist_threshold).any()]
    # Remove rows where all values are >= 5
    df_filtered = df_filtered_cols[(df_filtered_cols < dist_threshold).any(axis=1)]

    min_val = df_filtered.min().min()

    # Create a clustered heatmap with better label spacing
    g = sns.clustermap(df_filtered, 
                    cmap='Blues_r',  # Dark for low, white for high
                    # annot=False, 
                    # fmt=".2f", 
                    vmin=min_val,  
                    vmax=dist_threshold, 
                    cbar_kws={'label': 'Distance'},
                    method='ward',  # Clustering method
                    metric='euclidean',  # Distance metric
                    dendrogram_ratio=(0.1, 0.1),  # Reduce dendrogram size
                    xticklabels=True, 
                    yticklabels=True, 
                    figsize=(50,50))

    # Adjust label rotation and spacing
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right", fontsize=10)
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=10)

    plt.legend(loc='best') 

    # Adjust layout to prevent clipping
    # g.figure.subplots_adjust(bottom=0.3, left=0.3, top=0.95, right=0.95)

    plt.show()


if __name__ == '__main__':
    main()
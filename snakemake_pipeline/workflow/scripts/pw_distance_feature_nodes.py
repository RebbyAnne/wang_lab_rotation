from gensim.models import KeyedVectors
from scipy.spatial.distance import pdist, squareform
import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
import seaborn as sns

EMBEDDING_F = sys.argv[1]
GENE_FEATURES_MATRIX_F = sys.argv[2]
MICROBE_FEATURES_MATRIX_F = sys.argv[3]

DISTANCES_F = sys.argv[4]  

DIST_THRESHOLD = sys.argv[5]

def main():
    embedding_kv = KeyedVectors.load(EMBEDDING_F)

    gene_features_matrix = pd.read_csv(GENE_FEATURES_MATRIX_F, index_col=0)
    microbe_features_matrix = pd.read_csv(MICROBE_FEATURES_MATRIX_F, index_col=0)

    gene_features = gene_features_matrix.columns.tolist()
    microbe_features = microbe_features_matrix.columns.tolist()

    gene_feature_in_embeddings = [node for node in embedding_kv.index_to_key if node in gene_features]
    microbe_feature_in_embeddings = [node for node in embedding_kv.index_to_key if node in microbe_features]
    features_in_embeddings = gene_feature_in_embeddings + microbe_feature_in_embeddings

    embedding_matrix = np.array([embedding_kv[node] for node in features_in_embeddings])

    # Step 2: Calculate pairwise distances
    distance_vector = pdist(embedding_matrix, metric='euclidean')
    distance_matrix = squareform(distance_vector) # square matrix with pairwise distances

    distance_df = pd.DataFrame(distance_matrix, index=features_in_embeddings, columns=features_in_embeddings)
    distance_df.to_csv(DISTANCES_F, float_format='%.6f')    

    # df = distance_df.loc[gene_feature_in_embeddings, microbe_feature_in_embeddings]
    # # Remove columns where all values are >= 5
    # df_filtered_cols = df.loc[:, (df < dist_threshold).any()]
    # # Remove rows where all values are >= 5
    # df_filtered = df_filtered_cols[(df_filtered_cols < dist_threshold).any(axis=1)]

    # min_val = df_filtered.min().min()

    # plt.figure(figsize=(15, 15))
    # sns.heatmap(df_filtered, 
    #             cmap='Blues', 
    #             annot=False, 
    #             fmt=".2f", 
    #             vmin=0, 
    #             vmax=dist_threshold, 
    #             cbar_kws={'label': 'Correlation'},
    #             square=True)

    # plt.title("Heatmap of Microbe-Gene Node Distances")
    # plt.show()


if __name__ == '__main__':
    main()

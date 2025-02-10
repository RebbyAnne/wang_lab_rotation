from gensim.models import KeyedVectors
from scipy.spatial.distance import pdist, squareform
import numpy as np
import pandas as pd
import sys
from sklearn import decomposition
from itertools import combinations

EMBEDDING_F = sys.argv[1]
GENE_FEATURES_MATRIX_F = sys.argv[2]
MICROBE_FEATURES_MATRIX_F = sys.argv[3]
DISTANCES_F = sys.argv[4]  
PAIRWISE_DISTANCES_CSV = sys.argv[5]  # New argument for output CSV

N_COMPONENTS = int(sys.argv[6])

DIMENSIONS = int(sys.argv[7])
WALK_LEN = int(sys.argv[8])
N_WALKS = int(sys.argv[9])
PVAL = float(sys.argv[10])
QVAL = float(sys.argv[11])

def main():
    # Load embeddings
    embedding_kv = KeyedVectors.load(EMBEDDING_F)

    # Load feature matrices
    gene_features_matrix = pd.read_csv(GENE_FEATURES_MATRIX_F, index_col=0)
    microbe_features_matrix = pd.read_csv(MICROBE_FEATURES_MATRIX_F, index_col=0)

    # Get feature names
    gene_features = gene_features_matrix.columns.tolist()
    microbe_features = microbe_features_matrix.columns.tolist()

    # Identify features present in embeddings
    gene_feature_nodes = [node for node in embedding_kv.index_to_key if node in gene_features]
    microbe_feature_nodes = [node for node in embedding_kv.index_to_key if node in microbe_features]

    feature_nodes = gene_feature_nodes + microbe_feature_nodes

    # Extract feature vectors
    vectors = np.array([embedding_kv[node] for node in feature_nodes])

    # Perform PCA
    pca_model = decomposition.PCA(n_components=N_COMPONENTS)
    reduced_vectors = pca_model.fit_transform(vectors)

    # Compute pairwise distances
    distance_vector = pdist(reduced_vectors, metric='euclidean')
    distance_matrix = squareform(distance_vector)

    # Convert distance matrix to DataFrame
    distance_df = pd.DataFrame(distance_matrix, index=feature_nodes, columns=feature_nodes)
    distance_df.to_csv(DISTANCES_F, float_format='%.6f')    

    # Use itertools.combinations for unique pairs
    pairwise_data = []
    for feature1, feature2 in combinations(feature_nodes, 2):
        distance = distance_matrix[feature_nodes.index(feature1), feature_nodes.index(feature2)]
        pairwise_data.append([feature1, feature2, distance, DIMENSIONS, PVAL, QVAL, WALK_LEN, N_WALKS])

    pairwise_df = pd.DataFrame(pairwise_data, columns=[
        "feature1", "feature2", "distance", "dimensions", "p", "q", "walk_len", "n_walks"
    ])

    pairwise_df.to_csv(PAIRWISE_DISTANCES_CSV, index=False, float_format='%.6f')


if __name__ == '__main__':
    main()

from gensim.models import KeyedVectors
from sklearn import decomposition
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys

EMBEDDING_F = sys.argv[1]
GENE_FEATURES_MATRIX_F = sys.argv[2]
MICROBE_FEATURES_MATRIX_F = sys.argv[3]
PLOT_F = sys.argv[4]

N_COMPONENTS = int(sys.argv[5])
ALPHA = float(sys.argv[6])
SIZE = int(sys.argv[7])

DIMENSIONS = str(sys.argv[8])
WALK_LEN = str(sys.argv[9])
N_WALKS = str(sys.argv[10])
P = str(sys.argv[11])
Q = str(sys.argv[12])

def main():
    embedding_kv = KeyedVectors.load(EMBEDDING_F)

    gene_features_matrix = pd.read_csv(GENE_FEATURES_MATRIX_F, index_col=0)
    microbe_features_matrix = pd.read_csv(MICROBE_FEATURES_MATRIX_F, index_col=0)

    genes = gene_features_matrix.index.tolist()
    gene_features = gene_features_matrix.columns.tolist()

    microbes = microbe_features_matrix.index.tolist()
    microbe_features = microbe_features_matrix.columns.tolist()

    gene_nodes = [node for node in embedding_kv.index_to_key if node in genes]
    microbe_nodes = [node for node in embedding_kv.index_to_key if node in microbes]

    gene_feature_nodes = [node for node in embedding_kv.index_to_key if node in gene_features]
    microbe_feature_nodes = [node for node in embedding_kv.index_to_key if node in microbe_features]

    vectors = np.array([embedding_kv[node] for node in embedding_kv.index_to_key])

    pca_model = decomposition.PCA(n_components=N_COMPONENTS)
    pca_results = pca_model.fit_transform(vectors)

    plt.figure(figsize=(10, 10))
    x = pca_results[:, 0]
    y = pca_results[:, 1]

    colors = ["#540d6e", "#ee4266", "#ffd23f", "#3bceac"]

    # Plot gene nodes
    # gene_indices = [embedding_kv.index_to_key.index(node) for node in gene_nodes]
    # plt.scatter(x[gene_indices], y[gene_indices], 
    #             alpha=ALPHA, 
    #             linewidth=0, 
    #             c=colors[0], 
    #             label='Gene Nodes', 
    #             s=SIZE)

    # # Plot microbe nodes
    # microbe_indices = [embedding_kv.index_to_key.index(node) for node in microbe_nodes]
    # plt.scatter(x[microbe_indices], y[microbe_indices], 
    #             alpha=ALPHA, 
    #             linewidth=0, 
    #             c=colors[1], 
    #             label='Microbe Nodes',
    #             s=SIZE)

    # Plot gene feature nodes
    gene_feature_indices = [embedding_kv.index_to_key.index(node) for node in gene_feature_nodes]
    plt.scatter(x[gene_feature_indices], y[gene_feature_indices], 
                alpha=ALPHA, 
                linewidth=0, 
                c=colors[2], 
                label='Gene Feature Nodes', 
                s=SIZE)

    # Plot microbe feature nodes
    microbe_feature_indices = [embedding_kv.index_to_key.index(node) for node in microbe_feature_nodes]
    plt.scatter(x[microbe_feature_indices], y[microbe_feature_indices], 
                alpha=ALPHA, 
                linewidth=0, 
                c=colors[3], 
                label='Microbe Feature Nodes', 
                s=SIZE)

    plt.grid(True)
    plt.legend()

    plt.savefig(PLOT_F, format='png', dpi=150, bbox_inches='tight')

    plt.close()


if __name__ == '__main__':
    main()

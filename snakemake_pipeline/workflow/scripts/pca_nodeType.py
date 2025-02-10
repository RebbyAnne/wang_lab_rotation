from gensim.models import KeyedVectors
from sklearn import decomposition
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

EMBEDDING_F = "node2vec/noIsolates_embeddings.embeddings"
gene_features_matrix_f = "gene_features_matrix.csv"
microbe_features_matrix_f = "microbe_features_matrix.csv"
PLOT_F = "node2vec/noIsolates_pca.png"

N_COMPONENTS = 2
ALPHA = 0.7
SIZE = 20

def main():
    embedding_kv = KeyedVectors.load(EMBEDDING_F)

    gene_features_matrix = pd.read_csv(gene_features_matrix_f, index_col=0)
    microbe_features_matrix = pd.read_csv(microbe_features_matrix_f, index_col=0)

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

    # Plot microbe nodes
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

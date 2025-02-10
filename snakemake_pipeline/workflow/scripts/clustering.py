from gensim.models import KeyedVectors
from sklearn import decomposition
from sklearn import metrics
from sklearn.cluster import DBSCAN
import pandas as pd
import numpy as np
import sys
import csv
import matplotlib.pyplot as plt
import seaborn as sns

EMBEDDING_F = sys.argv[1]
GENE_FEATURES_MATRIX_F = sys.argv[2]
MICROBE_FEATURES_MATRIX_F = sys.argv[3]
PLOT_F = sys.argv[4]

N_COMPONENTS = int(sys.argv[5])
ALPHA = float(sys.argv[6])
SIZE = int(sys.argv[7])

MIN_SAMPLES = int(sys.argv[8])
EPSILON = float(sys.argv[9])

OUTPUT_F = sys.argv[10]


EMBEDDING_F = "workflow/out/embeddings/noIsolates_100Lw100Nw0.1p0.1q128k.embeddings"
GENE_FEATURES_MATRIX_F = "workflow/out/matrices/gene_features_matrix.csv"
MICROBE_FEATURES_MATRIX_F = "workflow/out/matrices/microbe_features_matrix.csv"

N_COMPONENTS = 2

MIN_SAMPLES = 3
EPSILON = 0.04

OUTPUT_F = "workflow/out/dbscan_clusters/noIsolates_100Lw100Nw0.1p0.1q128k_clusterInfo.csv"
PLOT_F = "workflow/out/dbscan_clusters_plots/noIsolates_100Lw100Nw0.1p0.1q128k_PCA_DBSCAN.png"

def main():
    vectors = extract_feature_vectors()

    reduced_dimension_vectors = reduce_dimensions(vectors)

    db, labels, n_clusters, n_noise_points, silhouette_score = cluster(reduced_dimension_vectors)

    write_file(n_clusters, n_noise_points, silhouette_score)

    plot_data(labels, reduced_dimension_vectors)


def extract_feature_vectors():

    embedding_kv = KeyedVectors.load(EMBEDDING_F)

    gene_features_matrix = pd.read_csv(GENE_FEATURES_MATRIX_F, index_col=0)
    microbe_features_matrix = pd.read_csv(MICROBE_FEATURES_MATRIX_F, index_col=0)

    gene_features = gene_features_matrix.columns.tolist()
    microbe_features = microbe_features_matrix.columns.tolist()

    gene_feature_nodes = [node for node in embedding_kv.index_to_key if node in gene_features]
    microbe_feature_nodes = [node for node in embedding_kv.index_to_key if node in microbe_features]

    feature_nodes = gene_feature_nodes + microbe_feature_nodes

    return np.array([embedding_kv[node] for node in feature_nodes])


def reduce_dimensions(vectors):
    pca_model = decomposition.PCA(n_components=N_COMPONENTS)
    return pca_model.fit_transform(vectors)


def cluster(reduced_dimension_vectors):
    db = DBSCAN(eps=EPSILON, min_samples=MIN_SAMPLES).fit(reduced_dimension_vectors)
    labels = db.labels_

    # Number of clusters in labels, ignoring noise if present.
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise_points = list(labels).count(-1)
    silhouette_score = metrics.silhouette_score(reduced_dimension_vectors, labels)

    return db, labels, n_clusters, n_noise_points, silhouette_score


def write_file(n_clusters, n_noise_points, silhouette_score):
    data = [
        ["n_clusters", "n_noise_points", "silhouette_score"],
        [n_clusters, n_noise_points, silhouette_score]
    ]

    with open(OUTPUT_F, mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerows(data)

 
def plot_data(labels, pca_results):
    unique_labels = sorted(set(labels))  # Ensure consistent order
    n_clusters = len(unique_labels)

    # Use a seaborn palette for 50+ colors
    palette = sns.color_palette("husl", n_colors=n_clusters)

    plt.figure(figsize=(12, 8))
    for k, color in zip(unique_labels, palette):
        if k == -1:
            # Black used for noise
            color = (0.0, 0.0, 0.0)

        # Filter points by cluster
        class_member_mask = (labels == k)
        xy = pca_results[class_member_mask]

        # Scatter plot for the cluster
        plt.scatter(xy[:, 0], xy[:, 1], c=[color], label=f"Cluster {k}" if k != -1 else "Noise", s=30)

    plt.title("DBSCAN Clustering Visualization (50+ Clusters)")
    plt.xlabel("PCA Component 1")
    plt.ylabel("PCA Component 2")
    plt.savefig(PLOT_F, format="png", dpi=150, bbox_inches="tight")
    plt.close()


    
if __name__ == '__main__':
    main()





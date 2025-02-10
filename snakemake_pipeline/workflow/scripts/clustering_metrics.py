from gensim.models import KeyedVectors
from sklearn import decomposition
from sklearn import metrics
from sklearn.cluster import DBSCAN
import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt
import seaborn as sns

EMBEDDING_F = sys.argv[1]
GENE_FEATURES_MATRIX_F = sys.argv[2]
MICROBE_FEATURES_MATRIX_F = sys.argv[3]

N_COMPONENTS = int(sys.argv[4])
ALPHA = float(sys.argv[5])
SIZE = int(sys.argv[6])

MIN_SAMPLES = int(sys.argv[7])
EPSILON = float(sys.argv[8])

DIMENSIONS = str(sys.argv[9])
WALK_LEN = str(sys.argv[10])
N_WALKS = str(sys.argv[11])
PVAL = str(sys.argv[12])
QVAL = str(sys.argv[13])

OUTPUT_F = sys.argv[14]
# PLOT_F = sys.argv[14]

def main():
    vectors, feature_nodes, node_types = extract_feature_vectors()

    reduced_dimension_vectors = reduce_dimensions(vectors)

    cluster_info_df, labels = cluster(reduced_dimension_vectors,feature_nodes, node_types)

    cluster_info_df.to_csv(OUTPUT_F, index=False)

    # plot_data(labels, reduced_dimension_vectors)

def extract_feature_vectors():
    embedding_kv = KeyedVectors.load(EMBEDDING_F)

    gene_features_matrix = pd.read_csv(GENE_FEATURES_MATRIX_F, index_col=0)
    microbe_features_matrix = pd.read_csv(MICROBE_FEATURES_MATRIX_F, index_col=0)

    gene_features = gene_features_matrix.columns.tolist()
    microbe_features = microbe_features_matrix.columns.tolist()

    gene_feature_nodes = [node for node in embedding_kv.index_to_key if node in gene_features]
    microbe_feature_nodes = [node for node in embedding_kv.index_to_key if node in microbe_features]

    feature_nodes = gene_feature_nodes + microbe_feature_nodes
    node_types = {node: "gene" if node in gene_feature_nodes else "microbe" for node in feature_nodes}

    vectors = np.array([embedding_kv[node] for node in feature_nodes])

    return vectors, feature_nodes, node_types

def reduce_dimensions(vectors):
    pca_model = decomposition.PCA(n_components=N_COMPONENTS)
    return pca_model.fit_transform(vectors)

def cluster(reduced_dimension_vectors, feature_nodes, node_types):
    db = DBSCAN(eps=EPSILON, min_samples=MIN_SAMPLES).fit(reduced_dimension_vectors)
    labels = db.labels_

    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise_points = list(labels).count(-1)
    overall_silhouette_score = metrics.silhouette_score(reduced_dimension_vectors, labels) if n_clusters > 1 else None

    cluster_info = []

    # cluster-specific metrics
    sample_silhouette_values = metrics.silhouette_samples(reduced_dimension_vectors, labels)

    for cluster_id in set(labels):
        if cluster_id == -1:
            continue  # Skip noise points

        cluster_silhouette_avg = sample_silhouette_values[labels == cluster_id].mean()

        # Get indices of nodes in the current cluster
        cluster_points = reduced_dimension_vectors[labels == cluster_id]
        cluster_indices = np.where(labels == cluster_id)[0]
        cluster_nodes = [feature_nodes[i] for i in cluster_indices]

        gene_nodes = [node for node in cluster_nodes if node_types[node] == "gene"]
        microbe_nodes = [node for node in cluster_nodes if node_types[node] == "microbe"]

        #  Node type composition
        # node_type_counts = {"gene": 0, "microbe": 0}

        # for node in cluster_nodes:
        #     node_type_counts[node_types[node]] += 1
        node_type_counts = {"gene": len(gene_nodes), "microbe": len(microbe_nodes)}

        # Cluster diameter
        cluster_diameter = np.max(metrics.pairwise_distances(cluster_points)) if len(cluster_points) > 1 else 0
    
        cluster_info.append({
            "Cluster ID": cluster_id,
            "Number of Nodes": len(cluster_nodes),
            "Gene Nodes": node_type_counts["gene"],
            "Microbe Nodes": node_type_counts["microbe"],
            "Gene Node List": gene_nodes,
            "Microbe Node List": microbe_nodes,
            "Silhouette Score": cluster_silhouette_avg,
            "Cluster Diameter": cluster_diameter,
            "min_samples": MIN_SAMPLES,
            "epsilon": EPSILON,
            "dimensions": DIMENSIONS,
            "walk_length": WALK_LEN,
            "n_walks": N_WALKS,
            "p": PVAL,
            "q": QVAL
        })

    return pd.DataFrame(cluster_info), labels


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
        plt.scatter(xy[:, 0], xy[:, 1], c=[color], label=f"Cluster {k}" if k != -1 else "Noise", s=SIZE, alpha = ALPHA)

    title = "DBSCAN params: " + str(MIN_SAMPLES) + " min samples, " + str(EPSILON) + " epsilon \n"  + \
            "Node2Vec params: " + str(DIMENSIONS) + " dimensions, " + str(WALK_LEN) + " walk length, " + str(N_WALKS) + " num walks, " + str(PVAL) + " p, " + str(QVAL) + " q"

    plt.title(title)
    plt.xlabel("PCA Component 1")
    plt.ylabel("PCA Component 2")
    plt.savefig(PLOT_F, format="png", dpi=150, bbox_inches="tight")
    plt.close()

if __name__ == '__main__':
    main()

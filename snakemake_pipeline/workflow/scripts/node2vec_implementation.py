import networkx as nx
import pandas as pd
from node2vec import Node2Vec

n_dimensions = 128
n_nodes_walk = 4
n_walks = 200
# hyperparams for DFS
p = 10
q = 0.1

adjacency_matrix = pd.read_csv("positive_correlation_and_feature_adjacency_matrix_thresholdFilter.csv", index_col=0)
EMBEDDING_FILENAME = "node2vec/noIsolates_embeddings.embeddings"
EMBEDDING_MODEL_FILENAME = "node2vec/noIsolates_model.model"

graph = nx.from_pandas_adjacency(adjacency_matrix)
graph.remove_nodes_from(list(nx.isolates(graph)))

node2vec = Node2Vec(graph, dimensions=n_dimensions, 
                    walk_length=n_nodes_walk, 
                    num_walks=n_walks, 
                    p=p, 
                    q=q)

model = node2vec.fit(window=10, min_count=1, batch_words=4)

# Save embeddings for later use
model.wv.save(EMBEDDING_FILENAME)

# Save model for later use
# model.save(EMBEDDING_MODEL_FILENAME)
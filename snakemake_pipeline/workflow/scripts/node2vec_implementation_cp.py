import networkx as nx
import pandas as pd
from node2vec import Node2Vec
import sys

ADJACENCY_MATRIX_F = sys.argv[1]
EMBEDDING_FILENAME = sys.argv[2]

n_dimensions = int(sys.argv[3])
n_nodes_walk = int(sys.argv[4])
n_walks = int(sys.argv[5])
# hyperparams for DFS
p = float(sys.argv[6])
q = float(sys.argv[7])

adjacency_matrix = pd.read_csv(ADJACENCY_MATRIX_F, index_col=0)

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

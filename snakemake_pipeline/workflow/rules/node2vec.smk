rule node2vec:
    input: join(config["matricesDir"], "positive_correlation_and_feature_adjacency_matrix_thresholdFilter.csv")
    output: join(config["embeddingsDir"], "noIsolates_{walk_length}Lw{n_walks}Nw{p}p{q}q{k}k.embeddings")
    params:
        dimensions="{k}",
        walk_length="{walk_length}",
        n_walks="{n_walks}",
        p="{p}",
        q="{q}"
    shell:
        """
        python3 workflow/scripts/node2vec_implementation_cp.py \
        {input} \
        {output} \
        {params.dimensions} \
        {params.walk_length} \
        {params.n_walks} \
        {params.p} \
        {params.q}
        """

rule PCA_ByNodeType:
    input: 
        embeddings = join(config["embeddingsDir"], "noIsolates_{walk_length}Lw{n_walks}Nw{p}p{q}q{k}k.embeddings"), 
        gene_features_matrix = join(config["matricesDir"], "gene_features_matrix.csv"),
        microbe_features_matrix = join(config["matricesDir"], "microbe_features_matrix.csv")
    output: join(config["plotsDir"], "noIsolates_{walk_length}Lw{n_walks}Nw{p}p{q}q{k}k_PCA.png")
    params:
        n_components=config["n_components"], 
        alpha=config["alpha"], 
        size=config["size"], 
        dimensions="{k}",
        walk_length="{walk_length}",
        n_walks="{n_walks}",
        p="{p}",
        q="{q}"
    shell:
        """
        python3 workflow/scripts/pca_nodeType_cp.py \
        {input.embeddings} \
        {input.gene_features_matrix} \
        {input.microbe_features_matrix} \
        {output} \
        {params.n_components} \
        {params.alpha} \
        {params.size} \
        {params.dimensions} \
        {params.walk_length} \
        {params.n_walks} \
        {params.p} \
        {params.q}
        """

rule PCA_ByNodeType_featureNodes:
    input: 
        embeddings = join(config["embeddingsDir"], "noIsolates_{walk_length}Lw{n_walks}Nw{p}p{q}q{k}k.embeddings"), 
        gene_features_matrix = join(config["matricesDir"], "gene_features_matrix.csv"),
        microbe_features_matrix = join(config["matricesDir"], "microbe_features_matrix.csv")
    output: join(config["plotsDir"], "noIsolates_{walk_length}Lw{n_walks}Nw{p}p{q}q{k}k_PCA_featureNodes.png")
    params:
        n_components=config["n_components"], 
        alpha=config["alpha"], 
        size=config["size"]
    shell:
        """
        python3 workflow/scripts/pca_nodeType_cp.py \
        {input.embeddings} \
        {input.gene_features_matrix} \
        {input.microbe_features_matrix} \
        {output} \
        {params.n_components} \
        {params.alpha} \
        {params.size}
        """
     
rule pw_featureDistances:
    input: 
        embeddings = join(config["embeddingsDir"], "noIsolates_{walk_length}Lw{n_walks}Nw{p}p{q}q{k}k.embeddings"), 
        gene_features_matrix = join(config["matricesDir"], "gene_features_matrix.csv"),
        microbe_features_matrix = join(config["matricesDir"], "microbe_features_matrix.csv")
    output: 
        dist_matrix=join(config["distMatrix"], "noIsolates_{walk_length}Lw{n_walks}Nw{p}p{q}q{k}k_featureDistanceMatrix.csv"),
        dist_csv=join(config["distDF"], "noIsolates_{walk_length}Lw{n_walks}Nw{p}p{q}q{k}k_featureDistanceData_long.csv")
    params:        
        n_components=config["n_components"], 
        dimensions="{k}",
        walk_length="{walk_length}",
        n_walks="{n_walks}",
        p="{p}",
        q="{q}"
    shell:
        """
        python3 workflow/scripts/pw_redDim_distance_feature_nodes.py \
        {input.embeddings} \
        {input.gene_features_matrix} \
        {input.microbe_features_matrix} \
        {output.dist_matrix} \
        {output.dist_csv} \
        {params.dimensions} \
        {params.walk_length} \
        {params.n_walks} \
        {params.p} \
        {params.q}
        """

rule cluster_evaluate:
    input: 
        embeddings = join(config["embeddingsDir"], "noIsolates_{walk_length}Lw{n_walks}Nw{p}p{q}q{k}k.embeddings"), 
        gene_features_matrix = join(config["matricesDir"], "gene_features_matrix.csv"),
        microbe_features_matrix = join(config["matricesDir"], "microbe_features_matrix.csv")
    output:
        join(config["clusteringMetrics"], "noIsolates_{walk_length}Lw{n_walks}Nw{p}p{q}q{k}k_{min_sample}ms{epsilon}eps_DBSCAN_cluster_info.csv")
    params:
        n_components=config["n_components"], 
        alpha=config["alpha"], 
        size=config["size"], 
        min_sample="{min_sample}", 
        epsilon="{epsilon}",
        dimensions="{k}",
        walk_length="{walk_length}",
        n_walks="{n_walks}",
        p="{p}",
        q="{q}"
    shell:
        """
        python3 workflow/scripts/clustering_metrics.py \
        {input.embeddings} \
        {input.gene_features_matrix} \
        {input.microbe_features_matrix} \
        {params.n_components} \
        {params.alpha} \
        {params.size} \
        {params.min_sample} \
        {params.epsilon} \
        {params.dimensions} \
        {params.walk_length} \
        {params.n_walks} \
        {params.p} \
        {params.q} \
        {output}
        """

rule cluster_viz:
    input: 
        embeddings = join(config["embeddingsDir"], "noIsolates_{walk_length}Lw{n_walks}Nw{p}p{q}q{k}k.embeddings"), 
        gene_features_matrix = join(config["matricesDir"], "gene_features_matrix.csv"),
        microbe_features_matrix = join(config["matricesDir"], "microbe_features_matrix.csv")
    output:
        join(config["clusterPlots"], "noIsolates_{walk_length}Lw{n_walks}Nw{p}p{q}q{k}k_{min_sample}ms{epsilon}eps_PCA_DBSCAN.png")
    params:
        n_components=config["n_components"], 
        alpha=config["alpha"], 
        size=config["size"], 
        min_sample="{min_sample}", 
        epsilon="{epsilon}",
        dimensions="{k}",
        walk_length="{walk_length}",
        n_walks="{n_walks}",
        p="{p}",
        q="{q}"
    shell:
        """
        python3 workflow/scripts/clustering_viz.py \
        {input.embeddings} \
        {input.gene_features_matrix} \
        {input.microbe_features_matrix} \
        {params.n_components} \
        {params.alpha} \
        {params.size} \
        {params.min_sample} \
        {params.epsilon} \
        {params.dimensions} \
        {params.walk_length} \
        {params.n_walks} \
        {params.p} \
        {params.q} \
        {output}
        """


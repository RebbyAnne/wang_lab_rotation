import pandas as pd
import itertools
import numpy as np

correlation_matrix_f = "../out/matrices/microbe_gene_correlation_matrix.csv"
pvalue_matrix_f = "../out/matrices/microbe_gene_pval_matrix.csv"
gene_features_matrix_f = "../out/matrices/gene_features_matrix.csv"
microbe_features_matrix_f = "../out/matrices/microbe_features_matrix.csv"

corr_threshold = 0
pval_threshold = 0.05

def main(): 
    correlation_matrix, pvalue_matrix, gene_features_matrix, microbe_features_matrix = read_dfs()

    correlation_matrix, pvalue_matrix, gene_features_matrix, microbe_features_matrix = get_common_indeces(correlation_matrix, pvalue_matrix, gene_features_matrix, microbe_features_matrix)

    threshold_correlation_matrix = threshold_corr_matrix(correlation_matrix, pvalue_matrix)

    # adjacency_matrix = merge_matrices(threshold_correlation_matrix, gene_features_matrix, microbe_features_matrix)

    feature_correlation_matrix(threshold_correlation_matrix, microbe_features_matrix, gene_features_matrix)

    # print(adjacency_matrix)

    # pos_adjacency_matrix = adjacency_matrix.clip(lower=0)
    # neg_adjacency_matrix = adjacency_matrix.clip(upper=0)

    # adjacency_matrix.to_csv("correlation_and_feature_adjacency_matrix_thresholdFilter.csv")
    # pos_adjacency_matrix.to_csv("positive_correlation_and_feature_adjacency_matrix_thresholdFilter.csv")

def read_dfs():
    correlation_matrix = pd.read_csv(correlation_matrix_f, index_col=0)
    pvalue_matrix = pd.read_csv(pvalue_matrix_f, index_col=0)
    gene_features_matrix = pd.read_csv(gene_features_matrix_f, index_col=0)
    microbe_features_matrix = pd.read_csv(microbe_features_matrix_f, index_col=0) 

    gene_features_matrix = gene_features_matrix.astype("float")
    microbe_features_matrix = microbe_features_matrix.astype("float")

    return correlation_matrix, pvalue_matrix, gene_features_matrix, microbe_features_matrix

def get_common_indeces(correlation_matrix, pvalue_matrix, gene_features_matrix, microbe_features_matrix):
    common_microbes = correlation_matrix.index.intersection(microbe_features_matrix.index)
    common_genes = correlation_matrix.columns.intersection(gene_features_matrix.index)

    combined = common_microbes.tolist() + common_genes.tolist()

    common_microbe_feature_matrix = microbe_features_matrix.loc[common_microbes]
    common_gene_feature_matrix = gene_features_matrix.loc[common_genes]
    common_corr_matrix = correlation_matrix.loc[combined, combined]
    common_pval_matrix = pvalue_matrix.loc[combined, combined]

    # set microbe-microbe and gene-gene pairs to 0
    common_corr_matrix.loc[common_microbes, common_microbes] = 0
    common_corr_matrix.loc[common_genes, common_genes] = 0

    return common_corr_matrix, common_pval_matrix, common_gene_feature_matrix, common_microbe_feature_matrix

def threshold_corr_matrix(correlation_matrix, pvalue_matrix):
    correlation_matrix[pvalue_matrix > pval_threshold] = 0
    correlation_matrix[correlation_matrix < corr_threshold] = 0

    return correlation_matrix


def merge_matrices(correlation_matrix, gene_features_matrix, microbe_features_matrix):
    correlation_with_gene_features = correlation_matrix.merge(
        gene_features_matrix, 
        left_index=True, 
        right_index=True, 
        how="left"
    ).fillna(0)

    correlation_with_gene_features_sq = correlation_with_gene_features.T.merge(
        gene_features_matrix, 
        left_index=True, 
        right_index=True, 
        how="left"
    ).fillna(0)

    final_adjacency_matrix = correlation_with_gene_features_sq.merge(
        microbe_features_matrix, 
        left_index=True, 
        right_index=True, 
        how="left"
    ).fillna(0)

    final_adjacency_matrix_sq = final_adjacency_matrix.T.merge(
        microbe_features_matrix, 
        left_index=True, 
        right_index=True, 
        how="left"
    ).fillna(0)

    return final_adjacency_matrix_sq.loc[final_adjacency_matrix_sq.index, final_adjacency_matrix_sq.index]

def feature_correlation_matrix(corr_matrix, microbe_feature_matrix, gene_feature_matrix):
    print(corr_matrix)
    print(corr_matrix.shape)

    print(microbe_feature_matrix)
    print(microbe_feature_matrix.shape)

    print(gene_feature_matrix)
    print(gene_feature_matrix.shape)

    # common_microbes = corr_matrix.index.intersection(microbe_feature_matrix.index)
    # common_genes = corr_matrix.columns.intersection(gene_feature_matrix.index)

    # common_microbe_feature_matrix = microbe_feature_matrix.loc[common_microbes]
    # common_gene_feature_matrix = gene_feature_matrix.loc[common_genes]
    # common_corr_matrix = corr_matrix.loc[common_microbes, common_genes]

    # Aggregate correlations for microbe features (transposed)
    # microbe_feature_correlation =  microbe_feature_matrix.T.abs() @ corr_matrix.abs()

    # Aggregate correlations for gene features
    # feature_correlation_matrix = microbe_feature_correlation @ common_gene_feature_matrix.abs()

    # return feature_correlation_matrix


if __name__ == "__main__":
    main()


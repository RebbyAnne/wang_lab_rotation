import pandas as pd
import scipy.stats as stats

directory = "../20240325_ExfoSeq_SPARC_IBD_Healthy"
# HOST data
sample_metadata_host = pd.read_csv(directory + "/sample.metadata.tsv", sep="\t")
gene_metadata = pd.read_csv(directory + "/gene.metadata.tsv", sep="\t")
housekeeping_genes = pd.read_csv(directory +  "/geneSymbol.house_keep.list", header=None)

data_count_P1_2 = pd.read_csv(directory + "/archive/output_merge.P1_P2.tsv", sep="\t")
data_count_P3 = pd.read_csv(directory + "/archive/output_merge.P3.tsv", sep="\t")

data_stat_P1_2 = pd.read_csv(directory + "/archive/output_merge.P1_P2.stat", sep="\t")
data_stat_P3 = pd.read_csv(directory + "/archive/output_merge.P3.stat", sep="\t")

# MICROBE data
sample_metadata_microbe = pd.read_csv(directory + "/archive/archive_20231126_16SV4_SPARC_IBD_DIVERS_H1_H7_H72/sample.metadata.tsv", sep="\t")
otu_taxonomy = pd.read_csv(directory + "/archive/archive_20231126_16SV4_SPARC_IBD_DIVERS_H1_H7_H72/rdp.taxonomy.tsv", sep="\t", header=None)
otu_table = pd.read_csv(directory + "/archive/archive_20231126_16SV4_SPARC_IBD_DIVERS_H1_H7_H72/otu_table_97.txt", sep="\t")

column_names = ['otu', 'genus', 'root', 'kingdom', 'phylum', 'class', 'order', 'family']
taxonomy_table = pd.read_csv(directory + "/archive/archive_20231126_16SV4_SPARC_IBD_DIVERS_H1_H7_H72/rdp.taxonomy.tsv", sep="\t", 
                             header=None, names=column_names)

kegg_metadata = pd.read_csv(directory + "/kegg_pathway_metadata.csv")

def main():    
    host_cpm_data = normalize_host_genes(housekeeping_genes)

    microbe_cpm_data = normalize_microbe_otu(otu_table, sample_metadata_microbe)

    # corr_matrix = get_correlation_matrix(host_cpm_data, microbe_cpm_data) 

    corr_matrix, pval_matrix = get_correlation_pval_matrix(host_cpm_data, microbe_cpm_data)

    corr_matrix.to_csv("microbe_gene_correlation_matrix.csv", index=True)
    pval_matrix.to_csv("microbe_gene_pval_matrix.csv", index=True) 

    microbe_feature_matrix = get_microbe_feature_matrix(taxonomy_table)
    microbe_feature_matrix.to_csv("microbe_features_matrix.csv", index=True)

    gene_feature_matrix = get_hostGene_feature_matrix(kegg_metadata)
    gene_feature_matrix.to_csv("gene_features_matrix.csv", index=True)

    feature_corr_matrix = gene_microbe_feature_correlation_matrix(corr_matrix,
                                                                  microbe_feature_matrix,
                                                                  gene_feature_matrix
    )

    feature_corr_matrix.to_csv("feature_corr_matrix.csv", index=True)

    normalized_matrix = normalize_matrix(feature_corr_matrix, microbe_feature_matrix, gene_feature_matrix)

    normalized_matrix.to_csv("microbe_gene_count_normalized_matrix.csv")

def normalize_host_genes(housekeeping_genes):

    housekeeping_genes_list = housekeeping_genes[0].tolist()
        # Combine gene counts and read stats
    gene_counts = pd.concat([data_count_P1_2, data_count_P3], ignore_index=True)
    read_stats = pd.concat([data_stat_P1_2, data_stat_P3], ignore_index=True)

    # Housekeeping gene counts
    hk_gene_counts = (
        gene_counts.merge(gene_metadata, on="gene")
        .query("geneSymbol in @housekeeping_genes_list")
        .groupby("sample")
        .agg(sum_umi_hkgenes=("num_of_umi", "sum"))
        .reset_index()
    )

    # Join gene metadata and compute CPM
    host_cpm_data = (
        gene_counts.merge(gene_metadata, on="gene")
        .merge(read_stats, on="sample")
        .merge(hk_gene_counts, on="sample")
        .assign(
            rna_cpm=lambda df: (df["num_of_umi"] / df["unique_UMI"]) * 1e6,
            cpm_hk_norm=lambda df: (df["num_of_umi"] / df["sum_umi_hkgenes"]) * 1e6
        )
        .merge(
            sample_metadata_host[["label", "patient", "disease_type"]],
            left_on="sample",
            right_on="label"
        )
        .groupby(["geneSymbol", "patient"])
        .agg(mean_hk_norm_cpm=("cpm_hk_norm", "mean"))
        .reset_index()
        .pivot(index="patient", columns="geneSymbol", values="mean_hk_norm_cpm")
        .reset_index()
    )

    return(host_cpm_data)

def normalize_microbe_otu(otu_table, sample_metadata_microbe):
    # Compute total reads
    microbe_total_counts = (
        otu_table.drop(columns=["OTU"])
        .sum(axis=0)
        .reset_index(name="total_16S_reads")
        .rename(columns={"index": "label"})
    )

    sample_metadata_microbe = sample_metadata_microbe.rename(columns={"patient_ID": "patient"})
    sample_metadata_microbe['patient'] = sample_metadata_microbe['patient'].astype(str)

    # Compute CPM
    microbe_cpm_data = (
        otu_table.loc[:, otu_table.columns.str.startswith("rRNA") | (otu_table.columns == "OTU")]
        .melt(id_vars=["OTU"], var_name="label", value_name="num_of_reads")
        .merge(sample_metadata_microbe, on="label")
        .merge(microbe_total_counts, on="label")
        .assign(
            otu_cpm=lambda df: (df["num_of_reads"] / df["total_16S_reads"]) * 1e6
        )
        .merge(
            sample_metadata_host[["patient", "disease_type"]],
            on="patient"
        )
        .groupby(["OTU", "patient"])
        .agg(mean_otu_cpm=("otu_cpm", "mean"))
        .reset_index()
        .pivot(index="patient", columns="OTU", values="mean_otu_cpm")
        .reset_index()
    )

    return(microbe_cpm_data)


def get_correlation_matrix(microbe_counts, host_counts):
    host_microbe_merged = pd.merge(microbe_counts, host_counts, on="patient").drop(columns=["patient"])

    return host_microbe_merged.corr(method="pearson")


def get_correlation_pval_matrix(microbe_counts, host_counts):
    host_microbe_merged = pd.merge(microbe_counts, host_counts, on="patient").drop(columns=["patient"])
    host_microbe_merged.fillna(0, inplace=True)

    # Initialize empty matrices for correlation and p-values
    corr_matrix = pd.DataFrame(index=host_microbe_merged.columns, columns=host_microbe_merged.columns)
    pval_matrix = pd.DataFrame(index=host_microbe_merged.columns, columns=host_microbe_merged.columns)

    # Calculate the correlation and p-value for each pair of columns
    for col1 in host_microbe_merged.columns:
        for col2 in host_microbe_merged.columns:
            corr, pval = stats.pearsonr(host_microbe_merged[col1], host_microbe_merged[col2])
            corr_matrix.at[col1, col2] = corr
            pval_matrix.at[col1, col2] = pval

    # Convert the pval_matrix to numeric, to handle any issues with string types
    corr_matrix.fillna(0, inplace=True) 
    pval_matrix.fillna(0, inplace=True) 

    pval_matrix = pval_matrix.apply(pd.to_numeric)

    return corr_matrix, pval_matrix


def get_microbe_feature_matrix(tax_table):
    # Melt the dataframe so that we have OTUs and taxonomy features
    melted_df = tax_table.melt(id_vars=['otu'], 
                               value_vars=['genus', 'phylum', 'class', 'order', 'family'],
                               var_name='feature_type', 
                               value_name='feature')


    # Create a binary matrix (1 if OTU belongs to a taxonomy feature, 0 otherwise)
    return pd.crosstab(melted_df['otu'], melted_df['feature'])


def get_hostGene_feature_matrix(kegg_metadata):
    kegg_metadata = kegg_metadata.drop(columns=["kegg_pathway_id", "kegg_orthology_id"])
    melted_df = kegg_metadata.melt(id_vars='gene_symbol',
                                   value_vars=['kegg_pathway_name', 'kegg_orthology_name'],
                                   var_name='feature_type',
                                   value_name='feature')
    
    melted_df['feature'] = melted_df['feature'].str.split('|')  # Split multi-valued features into lists
    melted_df = melted_df.explode('feature')  # Expand lists into separate rows

    # Step 3: Create a binary matrix (gene_symbol x unique features)
    return pd.crosstab(melted_df['gene_symbol'], melted_df['feature'])


def gene_microbe_feature_correlation_matrix(corr_matrix, 
                                            microbe_feature_matrix, 
                                            gene_feature_matrix):
    """
    Calculate the feature-feature correlation matrix by summing the absolute
    values of correlations between all microbes and genes associated with features.

    Args:
    - correlation_matrix (DataFrame): Microbe x Gene correlation matrix.
    - microbe_features (DataFrame): Microbe x Microbe Features binary matrix.
    - gene_features (DataFrame): Gene x Gene Features binary matrix.

    Returns:
    - feature_correlation_matrix (DataFrame): Microbe Features x Gene Features matrix.
    """

    common_microbes = corr_matrix.index.intersection(microbe_feature_matrix.index)
    common_genes = corr_matrix.columns.intersection(gene_feature_matrix.index)

    common_microbe_feature_matrix = microbe_feature_matrix.loc[common_microbes]
    common_gene_feature_matrix = gene_feature_matrix.loc[common_genes]
    common_corr_matrix = corr_matrix.loc[common_microbes, common_genes]

    # Aggregate correlations for microbe features (transposed)
    microbe_feature_correlation =  common_microbe_feature_matrix.T.abs() @ common_corr_matrix.abs()

    # Aggregate correlations for gene features
    feature_correlation_matrix = microbe_feature_correlation @ common_gene_feature_matrix.abs()

    return feature_correlation_matrix

def normalize_matrix(correlation_matrix, microbe_feature_matrix, gene_feature_matrix):
    microbe_feature_counts = microbe_feature_matrix.sum(axis=0)  # Column sums
    gene_feature_counts = gene_feature_matrix.sum(axis=0)

    gene_count_normalized_matrix = correlation_matrix / gene_feature_counts
    normalized_matrix = gene_count_normalized_matrix.T / microbe_feature_counts
    
    return normalized_matrix.T


if __name__ == "__main__":
    main()


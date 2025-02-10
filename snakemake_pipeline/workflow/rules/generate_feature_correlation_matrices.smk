rule make_corr_matrix:
    input:
        sample_metadata_host = config["sample_metadata_host"], 
        gene_metadata = config["gene_metadata"], 
        housekeeping_genes = config["housekeeping_genes"], 
        gene_count_data_1_2 = config["gene_count_1_2"], 
        gene_count_data_3 = config["gene_count_3"], 
        gene_stats_1_2 = config["gene_stats_1_2"], 
        gene_stats_3 = config["gene_stats_3"],
        sample_metadata_microbe = config["sample_metadata_microbe"], 
        otu_tax = config["otu_tax"], 
        otu_table = config["otu_table"],
        tax_table = config["tax_table"],
        kegg_metadata = config["kegg_metadata"]
    output:

    shell:
        """
        python3 make_corr_matrix_cp.py \
        {input.sample_metadata_host} \
        {input.gene_metadata} \
        {input.}
        """


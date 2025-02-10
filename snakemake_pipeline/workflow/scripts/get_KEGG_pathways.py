import pandas as pd
from Bio.KEGG import REST

directory = "../20240325_ExfoSeq_SPARC_IBD_Healthy"

gene_metadata = pd.read_csv(directory + "/gene.metadata.tsv", sep="\t")

def main():
    gene_metadata_nodups = gene_metadata.drop_duplicates(subset='geneSymbol')

    genes = gene_metadata_nodups['geneSymbol'].tolist()

    gene_pathways = {}

    for gene in genes[1]:
        print(gene)
        gene_pathways[gene] =  get_kegg_pathways(gene)
    

def get_kegg_pathways(gene):
    try:
        # Fetch KEGG data for the gene
        kegg_data = REST.kegg_get(gene).read()

        print(kegg_data)
        
        # Find all pathways associated with the gene
        pathways = []
        for line in kegg_data.splitlines():
            if line.startswith("PATHWAY"):
                pathway = line.split("\t")[1]
                pathways.append(pathway)
        return pathways
    
    except Exception as e:
        print(f"Error retrieving data for {gene}: {e}")
        return []


if __name__ == "__main__":
    main()

library(biomaRt)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Define the list of GO terms
go_terms <- c("GO:0002376", "GO:0019221", "GO:0060333", "GO:0070102", "GO:0007259", "GO:0006955")

# Retrieve genes related to the GO terms
genes_df <- getBM(
  attributes = c("ensembl_transcript_id","ensembl_gene_id", "external_gene_name", "go_id", "name_1006"),
  filters = "go_parent_term",
  values = go_terms,
  mart = ensembl
)

#unique(genes_df$external_gene_name)

write.table(genes_df, file = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\found_SNPs\\immune_genes.tsv", quote = FALSE, sep = "\t", row.names = FALSE)


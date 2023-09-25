library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)

snps = read.csv(sep = "\t", comment.char = "#", header = FALSE, file="C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\found_SNPs\\PBMC_TCell\\SNPs_in_highAcetylation_AND_immuneGenes_AND_STATsignals_AND_IFNorILinducedFC50.vcf")

unique_values <- c()

# Loop through each row in the data.frame
for (i in 1:nrow(snps)) {
  values <- unlist(strsplit(as.character(snps$V10[i]), split = ";"))
  # Update the vector with new unique values
  unique_values <- unique(c(unique_values, values))
}

lung_expression = read.csv(sep = "\t", comment.char = "#", file="C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\RNA-seq_Data\\lung\\normalized_counts_all.tsv")
kidney_expression = read.csv(sep = "\t", comment.char = "#", file="C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\RNA-seq_Data\\kidney\\normalized_counts_all.tsv")
tcell_expression = read.csv(sep = "\t", comment.char = "#", file="C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\RNA-seq_Data\\TCells\\normalized_counts_all.tsv")

lung_expression <- lung_expression[lung_expression$hgnc_symbol %in% unique_values,]
kidney_expression <- kidney_expression[kidney_expression$hgnc_symbol %in% unique_values,]
tcell_expression <- tcell_expression[tcell_expression$hgnc_symbol %in% unique_values,]

lung_expression_long <- lung_expression %>%
  gather(key = "condition", value = "value", -hgnc_symbol)
colnames(lung_expression_long)<- c('gene','condition','gene_counts')
lung_expression_long$condition <- gsub("_rep\\d*", "", lung_expression_long$condition)

kidney_expression_long <- kidney_expression %>%
  gather(key = "condition", value = "value", -hgnc_symbol)
colnames(kidney_expression_long)<- c('gene','condition','gene_counts')
kidney_expression_long$condition <- gsub("ID\\d*_", "", kidney_expression_long$condition)


tcell_expression_long <- tcell_expression %>%
  gather(key = "condition", value = "value", -hgnc_symbol)
colnames(tcell_expression_long)<- c('gene','condition','gene_counts')
tcell_expression_long$condition <- str_match(tcell_expression_long$condition, "^[^_]+_[^_]+_([^_]+)")[,2]

palette_12 <- c("#D32F2F", "#1976D2", "#388E3C", "#FBC02D", "#8E24AA", "#F57C00", 
                "#0288D1", "#7B1FA2", "#C2185B", "#7B1FA2", "#C2185B", "#0288D1")

ggplot(lung_expression_long, aes(x = gene, y = log10(gene_counts), fill = condition)) +
  geom_boxplot() + scale_fill_manual(values = palette_12)

ggplot(kidney_expression_long, aes(x = gene, y = log2(gene_counts), fill = condition)) +
  geom_boxplot() + scale_fill_manual(values = palette_12)

ggplot(tcell_expression_long, aes(x = gene, y = log2(gene_counts), fill = condition)) +
  geom_boxplot() + scale_fill_manual(values = palette_12)


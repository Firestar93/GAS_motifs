library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)

#lung
interesting_genes_lung = read.csv(file = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\RNA-seq_Data\\lung\\combined.tsv", sep = "\t")

interesting_genes_lung <- interesting_genes_lung %>% 
  select(-starts_with("IFNb_"))

interesting_genes_lung <- interesting_genes_lung %>% 
  select(-starts_with("SAEC_GH"))

interesting_genes_lung <- interesting_genes_lung %>% 
  select(-starts_with("SAEC_Bar"))

interesting_genes_lung <- interesting_genes_lung %>% 
  select(-starts_with("SAEC_Ruxo"))

interesting_genes_lung_long <- interesting_genes_lung %>%
  gather(key = "condition", value = "value", -X)
colnames(interesting_genes_lung_long)<- c('gene','condition','gene_counts')

interesting_genes_lung_long$condition <- gsub("_rep\\d*", "", interesting_genes_lung_long$condition)

palette_12 <- c("#D32F2F", "#1976D2", "#388E3C", "#FBC02D", "#8E24AA", "#F57C00", 
                "#0288D1", "#7B1FA2", "#C2185B", "#7B1FA2", "#C2185B", "#0288D1")


ggplot(interesting_genes_lung_long, aes(x = gene, y = log(gene_counts), fill = condition)) +
  geom_boxplot() + scale_fill_manual(values = palette_12)


#kidney
kidney_expression = read.csv(file = 'C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\RNA-seq_Data\\kidney\\salmon.merged.gene_counts.tsv', sep="\t")
kidney_expression <- select(kidney_expression,-gene_id)

kidney_expression <- kidney_expression %>%
  gather(key = "condition", value = "value", -gene_name)
colnames(kidney_expression)<- c('gene','condition','gene_counts')

kidney_expression$gene_counts <- sub("\\..*$", "", kidney_expression$gene_counts)

kidney_expression$condition <- gsub("ID\\d*_", "", kidney_expression$condition)

kidney_expression$gene_counts<-as.integer(kidney_expression$gene_counts)

value_list<- c("ATF4","BHLH40-AS1","INPP5K","ITPR1","MIEF1","MYO1C","RNF144A","RSAD2")

kidney_expression <- kidney_expression %>% 
  filter(gene %in% value_list)


palette_12 <- c("#D32F2F", "#1976D2", "#388E3C", "#FBC02D", "#8E24AA", "#F57C00", 
                "#0288D1", "#7B1FA2", "#C2185B", "#7B1FA2", "#C2185B", "#0288D1")


ggplot(kidney_expression, aes(x = gene, y = gene_counts, fill = condition)) +
  geom_boxplot() + scale_fill_manual(values = palette_12)

#Tcells

Tcell_expression = read.csv(file = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\RNA-seq_Data\\TCells\\GSE189997_genecount_table.tsv", sep = "\t")
Tcell_expression <- select(Tcell_expression,c(-gene_id, -biotype, - chromosome, -region, -Exons, -Exon_length))

Tcell_expression <- Tcell_expression %>%
  gather(key = "condition", value = "value", -name)
colnames(Tcell_expression)<- c('gene','condition','gene_counts')

value_list<- c("ATF4","BHLH40-AS1","INPP5K","ITPR1","MIEF1","MYO1C","RNF144A","RSAD2")
value_list<- c("RSAD2","RNF144A","GRASLND")

Tcell_expression <- Tcell_expression %>% 
  filter(gene %in% value_list)

Tcell_expression$condition <- gsub("_HD_GU\\d*", "", Tcell_expression$condition)
Tcell_expression$condition <- gsub("_HIV_GU\\d*", "", Tcell_expression$condition)

palette_16 <- hcl.colors(16, "Set3")

ggplot(Tcell_expression, aes(x = gene, y = gene_counts, fill = condition)) +
  geom_boxplot() + scale_fill_manual(values = palette_16)

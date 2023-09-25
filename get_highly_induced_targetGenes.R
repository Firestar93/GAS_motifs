library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)
library(DESeq2)

log2fc_cutoff = log2(50)

#lung
lung_gene_expression = read.csv(file = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\RNA-seq_Data\\lung\\combined.tsv", sep = "\t")

rownames(lung_gene_expression)<-lung_gene_expression$hgnc_symbol
lung_gene_expression$hgnc_symbol<-NULL

# # Create the DESeqDataSet from the original complete data frame
# all_columns <- colnames(lung_gene_expression)
# conditions_all <- ifelse(grepl("control", all_columns, ignore.case = TRUE), "Control", "Treatment")
# colData_all <- data.frame(row.names = all_columns, condition = as.factor(conditions_all))
# 
# dds_all <- DESeqDataSetFromMatrix(countData = lung_gene_expression,
#                                   colData = colData_all,
#                                   design = ~ condition)
# 
# # Run DESeq analysis
# dds_all <- DESeq(dds_all)
# 
# # Extract normalized counts
# normalized_counts_all <- counts(dds_all, normalized=TRUE)
# # Add a column for row names and set its name to 'hgnc_symbol'
# normalized_counts_all_with_rowname <- data.frame(hgnc_symbol = rownames(normalized_counts_all), normalized_counts_all)
# # Write this new data frame to a tab-separated file
# write.table(normalized_counts_all_with_rowname, file = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\RNA-seq_Data\\lung\\normalized_counts_all.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

control_columns <- grep(".*Control.*", colnames(lung_gene_expression), value = TRUE)

# Extract columns that do not contain the word 'control' in their names
non_control_columns <- setdiff(colnames(lung_gene_expression), control_columns)

# Extract unique groups by removing replicate information (assuming format is 'GroupName_rep#')
unique_groups <- unique(gsub("_rep[0-9]+$", "", non_control_columns))

# Initialize a list to hold the new data frames
list_of_grouped_dfs <- list()

# Iterate over unique non-control groups to create new data frames
for (group in unique_groups) {
  # Find all columns that belong to this group
  group_columns <- grep(group, non_control_columns, value = TRUE)
  
  # Combine control columns with this group's columns
  new_df_group <- lung_gene_expression[, c(control_columns, group_columns)]
  
  # Store new data frame in list
  list_of_grouped_dfs[[group]] <- new_df_group
}


final_df <- data.frame()

# Iterate through each grouped data frame
for (group in names(list_of_grouped_dfs)) {
  
  # Retrieve the data frame for the group
  df_group <- list_of_grouped_dfs[[group]]
  
  # Generate colData for this specific data frame
  all_columns_group <- colnames(df_group)
  conditions_group <- ifelse(grepl("control", all_columns_group, ignore.case = TRUE), "Control", group)
  colData_group <- data.frame(row.names = all_columns_group, condition = as.factor(conditions_group))
  
  # Create the DESeqDataSet
  dds <- DESeqDataSetFromMatrix(countData = df_group,
                                colData = colData_group,
                                design = ~ condition)
  
  # Run DESeq analysis
  dds <- DESeq(dds)
  
  # Extract results
  res <- results(dds)
  
  # Sort by p-value
  res <- res[order(res$pvalue),]
  
  # Extract genes with log2 fold-change higher than X or lower than -X
  significant_genes <- subset(res, abs(log2FoldChange) > log2fc_cutoff)
  
  # Extract log2 fold change values
  lfc_values <- significant_genes$log2FoldChange
  
  # Add these to the final data frame
  if (is.null(final_df)) {
    final_df <- data.frame(row.names = rownames(significant_genes))
  }
  SAEC_group <- paste0("SAEC_", group)
  final_df[[SAEC_group]] <- NULL
  final_df[rownames(significant_genes), SAEC_group] <- lfc_values
  
}

rm(list = setdiff(ls(), c("final_df","log2fc_cutoff")))

##START WITH KIDNEY NOW
kidney_gene_expression <- read.csv("C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\RNA-seq_Data\\kidney\\salmon.merged.gene_counts.tsv", sep="\t")
kidney_gene_expression$gene_id<-NULL
# Identify duplicated 'gene_name'
duplicated_genes <- duplicated(kidney_gene_expression$gene_name) | duplicated(kidney_gene_expression$gene_name, fromLast = TRUE)
# Subset data frame to only include duplicated 'gene_name'
duplicated_df <- kidney_gene_expression[duplicated_genes, ]
# Aggregate duplicated rows by 'gene_name', summing the other columns
aggregated_df <- aggregate(. ~ gene_name, data = duplicated_df, sum)
# Remove duplicated rows from the original data frame
kidney_gene_expression <- kidney_gene_expression[!duplicated_genes, ]
# Combine the aggregated data back to the original data frame
kidney_gene_expression <- rbind(kidney_gene_expression, aggregated_df)
# Set the row names of 'final_df' based on 'gene_name'
rownames(kidney_gene_expression) <- kidney_gene_expression$gene_name
kidney_gene_expression$gene_name<-NULL

rm(list = setdiff(ls(), c("final_df","kidney_gene_expression","log2fc_cutoff")))

# Identify the control and non-control columns
control_columns <- grep("untreated", colnames(kidney_gene_expression), value = TRUE)
non_control_columns <- setdiff(colnames(kidney_gene_expression), control_columns)

# Extract unique groups by removing the 'ID' followed by numbers and underscore
unique_groups <- unique(gsub("^ID[0-9]+_", "", non_control_columns))

# Loop over each unique group to perform DESeq2 analysis
for (group in unique_groups) {
  
  # Columns corresponding to the current group
  group_columns <- grep(paste0("^ID[0-9]+_", group, "$"), colnames(kidney_gene_expression), value = TRUE)
  
  # Subset DataFrame to include only the control and the current group
  sub_df <- kidney_gene_expression[, c(control_columns, group_columns)]
  
  # Create colData DataFrame
  col_names <- colnames(sub_df)
  conditions <- ifelse(col_names %in% control_columns, "Control", "Treatment")
  colData <- data.frame(row.names = col_names, condition = as.factor(conditions))
  
  sub_df <- floor(sub_df)
  
  # Create DESeqDataSet and run DESeq2
  dds <- DESeqDataSetFromMatrix(countData = sub_df, colData = colData, design = ~ condition)
  dds <- DESeq(dds)
  
  # Extract results and filter by log2 fold change
  res <- results(dds)
  significant_genes <- subset(res, abs(log2FoldChange) > log2fc_cutoff)
  
  lfc_values <- significant_genes$log2FoldChange
  
  # Add these to the final data frame
  if (is.null(final_df)) {
    final_df <- data.frame(row.names = rownames(significant_genes))
  }
  kidney_group <- paste0("KIDNEY_", group)
  final_df[[kidney_group]] <- NULL
  final_df[rownames(significant_genes), kidney_group] <- lfc_values
}

# normalized_counts_all <- counts(dds, normalized=TRUE)
# # Add a column for row names and set its name to 'hgnc_symbol'
# normalized_counts_all_with_rowname <- data.frame(hgnc_symbol = rownames(normalized_counts_all), normalized_counts_all)
# # Write this new data frame to a tab-separated file
# write.table(normalized_counts_all_with_rowname, file = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\RNA-seq_Data\\kidney\\normalized_counts_all.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

rm(list = setdiff(ls(), c("final_df","log2fc_cutoff")))

##TCELLS
TCell_gene_expression <- read.csv("C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\RNA-seq_Data\\TCells\\GSE189997_genecount_table.tsv", sep="\t")
TCell_gene_expression$gene_id<-NULL
TCell_gene_expression$biotype<-NULL
TCell_gene_expression$chromosome<-NULL
TCell_gene_expression$region<-NULL
TCell_gene_expression$Exon_length<-NULL
TCell_gene_expression$Exons<-NULL
rownames(TCell_gene_expression)<-TCell_gene_expression$name
TCell_gene_expression$name<-NULL

# Create the DESeqDataSet from the original complete data frame
all_columns <- colnames(TCell_gene_expression)
conditions_all <- ifelse(grepl("Unstim", all_columns, ignore.case = TRUE), "Unstim", "Treatment")
colData_all <- data.frame(row.names = all_columns, condition = as.factor(conditions_all))

# TCell_gene_expression[is.na(TCell_gene_expression)] <- 0
# 
# dds_all <- DESeqDataSetFromMatrix(countData = TCell_gene_expression,
#                                   colData = colData_all,
#                                   design = ~ condition)
# 
# # Run DESeq analysis
# dds_all <- DESeq(dds_all)
# 
# # Extract normalized counts
# normalized_counts_all <- counts(dds_all, normalized=TRUE)
# # Add a column for row names and set its name to 'hgnc_symbol'
# normalized_counts_all_with_rowname <- data.frame(hgnc_symbol = rownames(normalized_counts_all), normalized_counts_all)
# # Write this new data frame to a tab-separated file
# write.table(normalized_counts_all_with_rowname, file = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\RNA-seq_Data\\TCells\\normalized.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

HD_columns <- grepl("_HD_", colnames(TCell_gene_expression))
HIV_columns <- grepl("_HIV_", colnames(TCell_gene_expression))

df_HD <- TCell_gene_expression[, HD_columns]
df_HIV <- TCell_gene_expression[, HIV_columns]

# Extract column names that match each group or the control
Unstim_columns <- colnames(df_HD)[grepl("_Unstim_", colnames(df_HD))]
IFNa_columns <- colnames(df_HD)[grepl("_IFNa_", colnames(df_HD))]
IL27_columns <- colnames(df_HD)[grepl("_IL27_", colnames(df_HD))]
IL6_columns <- colnames(df_HD)[grepl("_IL6_", colnames(df_HD))]

# Create data frames for each group along with their respective control
df_IFNa <- df_HD[, c(Unstim_columns, IFNa_columns)]
df_IL27 <- df_HD[, c(Unstim_columns, IL27_columns)]
df_IL6 <- df_HD[, c(Unstim_columns, IL6_columns)]

# Create a list of data frames for each condition
list_of_dfs <- list(IFNa = df_IFNa, IL27 = df_IL27, IL6 = df_IL6)

# Loop through the list and apply DESeq2
for (group_name in names(list_of_dfs)) {
  
  # Subset DataFrame for this loop iteration
  sub_df <- list_of_dfs[[group_name]]
  
  # Create colData DataFrame for this loop, matching the columns in 'sub_df'
  conditions <- ifelse(grepl("_Unstim_", colnames(sub_df)), "control", "treatment")
  colData <- data.frame(row.names = colnames(sub_df), condition = conditions)
  
  # Replace all NA values with 0 in the data frame
  sub_df[is.na(sub_df)] <- 0
  
  
  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = sub_df,
                                colData = colData,
                                design = ~ condition)
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Extract results and filter by log2 fold change
  res <- results(dds)
  significant_genes <- subset(res, abs(log2FoldChange) > log2fc_cutoff)
  
  # Extract log2 fold changes and store them in a vector
  lfc_values <- significant_genes$log2FoldChange
  
  # Add these to the final data frame with prefix TCELL_
  if (is.null(final_df)) {
    final_df <- data.frame(row.names = rownames(significant_genes))
  }
  
  tcell_group <- paste0("TCELL_", group_name)
  final_df[[tcell_group]] <- NA
  final_df[rownames(significant_genes), tcell_group] <- lfc_values
}


rm(list = setdiff(ls(), c("final_df","log2fc_cutoff")))

final_df <- data.frame(hgnc_symbol = rownames(final_df), final_df)
# Write this new data frame to a tab-separated file
write.table(final_df, file = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\RNA-seq_Data\\foldChange_bigger50_interferons.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


final_df$hgnc_symbol = rownames(final_df)
final_df$hgnc_symbol = NULL


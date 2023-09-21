library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)

#####
#SNPS
#total number of SNPs in 3 spacers and 4 spacers
#allSNPs = read.csv(file = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\found_SNPs\\hg38_all_potential_SNPs.sorted.vcf", sep = "\t", comment.char = "#", header = FALSE)
allSNPs = read.csv(file = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\found_SNPs\\PBMC_TCell\\SNPs_in_highAcetylation.vcf", sep = "\t", comment.char = "#", header = FALSE)
allSNPs <- subset(allSNPs, select = c(V1,V8))


word_counts <- allSNPs %>% 
  group_by(V8) %>% 
  summarise(count = n())

new_max_value = max(word_counts$count) + 5000

ggplot(word_counts, aes(x = V8, y = count)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = function(x) formatC(x, format = "f", big.mark = ",", digits = 0), limits = c(0, new_max_value)) +
  xlab("Words") +
  ylab("Frequency") +
  ggtitle("Frequency of SNPs in 3-spacers and 4-spacers")

#number of SNPs per chromosome of 3spacers and 4spacers
word_counts_grouped <- allSNPs %>%
  group_by(V1, V8) %>%
  summarise(count = n())

word_counts_grouped <- word_counts_grouped %>% 
  mutate(is_numeric = ifelse(grepl("^\\d+$", V1), 1, 0),
         numeric_value = ifelse(is_numeric == 1, as.integer(V1), NA_integer_))
word_counts_grouped <- word_counts_grouped %>% 
  arrange(is_numeric, desc(is_numeric), numeric_value, V1) %>% 
  select(-is_numeric, -numeric_value)


ggplot(word_counts_grouped, aes(x = factor(V1, levels = unique(word_counts_grouped$V1)), y = count, fill = V8)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("Words") +
  ylab("Frequency") +
  ggtitle("Frequency of SNPs in 3-spacers and 4-spacers per chromosome")

word_counts_SNPs <- word_counts
word_counts_grouped_SNPs <- word_counts_grouped

##### 
#MOTIFS
#total number of GAS motifs in 3 spacers and 4 spacers

ALLmotifs = read.csv(file = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\ALL_motifs\\FIMO\\allMotifs.sorted.bed", sep = "\t", header = FALSE)

word_counts <- ALLmotifs %>% 
  group_by(V4) %>% 
  summarise(count = n())

ggplot(word_counts, aes(x = V4, y = count)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = function(x) formatC(x, format = "f", big.mark = ",", digits = 0)) +
  xlab("Words") +
  ylab("Frequency") +
  ggtitle("Frequency of GAS motifs in hg38")

#number of GAS motifs per chromosome of 3spacers and 4spacers
word_counts_grouped <- ALLmotifs %>%
  group_by(V1, V4) %>%
  summarise(count = n())

word_counts_grouped$V1 <- sub("^chr", "", word_counts_grouped$V1)


word_counts_grouped <- word_counts_grouped %>% 
  mutate(is_numeric = ifelse(grepl("^\\d+$", V1), 1, 0),
         numeric_value = ifelse(is_numeric == 1, as.integer(V1), NA_integer_))
word_counts_grouped <- word_counts_grouped %>% 
  arrange(is_numeric, desc(is_numeric), numeric_value, V1) %>% 
  select(-is_numeric, -numeric_value)


ggplot(word_counts_grouped, aes(x = factor(V1, levels = unique(word_counts_grouped$V1)), y = count, fill = V4)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("Words") +
  ylab("Frequency") +
  ggtitle("Frequency of GAS motifs in 3-spacers and 4-spacers per chromosome")

word_counts_motifs = word_counts
word_counts_grouped_motifs = word_counts_grouped

#####
#RATIO

colnames(word_counts_motifs)<-c("motif","count")
colnames(word_counts_SNPs)<-c("motif","count")

joined_df_allMotifs <-  word_counts_motifs %>%
  inner_join(word_counts_SNPs, by = c("motif")) %>%
  mutate(ratio = word_counts_SNPs$count / word_counts_motifs$count)

colnames(word_counts_grouped_motifs)<-c("chr","motif","count")
colnames(word_counts_grouped_SNPs)<-c("chr","motif","count")
word_counts_grouped_motifs$chr <- sub("^chr", "", word_counts_grouped_motifs$chr)
word_counts_grouped_motifs$chr <- as.character(word_counts_grouped_motifs$chr)
word_counts_grouped_SNPs$chr <- as.character(word_counts_grouped_SNPs$chr)

word_counts_grouped_motifs$count <- as.integer(word_counts_grouped_motifs$count)
word_counts_grouped_SNPs$count <- as.integer(word_counts_grouped_SNPs$count)

word_counts_grouped_SNPs <- ungroup(word_counts_grouped_SNPs)
word_counts_grouped_motifs <- ungroup(word_counts_grouped_motifs)


joined_df_chrMotifs <-  word_counts_grouped_motifs %>%
  inner_join(word_counts_grouped_SNPs, by = c("chr","motif")) %>%
  mutate(ratio = word_counts_grouped_SNPs$count / word_counts_grouped_motifs$count)

joined_df_chrMotifs <- joined_df_chrMotifs %>% 
  mutate(is_numeric = ifelse(grepl("^\\d+$", chr), 1, 0),
         numeric_value = ifelse(is_numeric == 1, as.integer(chr), NA_integer_))
joined_df_chrMotifs <- joined_df_chrMotifs %>% 
  arrange(is_numeric, desc(is_numeric), numeric_value, chr) %>% 
  select(-is_numeric, -numeric_value)

ggplot(joined_df_chrMotifs, aes(x = factor(chr, levels = unique(joined_df_chrMotifs$chr)), y = ratio, fill = motif)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(limits = c(0, 1.1)) +  # Extend y-axis to include 1.0
  xlab("Chromosome") +
  ylab("Ratio") +
  ggtitle("Grouped Barplot of Ratios by Chromosome")



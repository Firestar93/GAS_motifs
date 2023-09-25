import pandas as pd
import re

def extract_genes(text):
    # The regex pattern to match gene names enclosed between single quotes following the "gene=" pattern
    gene_pattern = r"gene='(.*?)'"
    genes = re.findall(gene_pattern, text)
    return ';'.join(genes)

def filter_genes(genes_str):
    genes_list = genes_str.split(';')
    filtered_genes = [gene for gene in genes_list if gene in hgnc_set]
    return ';'.join(filtered_genes)

output_file = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\found_SNPs\\PBMC_TCell\\SNPs_in_highAcetylation_AND_immuneGenes_AND_STATsignals_AND_IFNorILinducedFC50.vcf"

induced_genes = pd.read_csv('C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\RNA-seq_Data\\foldChange_bigger50_interferons.tsv', sep="\t", comment="#")
hgnc_set = set(induced_genes['hgnc_symbol'])


all_SNPs = pd.read_csv('C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\found_SNPs\\PBMC_TCell\\SNPs_in_highAcetylation_AND_immuneGenes_AND_STATsignals.vcf', sep="\t", comment="#", header=None)
all_SNPs = all_SNPs.rename(columns={0 : '#CHROM',1: 'POS',2  : 'ID', 3: 'REF',4 : 'ALT', 5 : 'QUAL',6 : 'FILTER', 7 : 'INFO', 8 : 'FORMAT',9 : 'NA19909'})

all_SNPs['genes'] = all_SNPs['INFO'].apply(extract_genes)

all_SNPs['filtered_genes'] = all_SNPs['genes'].apply(filter_genes)

all_SNPs = all_SNPs.loc[all_SNPs['filtered_genes'] != '']
all_SNPs = all_SNPs.drop('filtered_genes', axis=1)


with open(output_file, 'w') as f:
    f.write("##fileformat=VCFv4.3\n##fileDate=20090805\n##source=myImputationProgramV3.1\n##reference=file:///seq/references/1000GenomesPilot-NCBI36.fastas\n")

all_SNPs.to_csv(output_file, index=False, sep='\t', index_label=False, mode= "a")



print("X")




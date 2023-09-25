import pandas as pd

output_file = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\found_SNPs\\hg38_casein_region_SNPs"
all_SNPs = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\found_SNPs\\hg38_all_potential_SNPs.sorted.vcf"
df_snps = pd.read_csv(all_SNPs, sep='\t', comment="#", header=None)


chrom = "4"
start = 69871707
end = 70311230

df_snps[0] = df_snps[0].astype(str)
df_snps[1] = df_snps[1].astype(int)

filtered_df = df_snps[(df_snps[1] >= start) & (df_snps[1] <= end) & (df_snps[0] == chrom)]

filtered_df = filtered_df.rename(columns={0 : '#CHROM',1: 'POS',2  : 'ID', 3: 'REF',4 : 'ALT', 5 : 'QUAL',6 : 'FILTER', 7 : 'INFO', 8 : 'FORMAT',9 : 'NA19909'})

with open(output_file, 'w') as f:
    f.write("##fileformat=VCFv4.3\n##fileDate=20090805\n##source=myImputationProgramV3.1\n##reference=file:///seq/references/1000GenomesPilot-NCBI36.fastas\n")

filtered_df.to_csv(output_file, index=False, sep='\t', index_label=False, mode= "a")

print('X')

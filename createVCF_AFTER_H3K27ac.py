import pandas as pd

df_snps = pd.read_csv("/home/markus/Dropbox/UNI/NIH_Projects/GAS_motifs_SNPs/found_SNPs/SurvivingSNPs_AFTERH3K27ac.bed", index_col=False, sep='\t')

df_snps = df_snps.drop_duplicates(subset=['ID'])
df_snps = df_snps.drop(columns=["N", "N.1", "N.2", "N.3", "N CHRO", "START","END","ELEID","QUAL"])
df_snps['QUAL'] = '.'
df_snps['FILTER'] = '.'
df_snps['INFO'] = '.'
df_snps['FORMAT'] = '.'
df_snps['NA19909'] = '.'

df_snps.to_csv("/home/markus/Dropbox/UNI/NIH_Projects/GAS_motifs_SNPs/found_SNPs/SurvivingSNPs_AFTERH3K27ac.vcf", index=False, sep='\t', index_label=False)


print("X")


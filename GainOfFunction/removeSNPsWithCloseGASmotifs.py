import pandas as pd

all_SNPs_input = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\GainOfFunction\\SNPs_in_highAcetylation_AND_immuneGenes_AND_correctSignTarget_CLEANED_upstreamLess10KB.vcf"
too_close_SNPs = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\GainOfFunction\\upstream10KB_has_close_GAS_motif.vcf"

output_cleaned = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\GainOfFunction\\SNPs_in_highAcetylation_AND_immuneGenes_AND_correctSignTarget_CLEANED_upstreamLess10KB_noGAS_200bp.vcf"

df_allSNPs = pd.read_csv(all_SNPs_input, sep='\t', comment="#", header=None)
df_SNPs_with_tooCloseGAS = pd.read_csv(too_close_SNPs, sep='\t', comment="#", header=None)

df_allSNPs.columns = ['CHR','POS','RSID','REF','ALT','QUAL','FILTER','INFO','FORMAT']
df_SNPs_with_tooCloseGAS.columns = ['CHR','POS','RSID','REF','ALT','QUAL','FILTER','INFO','FORMAT']

# Strip whitespace from the RSID columns just in case there is any.
df_allSNPs['RSID'] = df_allSNPs['RSID'].str.strip()
df_SNPs_with_tooCloseGAS['RSID'] = df_SNPs_with_tooCloseGAS['RSID'].str.strip()

# Ensure that both RSID columns are strings.
df_allSNPs['RSID'] = df_allSNPs['RSID'].astype(str)
df_SNPs_with_tooCloseGAS['RSID'] = df_SNPs_with_tooCloseGAS['RSID'].astype(str)

df_unique_SNPs = df_allSNPs[~df_allSNPs['RSID'].isin(df_SNPs_with_tooCloseGAS['RSID'])]

with open(output_cleaned, 'w') as f:
    f.write("##fileformat=VCFv4.3\n##fileDate=20090805\n##source=myImputationProgramV3.1\n##reference=file:///seq/references/1000GenomesPilot-NCBI36.fastas\n#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA19909\n")
# Write vcf_data_filtered to a file without the index and header
df_unique_SNPs.to_csv(output_cleaned, sep='\t', index=False, header=False, mode='a')


print("X")



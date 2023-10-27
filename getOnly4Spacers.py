import pandas as pd

all_SNPs = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\found_SNPs\\PBMC_TCell\\SNPs_in_highAcetylation_AND_immuneGenes_AND_STATsignals.vcf"
output_file = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\found_SNPs\\PBMC_TCell\\SNPs_in_highAcetylation_AND_immuneGenes_AND_STATsignals_4SPACERSONLY.vcf"

df_snps = pd.read_csv(all_SNPs, sep='\t', comment="#", header=None)

df_snps_4spacers = df_snps[df_snps[7].str.startswith('4spacer')]

with open(output_file, 'w') as f:
    f.write("##fileformat=VCFv4.3\n##fileDate=20090805\n##source=myImputationProgramV3.1\n##reference=file:///seq/references/1000GenomesPilot-NCBI36.fastas\n")

df_snps_4spacers.to_csv(output_file, index=False, sep='\t', index_label=False, mode= "a")


print("X")
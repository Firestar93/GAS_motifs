import pandas as pd

output_cleaned_VCF = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\GAS_motifs\\GainOfFunction\\SNPs_in_highAcetylation_AND_immuneGenes_AND_correctSignTarget_CLEANED.vcf"

input_filtered_VCF = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\GAS_motifs\\GainOfFunction\\SNPs_in_highAcetylation_AND_immuneGenes_AND_correctSignTarget.vcf"
input_SNPs = pd.read_csv(input_filtered_VCF, comment='#', sep='\t', header=0)
input_SNPs.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']

input_SNPs = input_SNPs.drop(columns='POS')

input_correct_SNP_positions = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\GAS_motifs\\GainOfFunction\\SNPs_in_highAcetylation_AND_immuneGenes.vcf"
input_SNPs_correct_positions = pd.read_csv(input_correct_SNP_positions, comment='#', sep='\t', header=0)
input_SNPs_correct_positions.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']

merged_df = input_SNPs.merge(input_SNPs_correct_positions[['ID', 'POS']],
                                         on='ID',
                                         how='left')

column_order = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
merged_df = merged_df.reindex(columns=column_order)

with open(output_cleaned_VCF, 'w') as f:
    f.write("##fileformat=VCFv4.3\n##fileDate=20090805\n##source=myImputationProgramV3.1\n##reference=file:///seq/references/1000GenomesPilot-NCBI36.fastas\n#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA19909\n")
# Write vcf_data_filtered to a file without the index and header
merged_df.to_csv(output_cleaned_VCF, sep='\t', index=False, header=False, mode='a')


print("X")


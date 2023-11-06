import pandas as pd

vcf_file = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\GAS_motifs\\GainOfFunction\\SNPs_in_highAcetylation_AND_immuneGenes.vcf"
output_filtered_VCF = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\GAS_motifs\\GainOfFunction\\SNPs_in_highAcetylation_AND_immuneGenes_AND_correctSignTarget.vcf"

vcf_data = pd.read_csv(vcf_file, comment='#', sep='\t', header=0)
vcf_data.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']

output_filtered_VCF_data = pd.read_csv(output_filtered_VCF, comment='#', sep='\t', header=0)
output_filtered_VCF_data.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']

unique_in_output_filtered = output_filtered_VCF_data[~output_filtered_VCF_data['ID'].isin(vcf_data['ID'])]
unique_in_vcf_data = vcf_data[~vcf_data['ID'].isin(output_filtered_VCF_data['ID'])]
unique_rows = pd.concat([unique_in_output_filtered, unique_in_vcf_data])


print("X")

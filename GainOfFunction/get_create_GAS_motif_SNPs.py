import csv
import pandas as pd
import os
from datetime import datetime

GAS_file = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\GAS_motifs\ALL_motifs\\FIMO_gainOfFunction\\all_motifs.sorted.bed"

#GAS_dir = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\ALL_motifs\\FIMO_gainOfFunction"
dbSNP_dir = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\GAS_motifs\\01_INITIAL_WORK_Shreeti\\SNP\\SNP files"
#file_output_SNPs = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\GainOfFunction_SNPs\\potential_gainOfFunction_SNPs.tsv"
file_output_GAS = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\GainOfFunction\\all_potential_GAS_motifs.bed"
file_output_SNPs_AND_GAS = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\GainOfFunction\\potential_GAS_with_SNPs.csv"
file_output_only_SNPs = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\GainOfFunction\\SNPs_creating_GAS.vcf"

with open(file_output_only_SNPs, 'w') as f:
    f.write("##fileformat=VCFv4.3\n##fileDate=20090805\n##source=myImputationProgramV3.1\n##reference=file:///seq/references/1000GenomesPilot-NCBI36.fastas\n#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA19909\n")

df_GAS = pd.read_csv(GAS_file, sep='\t', comment="#", header=None)
df_GAS.columns = ['chr', 'start', 'end', 'motif_type']

output_data = []

for dirname in sorted(os.listdir(dbSNP_dir)):
    dir1 = os.path.join(dbSNP_dir, dirname)
    if os.path.isfile(dir1):
        print(dir1)

        if dirname.startswith("chr"):
            num = dirname[3:].split(".")[0]

            current_time = datetime.now()
            print("[START] Read in dbSNP @ " + str(current_time))

            snp_results = pd.read_csv(dir1, sep='\t', comment="#")

            current_time = datetime.now()
            print("[START] Check out SNPs @ " + str(current_time))

            df_GAS['chr'] = df_GAS['chr'].astype(str)
            df_chr_GAS_motifs =df_GAS[df_GAS['chr'] == num]

            if df_chr_GAS_motifs.empty:
                continue

            calculate_i = lambda row: (row['start']+1) if row['motif_type'].startswith("T1_GAS") else \
                (row['start'] + 2) if row['motif_type'].startswith("T2_GAS") else \
                    (row['start'] + 3) if row['motif_type'].startswith("C_GAS") else \
                        (row['end'] - 1) if row['motif_type'].startswith("G_GAS") else \
                            (row['end'] ) if row['motif_type'].startswith("A1_GAS") else \
                                (row['end'] + 1) if row['motif_type'].startswith("A2_GAS") else \
                                    None  # or some default value if none of the conditions match
            df_chr_GAS_motifs['SNP_position'] = df_chr_GAS_motifs.apply(calculate_i, axis=1)

            df_chr_GAS_motifs = df_chr_GAS_motifs.merge(snp_results[['RSID', 'CHR', 'REF/ALT', 'STOP']],
                                                        left_on='SNP_position',
                                                        right_on='STOP',
                                                        how='left')
            df_chr_GAS_motifs['start'] = df_chr_GAS_motifs['start'] + 1
            df_chr_GAS_motifs['end'] = df_chr_GAS_motifs['end'] + 1

            df_chr_GAS_motifs[['chr', 'start', 'end', 'motif_type']].to_csv(file_output_GAS, mode='a', header=False,
                                                                            index=False, sep='\t')

            df_chr_GAS_motifs = df_chr_GAS_motifs.dropna()

            df_chr_GAS_motifs.to_csv(file_output_SNPs_AND_GAS, mode='a', header=False,
                                    index=False, sep='\t')

            df_VCF = df_chr_GAS_motifs[['CHR', 'STOP', 'RSID', 'REF/ALT', 'motif_type']].rename(
                columns={'motif_type': 'INFO', 'STOP': 'POS'}).copy()

            split_columns = df_VCF['REF/ALT'].str.split('/', expand=True)
            df_VCF['REF'] = split_columns[0]
            df_VCF['ALT'] = split_columns[1]
            df_VCF = df_VCF.drop(columns=['REF/ALT'])

            condition = (df_VCF['REF'] == 'lengthTooLong') | (df_VCF['ALT'] == 'lengthTooLong')
            df_VCF.loc[condition, 'REF'] = 'T'
            df_VCF.loc[condition, 'ALT'] = 'A'

            df_VCF['QUAL'] = '.'
            df_VCF['FILTER'] = '.'
            df_VCF['FORMAT'] = '.'
            df_VCF['NA19909'] = '.'

            desired_order = ['CHR', 'POS', 'RSID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'NA19909']
            df_VCF = df_VCF[desired_order]
            df_VCF['POS'] = df_VCF['POS'].astype(int)

            df_VCF.to_csv(file_output_only_SNPs, mode='a', header=False, index=False, sep='\t')

            print("X")


print("X")



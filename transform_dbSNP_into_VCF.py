import os
import pandas as pd
import chardet

input_dir = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\GAS_motifs\\01_INITIAL_WORK_Shreeti\\SNP\\SNP files"


# iterate over files in
# that directory

for filename in sorted(os.listdir(input_dir)):
    f = os.path.join(input_dir, filename)
    # checking if it is a file
    if os.path.isfile(f):
        print(f)

        with open(f, "rb") as f_X:
            result = chardet.detect(f_X.read())

        df_snps = pd.read_csv(f, sep='\t', comment="#", encoding=result['encoding'])
        df_snps = df_snps.drop(columns=["START", "LOCATION"])
        df_snps = df_snps.rename(columns={'REF/ALT': 'REF_ALT'})
        df_snps = df_snps.join(pd.DataFrame([i.split('/', 1) for i in df_snps.REF_ALT], columns=['REF', 'ALT']))
        df_snps['ALT'] = 'N'
        df_snps = df_snps.drop(columns=["REF_ALT"])

        df_snps = df_snps.rename(columns={'STOP': 'POS'})
        df_snps = df_snps.rename(columns={'RSID': 'ID'})
        df_snps = df_snps.rename(columns={'CHR': '#CHROM'})
        df_snps['QUAL'] = '.'
        df_snps['FILTER'] = '.'
        df_snps['INFO'] = '.'

        df_snps['#CHROM'] = df_snps['#CHROM'].str.replace('chr', '', regex=True)

        df_snps.to_csv(
            "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\GAS_motifs\\01_INITIAL_WORK_Shreeti\\SNP\\SNP files\\allSNPs.vcf", index=False, sep='\t', index_label=False, mode='a', header=False)

        print("x")

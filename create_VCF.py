import os
import pandas as pd

input_files = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\GAS_motifs\\ALL_motifs\\FIMO"
output_file = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\GAS_motifs\\found_SNPs\\hg38_all_potential_SNPs.vcf"

df_snps_ALL = pd.DataFrame(columns=['#CHROM','POS', 'ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','NA19909'])

# iterate over files in
# that directory
for filename in sorted(os.listdir(input_files)):
    f = os.path.join(input_files, filename)
    # checking if it is a file
    if os.path.isdir(f):
        print(f)

        f_snps = os.path.join(f, "SNP_that_destroy_motifs.tsv")

        f_motif_GFF = os.path.join(f, "results", "fimo_adjusted.gff")

        df_snps = pd.read_csv(f_snps, sep='\t')
        df_snps = df_snps.drop(columns=["START","LOCATION"])

        df_snps=df_snps.iloc[:,[1,3,0,2]]

        df_snps = df_snps.rename(columns={'REF/ALT': 'REF_ALT'})
        df_snps = df_snps.join(pd.DataFrame([i.split('/', 1) for i in df_snps.REF_ALT], columns=['REF', 'ALT']))
        df_snps = df_snps.drop(columns=["REF_ALT"])
        df_snps['ALT'] = df_snps['ALT'].str.split('/').str.get(0)
        #df_snps['ALT'] = df_snps['ALT'].str.strip()
        #df_snps['REF'] = df_snps['REF'].str.strip()
        #df_snps['ALT'] = df_snps['ALT'].str.replace('lengthTooLong', 'A')
        #df_snps['REF'] = df_snps['REF'].str.replace('lengthTooLong', 'A')
        #df_snps['ALT'] = df_snps['ALT'].str.replace(' ', 'T')
        #df_snps['REF'] = df_snps['REF'].str.replace(' ', 'T')
        df_snps = df_snps.drop(df_snps[df_snps.isin(["lengthTooLong"]).any(axis=1)].index)
        df_snps.reset_index(drop=True, inplace=True)

        df_snps['QUAL'] = '.'
        df_snps['FILTER'] = '.'
        df_snps['INFO'] = '.'
        df_snps['FORMAT'] = '.'
        df_snps['NA19909'] = '.'

        #get information if 3spacer or 4spacer
        gff = pd.read_csv(f_motif_GFF, comment="#", sep="\t", header=None)
        for i in range(len(df_snps)):
            position = df_snps.loc[i, 'STOP']
            is_inside = (gff[3] <= position) & (gff[4] >= position)
            row_number = is_inside.idxmax()
            sequence = gff.loc[row_number, 13]
            if len(sequence)==19:
                df_snps.loc[i, 'INFO'] = "4spacer"
            if len(sequence) == 18:
                df_snps.loc[i, 'INFO'] = "3spacer"


        df_snps = df_snps.rename(columns={'STOP': 'POS'})
        df_snps = df_snps.rename(columns={'RSID': 'ID'})
        df_snps = df_snps.rename(columns={'CHR': '#CHROM'})

        data = [df_snps_ALL, df_snps]
        df_snps_ALL = pd.concat(data)

df_snps_ALL['#CHROM'] = df_snps_ALL['#CHROM'].str.replace('chr', '', regex=True)

with open(output_file, 'w') as f:
    f.write("##fileformat=VCFv4.3\n##fileDate=20090805\n##source=myImputationProgramV3.1\n##reference=file:///seq/references/1000GenomesPilot-NCBI36.fastas\n")

df_snps_ALL.to_csv(output_file, index=False, sep='\t', index_label=False, mode= "a")



print("X")


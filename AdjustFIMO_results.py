import os
import re
import csv
import pandas as pd


if __name__ == "__main__":

    input_dir = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\ALL_motifs\\FIMO"


    for filename in sorted(os.listdir(input_dir)):
        folder_chr = os.path.join(input_dir, filename)
        # checking if it is a file
        if os.path.isdir(folder_chr):
            print(filename + " START")

            folder_chr_results = os.path.join(folder_chr, "results")

            file_chr_results = os.path.join(folder_chr_results, "fimo.gff")
            file_fimo_adjusted = os.path.join(folder_chr_results, "fimo_adjusted.gff")


            all_GAS_motifs = pd.read_csv(file_chr_results, sep='\t', comment="#", header=None)

            split_columns = all_GAS_motifs[8].str.split(';', expand=True)
            all_GAS_motifs = pd.concat([all_GAS_motifs, split_columns], axis=1)
            all_GAS_motifs.columns = range(all_GAS_motifs.shape[1])

            #remove all that have mistakes
            pattern = r'TTC.{3,4}GAA'
            column_name = all_GAS_motifs.columns[13]
            all_GAS_motifs = all_GAS_motifs[all_GAS_motifs[column_name].str.contains(pattern, regex=True, na=False)]

            column_name = all_GAS_motifs.columns[3]
            all_GAS_motifs[column_name] = all_GAS_motifs[column_name] + 1

            column_name = all_GAS_motifs.columns[4]
            all_GAS_motifs[column_name] = all_GAS_motifs[column_name] + 1

            with open(file_fimo_adjusted, 'w') as f:
                f.write("##gff-version 3" + '\n')

            all_GAS_motifs.to_csv(file_fimo_adjusted, mode="a", sep="\t", header=False, index=False)


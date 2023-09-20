import os
import re
import csv
import pandas as pd
from datetime import datetime


if __name__ == "__main__":

    input_dir = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\ALL_motifs\\FIMO"
    output_file = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\ALL_motifs\\FIMO\\allMotifs.bed"

    df_MOTIFS_ALL = pd.DataFrame(columns=["CHR","START","END","NAME"])

    for filename in sorted(os.listdir(input_dir)):
        folder_chr = os.path.join(input_dir, filename)
        # checking if it is a file
        if os.path.isdir(folder_chr):
            current_time = datetime.now()
            print(filename + " START @ " + str(current_time))

            folder_chr_results = os.path.join(folder_chr, "results")
            gff_file = os.path.join(folder_chr_results, "fimo_adjusted.gff")

            all_GAS_motifs = pd.read_csv(gff_file, sep='\t', comment="#", header=None)
            split_columns = all_GAS_motifs[8].str.split(';', expand=True)
            all_GAS_motifs = pd.concat([all_GAS_motifs, split_columns], axis=1)
            all_GAS_motifs.columns = range(all_GAS_motifs.shape[1])

            for i in range(len(all_GAS_motifs)):

                if len(all_GAS_motifs.loc[i,13]) == 18:
                    all_GAS_motifs.loc[i, 13] = "3spacer"
                if len(all_GAS_motifs.loc[i, 13]) == 19:
                    all_GAS_motifs.loc[i, 13] = "4spacer"

            all_GAS_motifs = all_GAS_motifs.reindex(columns=all_GAS_motifs.columns)

            indices_to_keep = [0,3,4,13]
            all_GAS_motifs = all_GAS_motifs.iloc[:, indices_to_keep]

            new_column_names = {
                0: 'CHR',
                3: 'START',
                4: 'END',
                13: 'NAME'
            }

            all_GAS_motifs.rename(columns=new_column_names, inplace=True)

            data = [df_MOTIFS_ALL, all_GAS_motifs]
            df_MOTIFS_ALL = pd.concat(data)

            current_time = datetime.now()
            print(filename + " END @ " + str(current_time))

    #write out GFF
    df_MOTIFS_ALL.to_csv(output_file, index=False, header=False, sep="\t")




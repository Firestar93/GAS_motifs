import os
import pandas as pd

input_dir = "/home/markus/Dropbox/UNI/NIH_Projects/GAS_motifs_SNPs/found_SNPs/bed_files"

df_beds_ALL = pd.DataFrame()

# iterate over files in
# that directory
for filename in sorted(os.listdir(input_dir)):
    f = os.path.join(input_dir, filename)
    # checking if it is a file
    if os.path.isfile(f):
        print(f)

        df_bed = pd.read_csv(f, sep='\t', header = None)
        df_bed[3] = 'NAME'

        data = [df_beds_ALL, df_bed]
        df_beds_ALL = pd.concat(data)


df_beds_ALL.to_csv(input_dir+"/combined_bed.bed", index=False, sep='\t', index_label=False)


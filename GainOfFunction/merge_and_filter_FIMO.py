import pandas as pd
from itertools import product
import os

def expand_motif(motif, substitutions):
    positions = []
    options = []
    for i, letter in enumerate(motif):
        if letter in substitutions:
            positions.append(i)
            options.append(substitutions[letter])
        else:
            options.append(letter)

    for p in product(*options):
        yield ''.join(p)


def generate_replacement_dict(motifs, motif_to_id, substitutions):
    replacement_dict = {}

    for motif, motif_id in motif_to_id.items():
        options = []

        for letter in motif:
            options.append(substitutions.get(letter, letter))

        for item in product(*options):
            expanded_motif = ''.join(item)
            replacement_dict[expanded_motif] = motif_id

    return replacement_dict

GAS_dir = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\GAS_motifs\\ALL_motifs\\FIMO_gainOfFunction"
GAS_out_file = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\GAS_motifs\\ALL_motifs\\FIMO_gainOfFunction\\all_motifs.bed"

for dirname in sorted(os.listdir(GAS_dir)):
    dir1 = os.path.join(GAS_dir, dirname)
    if dirname!="chr1":
        continue
    if os.path.isdir(dir1):
        print(dir1)
        for dirname2 in sorted(os.listdir(dir1)):
            dir2 = os.path.join(dir1, dirname2)
            if os.path.isdir(dir2):
                print(dir2)
                GAS_file = os.path.join(dir2, "fimo.tsv")
                print(GAS_file)
                df_GAS = pd.read_csv(GAS_file, sep='\t', comment="#")
                num_rows = len(df_GAS)
                print("Number of rows:", str(num_rows))

                # Define the nucleotide substitutions
                substitutions = {
                    'V': ['A', 'C', 'G'],
                    'D': ['A', 'T', 'G'],
                    'H': ['A', 'T', 'C'],
                    'B': ['G', 'T', 'C'],
                    'N': ['A', 'G', 'T', 'C']
                }

                # Define the motifs
                motifs = [
                    "VTCNNNGAA",
                    "TVCNNNGAA",
                    "TTDNNNGAA",
                    "TTCNNNHAA",
                    "TTCNNNGBA",
                    "TTCNNNGAB"
                ]

                # Create the hash set to store unique expanded motifs
                motif_set = set()

                # Expand and add each motif to the set
                for motif in motifs:
                    for expanded_motif in expand_motif(motif, substitutions):
                        motif_set.add(expanded_motif)

                df_GAS = df_GAS[df_GAS['matched_sequence'].isin(motif_set)]
                num_rows = len(df_GAS)
                print("Number of rows:", str(num_rows))

                motif_to_id = {
                    'VTCNNNGAA': 'T1_GAS',
                    'TVCNNNGAA': 'T2_GAS',
                    'TTDNNNGAA': 'C_GAS',
                    'TTCNNNHAA': 'G_GAS',
                    'TTCNNNGBA': 'A1_GAS',
                    'TTCNNNGAB': 'A2_GAS'
                }

                replacement_dict = generate_replacement_dict(motif_to_id.keys(), motif_to_id, substitutions)

                # make it a bed file format
                df_GAS['motif_alt_id'] = df_GAS['matched_sequence'].replace(replacement_dict)
                df_GAS.drop(['motif_id', 'strand', 'score', 'q-value', 'p-value'], axis=1, inplace=True)
                df_GAS['info'] = df_GAS['motif_alt_id'].astype(str) + '-' + df_GAS['matched_sequence'].astype(str)
                df_GAS.drop(['motif_alt_id', 'matched_sequence'], axis=1, inplace=True)
                column_order = [col for col in df_GAS.columns if col != 'info'] + ['info']
                df_GAS = df_GAS[column_order]
                df_GAS['sequence_name'] = df_GAS['sequence_name'].str.replace('chr', '', regex=False)
                df_GAS = df_GAS.drop_duplicates()

                df_GAS.to_csv(GAS_out_file, mode='a', header=False, index=False, sep="\t")

print("X")


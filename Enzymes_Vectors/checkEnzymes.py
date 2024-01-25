import csv
import os
import re

def read_and_process_fasta(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()[1:]  # Skip the first line
        return ''.join(line.strip() for line in lines)

def read_csv_to_dict(file_path):
    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        return {rows[0]: convert_to_regex(rows[1]) for rows in reader}

def convert_to_regex(sequence):
    return re.sub('n', '[ATCG]', sequence)

def find_enzymes_with_no_hits(enzyme_sequence, sequences):
    enzymes_with_no_hits = []

    for enzyme, regex_pattern in enzyme_sequence.items():
        has_hit = False
        for sequence in sequences.values():
            if re.search(regex_pattern, sequence):
                has_hit = True
                break  # Break out of the loop if a match is found

        if not has_hit:
            enzymes_with_no_hits.append(enzyme)

    return enzymes_with_no_hits

def save_sequence_to_fasta(sequence, key, directory):
    filename = f"{key}_sequence_including_enzymeSite.fa"
    file_path = os.path.join(directory, filename)

    with open(file_path, 'w') as file:
        file.write(f"> {key}\n")  # FASTA header
        for i in range(0, len(sequence), 50):
            file.write(sequence[i:i+50] + "\n")
directory_pre ="C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\GainOfFunction\\presentations\\Luciferase"
directory = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\GainOfFunction\\presentations\\Luciferase\\Sequence_changedGAS"  # Replace with your directory path
directory_out = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\GainOfFunction\\presentations\\Luciferase\\including_restriction_enzyme"

file_names = ["IL7R_DNA_sequence.fa", "IRF3_DNA_sequence.fa", "JAK2_DNA_sequence.fa", "JAK3_DNA_sequence.fa", "PTPN2_DNA_sequence.fa","SOCS1_DNA_sequence.fa"]  # Replace with your file names
csv_filename = "restriction_enzyme_sequences_and_productURL.csv"  # Replace with your CSV filename

csv_file_path = os.path.join(directory_pre, csv_filename)
enzyme_sequence = read_csv_to_dict(csv_file_path)

sequences = {}

for file_name in file_names:
    key = file_name.split('_')[0]  # Extract key from filename
    file_path = os.path.join(directory, file_name)
    sequences[key] = read_and_process_fasta(file_path)


enzymes_with_no_hits = find_enzymes_with_no_hits(enzyme_sequence, sequences)

print("Enzymes with no hits:")
for enzyme in enzymes_with_no_hits:
    print(enzyme)


prefix = "AATTGCTAGC"
suffix = "GCTAGCGGAT"

for key in sequences:
    sequences[key] = prefix + sequences[key] + suffix


for key, sequence in sequences.items():
    save_sequence_to_fasta(sequence, key, directory_out)

print("X")
import os
import re
import pandas as pd

def count_N_around_chromosome(sequence):
    pattern = r"(N+)?([ATGC]+)(N+)?"
    match = re.search(pattern, sequence)
    if match:
        before_chromosome = match.group(1)
        after_chromosome = match.group(3)
        len_before = len(before_chromosome) if before_chromosome else 0
        len_after = len(after_chromosome) if after_chromosome else 0
        return len_before, len_after
    else:
        return 0, 0


def sequence_to_fasta(identifier, sequence, line_length=80):
    # Initialize an empty list to store the FASTA lines
    fasta_lines = []

    # Add the identifier line
    fasta_lines.append(identifier)

    # Add sequence lines
    for i in range(0, len(sequence), line_length):
        subseq = sequence[i:i + line_length]
        fasta_lines.append(subseq)

    # Join the lines to form the FASTA string
    fasta_str = "\n".join(fasta_lines)
    return fasta_str

#input_file = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\ALL_motifs\\FIMO\\chr1.fa"
#output_file = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\ALL_motifs\\FIMO\\3gapper\\chr1_modified.fa"

input_dir = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\01_INITIAL_WORK_Shreeti\\Genome\\Genome sequence"
output_dir = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\ALL_motifs\\FIMO"

for filename in sorted(os.listdir(input_dir)):
    f = os.path.join(input_dir, filename)
    # checking if it is a file
    if os.path.isfile(f):
        chr_name = filename.split(".fa")[0]
        header = ">"+chr_name+":"

        folder_output = output_dir + "\\" + chr_name
        file_output = folder_output + "\\" + chr_name + ".fa"

        if not os.path.exists(folder_output):
            os.makedirs(folder_output)

        with open(f, 'r') as sequence_file:
            seq = sequence_file.read().replace('\n', '')

        seq = seq.replace(">" + chr_name, "")

        len_before, len_after = count_N_around_chromosome(seq)

        stop = len(seq) - len_after

        seq = seq[len_before:stop]

        header = header + str(len_before) + "-" + str(stop)

        fasta_format = sequence_to_fasta(header, seq)

        with open(file_output, 'w') as file:
            file.write(fasta_format)

        print("X")

print("X")

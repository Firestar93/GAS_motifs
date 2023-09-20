import re
import os
from pathlib import Path


# Function to find motif alignments in a sequence and write to BED file
def process_genome_sequence(sequence_file, chromosome_name, output_file):
    with open(sequence_file, 'r') as sequence_file:
        seq = sequence_file.read().replace('\n','')

    #motif = "TTC...GAA|TTC....GAA|AAG...CTT|AAG....CTT"

    motif = "TTC....GAA"


    alignments = re.findall(motif, seq, re.IGNORECASE)

    # Create a list to store position information and nucleotides
    bed_entries = []

    for alignment in alignments:
        start_pos = seq.index(alignment)
        end_pos = start_pos + len(alignment)

        bed_entry = f"{chromosome_name}\t{start_pos}\t{end_pos}"
        bed_entries.append(bed_entry)

    with open(output_file, 'a') as output_bed_file:
        for entry in bed_entries:
            output_bed_file.write(entry + '\n')

    print(chromosome_name)


def main_files(input_dir_genome, output_file):
    genome_dir = Path(input_dir_genome)
    #output_dir = Path(output_dir)

    for genome_file in genome_dir.glob('*.fa'):
        chr_name = genome_file.stem
        #output_file = output_dir / (chr_name + '_output_motif_alignment.bed')

        process_genome_sequence(genome_file, chr_name, output_file)

if __name__ == "__main__":
    #input_dir_genome = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\01_INITIAL_WORK_Shreeti\\Genome\\Genome sequence"
    input_dir_genome = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\enhancer_sequence"
    #output_dir = "C:\\Users\\hoffmannmd\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\ALL_motifs\\All_Motifs\\all_motifs.bed"
    output_dir = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\enhancer_sequence\\test.bed"


    main_files(input_dir_genome, output_dir)
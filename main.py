import csv
import re
# import os
from pathlib import Path


def snp_dict(file_path):
    information_dict = {}

    with open(file_path, 'r') as SNP_file:
        tsv_reader = csv.DictReader(SNP_file, delimiter='\t')

        for snp in tsv_reader:
            position = int(snp['STOP'])
            information_dict[position] = {
                'CHR': snp['CHR'],
                'RSID': snp['RSID'],
                'REF/ALT': snp['REF/ALT'],
                'START': snp['START'],
                'STOP': snp['STOP'],
                'LOCATION': snp['LOCATION']
            }

    return information_dict


def process_genome_sequence(genome_file, snp_results, output_file):
    with open(genome_file, 'r') as file:
        seq = file.read()
        seq = seq.replace('\n', '')

    motif = "TTC...GAA|TTC....GAA"
    alignments = re.findall(motif, seq)

    output_data = []

    for alignment in alignments:
        start_pos = seq.index(alignment)
        end_pos = start_pos + len(alignment)

        pos_to_check = list(range(start_pos + 1, start_pos + 4)) + list(range(end_pos - 2, end_pos + 1))

        for i in pos_to_check:
            actual_position = i
            if actual_position in snp_results:
                rs_id = snp_results[actual_position]['RSID']
                if not any(d['RSID'] == rs_id for d in output_data):
                    output_data.append({
                        'RSID': rs_id,
                        'CHR': snp_results[actual_position]['CHR'],
                        'REF/ALT': snp_results[actual_position]['REF/ALT'],
                        'START': snp_results[actual_position]['START'],
                        'STOP': snp_results[actual_position]['STOP'],
                        'LOCATION': snp_results[actual_position]['LOCATION']
                    })

    with open(output_file, 'w', newline='') as tsv_output_file:
        fieldnames = ['RSID', 'CHR', 'REF/ALT', 'START', 'STOP', 'LOCATION']
        tsv_writer = csv.DictWriter(tsv_output_file, fieldnames=fieldnames, delimiter='\t')
        tsv_writer.writeheader()
        tsv_writer.writerows(output_data)


def main_files(input_dir_dbsnp, input_dir_genome, output_dir):
    dbsnp_dir = Path(input_dir_dbsnp)
    genome_dir = Path(input_dir_genome)
    output_dir = Path(output_dir)

    for dbsnp_file in dbsnp_dir.glob('*.tsv'):
        chr_name = dbsnp_file.stem
        genome_file = genome_dir / (chr_name + '.fa')
        output_file = output_dir / (chr_name + '_snp_destroy_motif.tsv')

        if genome_file.exists():
            snp_results = snp_dict(dbsnp_file)
            process_genome_sequence(genome_file, snp_results, output_file)

if __name__ == "__main__":
    input_dir_dbsnp = "/Users/chhatralasd/Downloads/SNP"
    input_dir_genome = "/Users/chhatralasd/Downloads/Genome"
    output_dir = "/Users/chhatralasd/Downloads/found_snp_motif"

    main_files(input_dir_dbsnp, input_dir_genome, output_dir)
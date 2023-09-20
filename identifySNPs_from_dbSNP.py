import os
import re
import csv
import pandas as pd
from datetime import datetime


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

if __name__ == "__main__":

    input_dir = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\ALL_motifs\\FIMO"

    for filename in sorted(os.listdir(input_dir)):
        folder_chr = os.path.join(input_dir, filename)

        # checking if it is a file
        if os.path.isdir(folder_chr):

            current_time = datetime.now()
            print(filename + " START @ " + str(current_time))

            folder_chr_results = os.path.join(folder_chr, "results")

            file_output_SNPs = os.path.join(folder_chr, "SNP_that_destroy_motifs.tsv")
            file_chr_results = os.path.join(folder_chr_results, "fimo_adjusted.gff")

            all_GAS_motifs = pd.read_csv(file_chr_results, sep='\t', comment="#", header=None)
            split_columns = all_GAS_motifs[8].str.split(';', expand=True)
            all_GAS_motifs = pd.concat([all_GAS_motifs, split_columns], axis=1)

            #get SNP information
            dbSNP_inputFile_dir = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\01_INITIAL_WORK_Shreeti\\SNP\\SNP files"
            dbSNP_inputFile = os.path.join(dbSNP_inputFile_dir, filename + '.tsv')

            current_time = datetime.now()
            print("Started to read in dbSNP @ " + str(current_time))

            snp_results = []
            snp_results = snp_dict(dbSNP_inputFile)

            current_time = datetime.now()
            print("Read in dbSNP @ " + str(current_time))

            output_data = []

            for index, row in all_GAS_motifs.iterrows():
                start = int(row.iloc[3])
                end = int(row.iloc[4])

                for i in range(start,start+3):
                    if i in snp_results:
                        rs_id = snp_results[i]['RSID']
                        if not any(d['RSID'] == rs_id for d in output_data):
                            output_data.append({
                                'RSID': rs_id,
                                'CHR': snp_results[i]['CHR'],
                                'REF/ALT': snp_results[i]['REF/ALT'],
                                'START': snp_results[i]['START'],
                                'STOP': snp_results[i]['STOP'],
                                'LOCATION': snp_results[i]['LOCATION']
                            })

                for i in range(end-2,end+1):
                    if i in snp_results:
                        rs_id = snp_results[i]['RSID']
                        if not any(d['RSID'] == rs_id for d in output_data):
                            output_data.append({
                                'RSID': rs_id,
                                'CHR': snp_results[i]['CHR'],
                                'REF/ALT': snp_results[i]['REF/ALT'],
                                'START': snp_results[i]['START'],
                                'STOP': snp_results[i]['STOP'],
                                'LOCATION': snp_results[i]['LOCATION']
                            })

            with open(file_output_SNPs, 'w', newline='') as tsv_output_file:
                fieldnames = ['RSID', 'CHR', 'REF/ALT', 'START', 'STOP', 'LOCATION']
                tsv_writer = csv.DictWriter(tsv_output_file, fieldnames=fieldnames, delimiter='\t')
                tsv_writer.writeheader()
                tsv_writer.writerows(output_data)

            print(filename + " END @ " + str(current_time))



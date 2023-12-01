import requests
import mygene
from itertools import product
import pandas as pd
from Bio import Entrez
import time
from datetime import datetime

import re

ENSEMBL_REST_API = "http://rest.ensembl.org"
HEADERS = {'Content-Type': 'application/json'}

def get_gene_info(gene_symbol, email):

    Entrez.email = email

    max_attempts = 100
    attempt = 0
    search_results = None

    while attempt < max_attempts and search_results is None:
        try:
            handle = Entrez.esearch(db="gene", term=f"{gene_symbol}[Gene Name] AND Homo sapiens[Organism]")
            search_results = Entrez.read(handle)
            handle.close()
            break  # Success, so break out of the loop
        except Exception as e:  # Can specify the exact exception if known
            print(f"Attempt {attempt + 1} failed with error: {e} for gene {gene_symbol}")
            attempt += 1
            time.sleep(5)

    # Validate that the search results contain the expected gene symbol
    for gene_id in search_results['IdList']:

        max_attempts_inner = 100
        attempt_inner = 0
        gene_records = None

        while attempt_inner < max_attempts_inner and gene_records is None:
            try:
                handle = Entrez.efetch(db="gene", id=gene_id, retmode="xml")
                gene_records = Entrez.read(handle)
                handle.close()
                break  # Success, so break out of the loop
            except Exception as e:  # Can specify the exact exception if known
                print(f"Attempt {attempt_inner + 1} failed with error: {e}")
                attempt_inner += 1
                time.sleep(1)

        # Check each gene record for the correct gene symbol and extract the required information
        for gene_record in gene_records:
            # Check the Gene-ref for the correct gene symbol
            if gene_record.get('Entrezgene_gene') and gene_record['Entrezgene_gene'].get('Gene-ref'):
                if gene_record['Entrezgene_gene']['Gene-ref'].get('Gene-ref_locus') == gene_symbol:
                    # Navigate through the nested structure safely
                    if gene_record.get('Entrezgene_locus'):
                        for locus in gene_record['Entrezgene_locus']:
                            if locus.get('Gene-commentary_seqs'):
                                for seq in locus['Gene-commentary_seqs']:
                                    if seq.get('Seq-loc_int'):
                                        seq_loc_int = seq['Seq-loc_int']
                                        if seq_loc_int.get('Seq-interval'):
                                            seq_interval = seq_loc_int['Seq-interval']
                                            start = int(seq_interval['Seq-interval_from'])
                                            end = int(seq_interval['Seq-interval_to'])
                                            strand_symbol = seq_interval['Seq-interval_strand']['Na-strand']
                                            strand = -1 if strand_symbol.attributes['value'] == 'minus' else 1
                                            # Return the extracted values
                                            if end > start and strand == -1:
                                                temp = end
                                                end = start
                                                start = temp
                                            return start, end, strand
    # If the correct gene symbol was not found, return None for all
    return None, None, None,

def format_motif_data_with_aggregation(motif_data, chromosome, upstream_start, strand):
    formatted_strings = {}
    for motif, positions in motif_data.items():
        # Initialize an empty list for each motif if not already present
        if motif not in formatted_strings:
            formatted_strings[motif] = []

        # Create formatted strings for each position and add to the list for the motif
        for position in positions:
            Y = position + upstream_start
            Z = Y + len(motif)  # Assuming the length of the motif is to be added
            position_string = f"{chromosome}:{Y}-{Z}:{strand}"
            formatted_strings[motif].append(position_string)

    # Aggregate position strings for each motif
    aggregated_strings = []
    for motif, position_strings in formatted_strings.items():
        aggregated_string = f"{motif}_{','.join(position_strings)}"
        aggregated_strings.append(aggregated_string)

    return "/".join(aggregated_strings)

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

def expand_motif(motif, substitutions):
    """ Recursively expand the motif with the given substitutions. """
    if not motif:
        return ['']

    # Expand the first character if it is in the substitutions
    first_char = motif[0]
    if first_char in substitutions:
        expanded_first_chars = substitutions[first_char]
    else:
        expanded_first_chars = [first_char]

    # Recursively expand the rest of the motif
    expanded_rest = expand_motif(motif[1:], substitutions)

    # Combine the expansions of the first character with the rest
    expanded_motif = []
    for char in expanded_first_chars:
        for suffix in expanded_rest:
            expanded_motif.append(char + suffix)

    return expanded_motif

def generate_sequences(variable, motifs, substitutions):
    """ Generate possible sequences based on the given variable and motif. """
    motif_dict = {
        "T1_GAS": "VTCNNNGAA",
        "T2_GAS": "TVCNNNGAA",
        "C_GAS": "TTDNNNGAA",
        "G_GAS": "TTCNNNHAA",
        "A1_GAS": "TTCNNNGBA",
        "A2_GAS": "TTCNNNGAB"
    }

    # Get the motif corresponding to the variable
    motif = motif_dict.get(variable)
    if not motif:
        return []

    # Expand the motif with substitutions
    return expand_motif(motif, substitutions)

# Define the nucleotide substitutions
substitutions = {
    'V': ['A', 'C', 'G'],
    'D': ['A', 'T', 'G'],
    'H': ['A', 'T', 'C'],
    'B': ['G', 'T', 'C'],
    'N': ['A', 'G', 'T', 'C']
}

def find_all_motif_occurrences(motif_set, sequence):
    found_motifs = {}
    for motif in motif_set:
        # Find all start positions of the motif in the sequence
        start = 0
        positions = []
        while start < len(sequence):
            pos = sequence.find(motif, start)
            if pos == -1:
                break
            positions.append(pos)
            start = pos + 1  # Increment start to find next occurrence
        if positions:
            found_motifs[motif] = positions
    return found_motifs

def get_human_genes_and_sequence(chromosome, start, end):
    """
    Get genes and sequence for a specified region in the human genome.
    """

    region = f"{chromosome}:{start}-{end}"
    sequence_url = f"{ENSEMBL_REST_API}/sequence/region/human/{region}?coord_system_version=GRCh38"
    genes_url = f"{ENSEMBL_REST_API}/overlap/region/human/{region}?feature=gene"

    attempts = 0
    max_attempts = 100
    sleep_time = 5  # time to wait before retrying in seconds
    success = False
    sequence_response = ""

    # Retry connection loop
    while not success and attempts < max_attempts:
        try:
            # Fetch the sequence for the region
            sequence_response = requests.get(sequence_url, headers={'Content-Type': 'text/plain'})
            success = True
        except Exception as e:
            print(f"Sequence Connection attempt {attempts + 1} failed: {e}")
            time.sleep(sleep_time)  # wait before retrying
            attempts += 1

    if sequence_response != "":
        if not sequence_response.ok:
            sequence_response.raise_for_status()
    sequence = sequence_response.text

    attempts_gene = 0
    max_attempts_genes = 100
    success = False
    genes_response = ""

    # Fetch the genes in the region
    while not success and attempts_gene < max_attempts_genes:
        try:
            genes_response = requests.get(genes_url, headers={'Content-Type': 'application/json'})
            success = True
        except Exception as e:
            print(f"Gene Connection attempt {attempts + 1} failed: {e}")
            time.sleep(sleep_time)  # wait before retrying
            attempts_gene += 1
    if genes_response != "":
        if not genes_response.ok:
            genes_response.raise_for_status()

    genes_data = genes_response.json()

    # Extract only the external names of the genes
    gene_names = [gene['external_name'] for gene in genes_data if 'external_name' in gene]

    return sequence, gene_names

def get_mouse_homolog_for_human_gene(human_gene_symbol):
    """
    Get the mouse homolog for a given human gene symbol.
    """

    attempts = 0
    max_attempts = 100
    sleep_time = 5  # time to wait before retrying in seconds
    success = False
    response = ""



    endpoint = f"{ENSEMBL_REST_API}/homology/symbol/human/{human_gene_symbol}?target_species=mouse"

    while not success and attempts < max_attempts:
        try:
            response = requests.get(endpoint, headers=HEADERS)
            success = True
        except Exception as e:
            print(f"Sequence Connection attempt {attempts + 1} failed: {e}")
            time.sleep(sleep_time)  # wait before retrying
            attempts += 1
    if response != "":
        if response.ok:
            data = response.json()

            try:
                ensembl_mouse_ids = [homology['target']['id'] for homology in data['data'][0]['homologies'] if
                                     homology['target']['species'] == 'mus_musculus']
            except Exception as e:
                print('Could not find mgi symbol for hgnc symbol skipping: ' + rs_ID)
                return None

            mg = mygene.MyGeneInfo()
            gene_info = mg.querymany(ensembl_mouse_ids, scopes='ensemblgene', fields='symbol', species='mouse')
            mgi_symbols = {hit.get('symbol') for hit in gene_info if 'symbol' in hit}
            mgi_symbols_list = list(mgi_symbols)
            first_mgi_symbol = mgi_symbols_list[0] if mgi_symbols_list else None

            return first_mgi_symbol
        else:
            print(f"Error fetching homolog for human gene {human_gene_symbol}")
            return None

def get_gene_coordinates_and_sequence(mgi_symbol, kb_upstream_length):
    gene_symbol = mgi_symbol
    upstream_length = kb_upstream_length  # Length of the sequence upstream of the gene to retrieve

    # Ensembl REST API endpoints
    lookup_url = f"http://rest.ensembl.org/lookup/symbol/mus_musculus/{gene_symbol}"
    sequence_url = "http://rest.ensembl.org/sequence/region/mus_musculus/"

    attempts = 0
    max_attempts = 100
    sleep_time = 5  # time to wait before retrying in seconds
    success = False
    lookup_response = ""

    while not success and attempts < max_attempts:
        try:
            # Step 1: Get gene information
            lookup_response = requests.get(lookup_url, headers={'Content-Type': 'application/json'})
            success = True
        except Exception as e:
            print(f"Sequence Connection attempt {attempts + 1} failed: {e}")
            time.sleep(sleep_time)  # wait before retrying
            attempts += 1

    if lookup_response != "":
        if lookup_response.ok:
            gene_info = lookup_response.json()
            chromosome = gene_info['seq_region_name']
            gene_start = gene_info['start']
            gene_end = gene_info['end']
            strand = gene_info['strand']

            # Step 2: Calculate the upstream region start and end based on the strand
            if strand == 1:
                upstream_start = max(1, gene_start - upstream_length)  # Ensure start is not less than 1
                upstream_end = gene_start - 1
            else:  # gene is on the negative strand
                upstream_start = gene_end + 1
                upstream_end = gene_end + upstream_length

            if upstream_start>upstream_end:
                temp_pos = upstream_end
                upstream_end = upstream_start
                upstream_start = temp_pos

            region = f"{chromosome}:{upstream_start}-{upstream_end}"

            sequence_url = f"{ENSEMBL_REST_API}/sequence/region/mouse/{region}?coord_system_version=GRCm38"

            attempts_seq = 0
            max_attempts_seq = 100
            success_seq = False
            sequence_response = ""

            while not success_seq and attempts_seq < max_attempts_seq:
                try:
                    sequence_response = requests.get(sequence_url, headers={'Content-Type': 'text/plain'})
                    success_seq = True
                except Exception as e:
                    print(f"Sequence Connection attempt {attempts_seq + 1} failed: {e}")
                    time.sleep(sleep_time)  # wait before retrying
                    attempts_seq += 1

            if sequence_response != "":
                if sequence_response.ok:
                    upstream_sequence = sequence = sequence_response.text
                    return chromosome, upstream_start, upstream_end, strand, upstream_sequence
            else:
                return None, None, None, None, None
    else:
        return None, None, None, None, None

email = "markus.hoffmann@nih.gov"

input_less10kb_VCF = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\GAS_motifs\\GainOfFunction\\SNPs_in_highAcetylation_AND_immuneGenes_AND_correctSignTarget_CLEANED_upstreamLess10KB.vcf"
output_less10kb_VCF = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\GAS_motifs\\GainOfFunction\\SNPs_in_highAcetylation_AND_immuneGenes_AND_correctSignTarget_CLEANED_upstreamLess10KB_mouseGASfound.vcf"
output_less10kb_mouse_bed = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\GAS_motifs\\GainOfFunction\\SNPs_in_highAcetylation_AND_immuneGenes_AND_correctSignTarget_CLEANED_upstreamLess10KB_mouseGASfound.bed"

input_SNPs = pd.read_csv(input_less10kb_VCF, comment='#', sep='\t', header=0)
input_SNPs.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']

bedFile_df = pd.DataFrame(columns=['chr', 'start', 'end', 'name'])


rows_to_drop = []

counter = 0

for index, row in input_SNPs.iterrows():
    if counter % 10 == 0:
        current_time = datetime.now()
        print(f"Checked {counter} rsID ..." + str(current_time))

    counter = counter + 1

    # Extract the gene name from the INFO column
    #
    info = row['INFO']
    rs_ID = row['ID']
    chromosome = row['CHROM']

    # Split the string into segments based on ';'
    segments = info.split(';')
    # Initialize a dictionary to store gene names and distances for genes with position='upstream'
    upstream_gene_distances = {}

    # Process each segment
    for segment in segments:
        if "position='upstream'" in segment:
            # Extract the gene name
            gene_start = segment.find("gene='") + len("gene='")
            gene_end = segment.find("'", gene_start)
            gene_name = segment[gene_start:gene_end]
            # Extract the distance
            distance_start = segment.find("distance=") + len("distance=")
            distance_end = segment.find(",", distance_start)
            distance = segment[distance_start:distance_end]

            # Store the gene name and distance
            upstream_gene_distances[gene_name] = distance

    for upstream_distance in upstream_gene_distances:

        start_gene, end_gene, strand = get_gene_info(upstream_distance, email)

        if strand == 1:
            start_search = start_gene - int(upstream_gene_distances[upstream_distance]) - 1000
            end_search = start_gene+100
        else:
            start_search = start_gene + int(upstream_gene_distances[upstream_distance]) + 1000
            end_search = start_gene - 100

        gas_type = info.split("-")[0]
        sequence_after_gas = info.split(gas_type + "-")[1].split(";")[0]

        motifs = [
            "VTCNNNGAA",
            "TVCNNNGAA",
            "TTDNNNGAA",
            "TTCNNNHAA",
            "TTCNNNGBA",
            "TTCNNNGAB"
        ]

        found_potential_motifs = False

        if start_search > end_search:
            temp = start_search
            start_search = end_search
            end_search = temp

        # Fetch human genes and sequence
        human_sequence, human_genes = get_human_genes_and_sequence(chromosome, start_search, end_search)

        mouse_gene_details = {}

        for human_gene_symbol in human_genes:
            mouse_ensembl_id = get_mouse_homolog_for_human_gene(human_gene_symbol)
            if mouse_ensembl_id is not None:
                try:
                    chromosome, upstream_start, upstream_end, strand, mouse_sequence = get_gene_coordinates_and_sequence(
                        mouse_ensembl_id, 10000)
                except Exception as e:
                    print('Unkown exception: ' + rs_ID)

                is_in_human = sequence_after_gas in human_sequence

                if not is_in_human:
                    print("PROBLEM OF NOT FINDING PREDICTED GAS OCCURRED WITH: " + rs_ID)

                motif_set = generate_sequences(gas_type, motifs, substitutions)

                mouse_motifs_found = find_all_motif_occurrences(motif_set, mouse_sequence)

                if len(mouse_motifs_found) == 0:
                    continue

                found_potential_motifs = True

                formatted_motif_string = format_motif_data_with_aggregation(mouse_motifs_found, chromosome,
                                                                            upstream_start, strand)

                coordinate_sets = formatted_motif_string.split('/')

                for coordinate in coordinate_sets:
                    # Splitting to extract chr, start, and end
                    parts = coordinate.split('_')
                    chr_and_positions = parts[1].split(':')
                    chr_inside = chr_and_positions[0]
                    start_end = chr_and_positions[1].split('-')
                    start = start_end[0]
                    end = start_end[1]

                    # Adding to DataFrame
                    bedFile_df=bedFile_df._append({'chr': chr_inside, 'start': start, 'end': end, 'name' : rs_ID}, ignore_index=True)

                info = info + ';mouse_GrCm37_' + str(chromosome) + ':' + str(upstream_start) + '-' + str(upstream_end) + ':' + str(strand) + '(' + formatted_motif_string + ")"
            else:
                rows_to_drop.append(index)
                continue
        input_SNPs.at[index, 'INFO'] = info

        if not found_potential_motifs:
            rows_to_drop.append(index)
        print("X")

input_SNPs = input_SNPs.drop(rows_to_drop)

with open(output_less10kb_VCF, 'w') as f:
    f.write("##fileformat=VCFv4.3\n##fileDate=20090805\n##source=myImputationProgramV3.1\n##reference=file:///seq/references/1000GenomesPilot-NCBI36.fastas\n#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA19909\n")
# Write vcf_data_filtered to a file without the index and header
input_SNPs.to_csv(output_less10kb_VCF, sep='\t', index=False, header=False, mode='a')

bedFile_df.to_csv(output_less10kb_mouse_bed, sep='\t', index=False)


# Example usage:
#26826255
#rs_ID = "rs1"
#chromosome = '1'
#start = 26816255
#end = 26836255
#info = "A1_GAS-TTCAAGGTA;(gene='ZDHHC18',position='upstream',distance=432,GeneOntology='immune system process[GO:0002376]');('T', 'T'):418925/('T', 'A'):7/totalNumberStudies:6/totalNumberSamplesStudies:837864/totalNumberNucleotideCovered:418932"
#gas_type = info.split("-")[0]
#sequence_after_gas = info.split(gas_type+"-")[1].split(";")[0]

'''
#get potential GAS motifs
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

'''


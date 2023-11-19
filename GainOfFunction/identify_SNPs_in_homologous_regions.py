import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Function to extract gene names and SNP positions from the VCF file
def extract_gene_names_and_snps(vcf_file):
    snps = {}
    with open(vcf_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            parts = line.split('\t')
            pos = parts[1]
            info_field = parts[7]
            gene_field = next((field for field in info_field.split(';') if field.startswith('gene=')), None)
            if gene_field:
                gene_name = gene_field.split('=')[1]
                snps[pos] = gene_name
    return snps

# Function to find homologous mouse genes using the Ensembl REST API
def find_homologous_mouse_genes(snps):
    homologous_regions = {}
    for pos, human_gene in snps.items():
        response = requests.get(f'http://rest.ensembl.org/homology/symbol/human/{human_gene}?target_species=mouse', headers={ "Content-Type" : "application/json"})
        if response.ok:
            data = response.json()
            for homology in data['data'][0]['homologies']:
                if homology['species'] == 'mus_musculus':
                    mouse_gene_id = homology['target']['id']
                    homologous_regions[pos] = mouse_gene_id
    return homologous_regions

# Function to retrieve mouse genomic sequences
def retrieve_mouse_sequences(homologous_regions):
    sequences = {}
    for pos, mouse_gene_id in homologous_regions.items():
        response = requests.get(f'http://rest.ensembl.org/sequence/id/{mouse_gene_id}?type=genomic', headers={ "Content-Type" : "application/json"})
        if response.ok:
            sequence_data = response.json()
            sequences[pos] = sequence_data['seq']
    return sequences

# Function to filter the VCF file
def filter_vcf_file(vcf_file, sequences):
    with open(vcf_file, 'r') as file, open('filtered_snps.vcf', 'w') as output_file:
        for line in file:
            if line.startswith('#'):
                output_file.write(line)
            else:
                parts = line.split('\t')
                pos = parts[1]
                if pos in sequences:
                    output_file.write(line)

# Path to your VCF file
vcf_file_path = 'C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\GainOfFunction\\SNPs_in_highAcetylation_AND_immuneGenes_AND_correctSignTarget_CLEANED_upstreamLess10KB_noGAS_200bp.vcf'

# Extract SNP information and associated human gene names
snps = extract_gene_names_and_snps(vcf_file_path)

# Find homologous mouse genes for the SNPs
homologous_regions = find_homologous_mouse_genes(snps)

# Retrieve mouse genomic sequences for the identified homologous genes
mouse_sequences = retrieve_mouse_sequences(homologous_regions)

# Filter the original VCF file to only include SNPs with homologous sequences in mouse
filter_vcf_file(vcf_file_path, mouse_sequences)

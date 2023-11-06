from Bio import Entrez
import pandas as pd
import re
import time

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


def determine_position_relative_to_gene(gene_start, gene_end, strand, nucleotide_position):
    # Determine position relative to gene
    if strand == 1:  # Positive strand
        if nucleotide_position < gene_start and nucleotide_position < gene_end:
            return 'upstream', gene_start - nucleotide_position
        elif nucleotide_position > gene_end and nucleotide_position > gene_start:
            return 'downstream', nucleotide_position - gene_end
        else:
            return 'within', 0
    else:  # Negative strand
        if nucleotide_position < gene_start and nucleotide_position < gene_end:
            return 'downstream', gene_end - nucleotide_position
        elif nucleotide_position > gene_start and nucleotide_position > gene_end:
            return 'upstream', nucleotide_position - gene_start
        else:
            return 'within', 0


email = "markus.hoffmann@nih.gov"

output_less10kb_VCF = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\GAS_motifs\\GainOfFunction\\SNPs_in_highAcetylation_AND_immuneGenes_AND_correctSignTarget_CLEANED_upstreamLess10KB.vcf"
output_more10kb_VCF = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\GAS_motifs\\GainOfFunction\\SNPs_in_highAcetylation_AND_immuneGenes_AND_correctSignTarget_CLEANED_upstreamLess10KB.vcf"
output_downstream_VCF = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\GAS_motifs\\GainOfFunction\\SNPs_in_highAcetylation_AND_immuneGenes_AND_correctSignTarget_CLEANED_upstreamLess10KB.vcf"
output_within_VCF = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\GAS_motifs\\GainOfFunction\\SNPs_in_highAcetylation_AND_immuneGenes_AND_correctSignTarget_CLEANED_upstreamLess10KB.vcf"

input_filtered_VCF = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\GAS_motifs\\GainOfFunction\\SNPs_in_highAcetylation_AND_immuneGenes_AND_correctSignTarget_CLEANED.vcf"
input_SNPs = pd.read_csv(input_filtered_VCF, comment='#', sep='\t', header=0)
input_SNPs.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']

columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
df_SNPs_upstream_less10kb = pd.DataFrame(columns=columns)
df_SNPs_upstream_more10kb = pd.DataFrame(columns=columns)
df_SNPs_downstream = pd.DataFrame(columns=columns)
df_SNPs_within = pd.DataFrame(columns=columns)

gene_pattern = re.compile(r"gene='(.*?)'")

counter = 0

for index, row in input_SNPs.iterrows():
    # Extract the gene name from the INFO column
    gene_search = gene_pattern.search(row['INFO'])
    gene_name = gene_search.group(1) if gene_search else 'Unknown'

    # Get the POS value as a variable
    position = row['POS']
    start, end, strand = get_gene_info(gene_name, email)


    if None in (start, end, strand):
        continue

    relative_position, distance = determine_position_relative_to_gene(start, end, strand, position)

    if None in (relative_position, distance):
        continue

    if relative_position == 'upstream':
        if distance <= 10000:
            df_SNPs_upstream_less10kb = df_SNPs_upstream_less10kb._append(row)
        else:
            df_SNPs_upstream_more10kb = df_SNPs_upstream_more10kb._append(row)
    elif relative_position == 'downstream':
        df_SNPs_downstream = df_SNPs_downstream._append(row)
    elif relative_position == 'within':
        df_SNPs_within = df_SNPs_within._append(row)

    counter += 1
    # Check if the counter is divisible by 10
    if counter % 10 == 0:
        print(f"Checked {counter} rsIDs...")

df_SNPs_upstream_less10kb.reset_index(drop=True, inplace=True)
df_SNPs_upstream_more10kb.reset_index(drop=True, inplace=True)
df_SNPs_downstream.reset_index(drop=True, inplace=True)
df_SNPs_within.reset_index(drop=True, inplace=True)

with open(output_less10kb_VCF, 'w') as f:
    f.write("##fileformat=VCFv4.3\n##fileDate=20090805\n##source=myImputationProgramV3.1\n##reference=file:///seq/references/1000GenomesPilot-NCBI36.fastas\n#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA19909\n")
# Write vcf_data_filtered to a file without the index and header
df_SNPs_upstream_less10kb.to_csv(output_less10kb_VCF, sep='\t', index=False, header=False, mode='a')

with open(output_more10kb_VCF, 'w') as f:
    f.write("##fileformat=VCFv4.3\n##fileDate=20090805\n##source=myImputationProgramV3.1\n##reference=file:///seq/references/1000GenomesPilot-NCBI36.fastas\n#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA19909\n")
# Write vcf_data_filtered to a file without the index and header
df_SNPs_upstream_more10kb.to_csv(output_more10kb_VCF, sep='\t', index=False, header=False, mode='a')

with open(output_downstream_VCF, 'w') as f:
    f.write("##fileformat=VCFv4.3\n##fileDate=20090805\n##source=myImputationProgramV3.1\n##reference=file:///seq/references/1000GenomesPilot-NCBI36.fastas\n#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA19909\n")
# Write vcf_data_filtered to a file without the index and header
df_SNPs_downstream.to_csv(output_downstream_VCF, sep='\t', index=False, header=False, mode='a')

with open(output_within_VCF, 'w') as f:
    f.write("##fileformat=VCFv4.3\n##fileDate=20090805\n##source=myImputationProgramV3.1\n##reference=file:///seq/references/1000GenomesPilot-NCBI36.fastas\n#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA19909\n")
# Write vcf_data_filtered to a file without the index and header
df_SNPs_within.to_csv(output_within_VCF, sep='\t', index=False, header=False, mode='a')

print("X")

'''
start, end, strand = get_gene_info("ARID5A", email)
print(start, end, strand)
position = 96594595   # Replace with the actual nucleotide position
relative_position, distance = determine_position_relative_to_gene(start, end, strand, position)
print(relative_position)
'''
from Bio import Entrez
import pandas as pd
import re
import time


def modify_row(row_inner, gene_to_modify, new_position, distance_value):
    # Extract the INFO string from the row
    info_str = row_inner['INFO']
    # Split the INFO field into gene blocks
    gene_blocks = info_str.split(';')

    # Process each gene block
    modified_blocks = []
    for block in gene_blocks:
        # Use a regex pattern to check if this block contains the gene we want to modify
        if re.search(rf"\(gene='{re.escape(gene_to_modify)}'", block):
            # Replace the position value and insert the distance value right after
            pattern = rf"(position=')[^']+(')"
            replacement = f"\\1{new_position}\\2,distance={distance_value}"
            block = re.sub(pattern, replacement, block)
        modified_blocks.append(block)

    # Reassemble the INFO field and update the row
    row_inner['INFO'] = ';'.join(modified_blocks)
    return row_inner

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
output_more10kb_VCF = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\GAS_motifs\\GainOfFunction\\SNPs_in_highAcetylation_AND_immuneGenes_AND_correctSignTarget_CLEANED_upstreamMore10KB.vcf"
output_downstream_VCF = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\GAS_motifs\\GainOfFunction\\SNPs_in_highAcetylation_AND_immuneGenes_AND_correctSignTarget_CLEANED_downstream.vcf"
output_within_VCF = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\GAS_motifs\\GainOfFunction\\SNPs_in_highAcetylation_AND_immuneGenes_AND_correctSignTarget_CLEANED_within.vcf"
output_multiple_VCF = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\GAS_motifs\\GainOfFunction\\SNPs_in_highAcetylation_AND_immuneGenes_AND_correctSignTarget_CLEANED_multipleTargets.vcf"

input_filtered_VCF = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\GAS_motifs\\GainOfFunction\\SNPs_in_highAcetylation_AND_immuneGenes_AND_correctSignTarget_CLEANED.vcf"
input_SNPs = pd.read_csv(input_filtered_VCF, comment='#', sep='\t', header=0)
input_SNPs.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']

columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
df_SNPs_upstream_less10kb = pd.DataFrame(columns=columns)
df_SNPs_upstream_more10kb = pd.DataFrame(columns=columns)
df_SNPs_downstream = pd.DataFrame(columns=columns)
df_SNPs_within = pd.DataFrame(columns=columns)
df_SNPs_multiple = pd.DataFrame(columns=columns)

gene_pattern = re.compile(r"gene='(.*?)'")

counter = 0

for index, row in input_SNPs.iterrows():
    # Extract the gene name from the INFO column
    #
    info_str = row['INFO']
    gene_blocks = info_str.split(';')
    matches = re.findall(gene_pattern, row['INFO'])
    df_SNPs_temp = pd.DataFrame(
        columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'GENE', 'POSITION',
                 'DISTANCE'])

    for gene_block in gene_blocks:
        gene_search = gene_pattern.search(gene_block)

        if gene_search:
            for group in gene_search.groups():
                gene_name = group

                # Get the POS value as a variable
                position = row['POS']
                start, end, strand = get_gene_info(gene_name, email)

                if None in (start, end, strand):
                    continue

                relative_position, distance = determine_position_relative_to_gene(start, end, strand, position)

                if None in (relative_position, distance):
                    continue

                print(str(len(matches)))

                if len(matches) > 1:
                    row_multiple = row
                    row_multiple['GENE'] = gene_name
                    row_multiple['POSITION'] = relative_position
                    row_multiple['DISTANCE'] = distance

                    df_SNPs_temp = df_SNPs_temp._append(row_multiple)
                    print("FOUND MORE THAN ONE IMMUNE GENE FOR ONE SNP")
                    continue

                row_modified = modify_row(row_inner=row, gene_to_modify=gene_name, new_position=relative_position,
                                          distance_value=distance)

                if relative_position == 'upstream':
                    if distance <= 10000:
                        df_SNPs_upstream_less10kb = df_SNPs_upstream_less10kb._append(row_modified)
                    else:
                        df_SNPs_upstream_more10kb = df_SNPs_upstream_more10kb._append(row_modified)
                elif relative_position == 'downstream':
                    df_SNPs_downstream = df_SNPs_downstream._append(row_modified)
                elif relative_position == 'within':
                    df_SNPs_within = df_SNPs_within._append(row_modified)

                counter += 1
                # Check if the counter is divisible by 10
                if counter % 10 == 0:
                    print(f"Checked {counter} gene - rsID combinations...")

    if not df_SNPs_temp.empty:
        df_SNPs_temp.reset_index(drop=True, inplace=True)
        for index_temp, row_temp in df_SNPs_temp.iterrows():
            gene_name = row_temp['GENE']
            relative_position = row_temp['POSITION']
            distance = row_temp['DISTANCE']

            row_modified = modify_row(row_inner=row_temp, gene_to_modify=gene_name, new_position=relative_position,
                                      distance_value=distance)
            df_SNPs_temp['INFO'] = row_modified['INFO']
        df_SNPs_temp = df_SNPs_temp.drop('GENE', axis=1)
        df_SNPs_temp = df_SNPs_temp.drop('POSITION', axis=1)
        df_SNPs_temp = df_SNPs_temp.drop('DISTANCE', axis=1)
        last_row = df_SNPs_temp.iloc[-1]
        df_SNPs_multiple = df_SNPs_multiple._append(last_row)
        df_SNPs_temp = pd.DataFrame(columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT','GENE','POSITION','DISTANCE'])



df_SNPs_upstream_less10kb.reset_index(drop=True, inplace=True)
df_SNPs_upstream_more10kb.reset_index(drop=True, inplace=True)
df_SNPs_downstream.reset_index(drop=True, inplace=True)
df_SNPs_within.reset_index(drop=True, inplace=True)
df_SNPs_multiple.reset_index(drop=True, inplace=True)

with open(output_multiple_VCF, 'w') as f:
    f.write("##fileformat=VCFv4.3\n##fileDate=20090805\n##source=myImputationProgramV3.1\n##reference=file:///seq/references/1000GenomesPilot-NCBI36.fastas\n#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA19909\n")
# Write vcf_data_filtered to a file without the index and header
df_SNPs_multiple.to_csv(output_less10kb_VCF, sep='\t', index=False, header=False, mode='a')

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
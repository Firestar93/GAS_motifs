import requests
import pandas as pd
from time import sleep
import re

def create_info_element(row):
    changes_str = '/'.join([f"{key}:{value}" for key, value in row['Nucleotide Changes'].items()])
    new_info_part = (
        f"{changes_str}/"
        f"totalNumberStudies:{row['Number of Studies']}/"
        f"totalNumberSamplesStudies:{row['Number of Samples']}/"
        f"totalNumberNucleotideCovered:{row['Mutation Samples']}"
    )
    return new_info_part

def parse_gas_pattern(info):
    # This regular expression will match the GAS pattern and the surrounding context
    gas_pattern = re.compile(r"([ATCG]\d*)_GAS")
    match = gas_pattern.search(info)
    if match:
        gas_type = match.group(1)
        if gas_type.startswith('A'):
            target_nucleotide = 'A'
        elif gas_type.startswith('T'):
            target_nucleotide = 'T'
        elif gas_type.startswith('C'):
            target_nucleotide = 'C'
        elif gas_type.startswith('G'):
            target_nucleotide = 'G'
        else:
            return None, None  # Return None if the GAS type is not recognized

        # Extract the number to determine if the change is significant
        number = int(gas_type[1:]) if len(gas_type) > 1 else 0
        return target_nucleotide, number > 1
    return None, None

def count_mutations(allele_annotations):
    mutation_counts = {}

    # Iterate over each allele annotation
    for annotation in allele_annotations:
        # Check for 'frequency' field in annotation
        if 'frequency' in annotation:
            # Iterate over each frequency entry
            for freq in annotation['frequency']:
                # Extract the sequences of interest
                deleted_seq = freq['observation']['deleted_sequence']
                inserted_seq = freq['observation']['inserted_sequence']
                allele_count = freq['allele_count']

                # Use the deleted and inserted sequences as a key in the dictionary
                key = (deleted_seq, inserted_seq)

                # Add the allele count to the total for this mutation
                if key in mutation_counts:
                    mutation_counts[key] += allele_count
                else:
                    mutation_counts[key] = allele_count

    return mutation_counts

def get_snp_info(rsid, info):
    # Extract just the numeric part of the rsID if it starts with "rs"
    numeric_rsid = ''.join(filter(str.isdigit, rsid))
    numeric_rsid = int(numeric_rsid)
    url = f"https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/{numeric_rsid}"

    attempt = 0
    while attempt < 100:
        try:
            response = requests.get(url, timeout=10)  # Set a timeout for the request
            if response.ok:
                data = response.json()

                # Extract allele annotations
                allele_annotations = data['primary_snapshot_data']['allele_annotations']
                mutation_counts = count_mutations(allele_annotations)

                target_nucleotide, change_needed = parse_gas_pattern(info)
                # Check if we need a significant change
                if change_needed:
                    # Filter mutation_counts based on target_nucleotide and the significant change
                    significant_mutations = {
                        mutation: count for mutation, count in mutation_counts.items()
                        if mutation[1] == target_nucleotide and count > 1
                    }

                    # If there are no significant mutations, return None
                    if significant_mutations:
                        return None, None, None, None

                # Initialize containers for allele frequency data
                num_studies = 0
                num_samples = 0
                mutation_samples = 0  # Count of samples with the mutation


                # Iterate over allele_annotations to sum up studies and samples
                for annotation in allele_annotations:
                    for frequency in annotation.get('frequency', []):
                        num_studies += 1
                        num_samples += frequency.get('total_count', 0)
                        mutation_samples += frequency.get('allele_count', 0)



                return mutation_counts, num_studies, num_samples, mutation_samples
            else:
                print(f"Failed to retrieve data for {rsid}")
                return None, None, None
        except requests.Timeout:
            print(f"Timeout occurred for {rsid}, attempt {attempt + 1}/100")
        except requests.RequestException as e:
            print(f"Request failed for {rsid}, attempt {attempt + 1}/100. Error: {e}")
        except Exception as e:
            print(f"An unexpected error occurred for {rsid}, attempt {attempt + 1}/100. Error: {e}")

        sleep(5)
        attempt += 1

    print(f"Failed to retrieve data after 100 attempts for {rsid}")
    return None, None, None, None

vcf_file = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\GAS_motifs\\GainOfFunction\\SNPs_in_highAcetylation_AND_immuneGenes.vcf"
output_SNP_info = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\GAS_motifs\\GainOfFunction\\SNPs_in_highAcetylation_AND_immuneGenes_AND_correctSignTarget.tsv"
output_filtered_VCF = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\GAS_motifs\\GainOfFunction\\SNPs_in_highAcetylation_AND_immuneGenes_AND_correctSignTarget.vcf"

# Read the VCF file, skipping the initial '##' lines
vcf_data = pd.read_csv(vcf_file, comment='#', sep='\t', header=0)
vcf_data.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']

# VCF files typically have a single header line starting with '#'
# The header line is necessary to correctly assign column names
vcf_data.columns = vcf_data.columns.str.lstrip('#')

# Filter rows where 'ID' starts with 'rs'
snp_rows = vcf_data[vcf_data['ID'].astype(str).str.startswith('rs')]

# Initialize an empty list to collect data
snp_info_list = []

# Initialize a counter before the loop starts
counter = 0
# Iterate over the filtered DataFrame
for index, row in snp_rows.iterrows():
    rsid = row['ID']
    info = row['INFO']
    changes, studies, samples, m_samples = get_snp_info(rsid, info)
    if changes is not None:
        snp_info_list.append({
            'rsID': rsid,
            'Nucleotide Changes': changes,
            'Number of Studies': studies,
            'Number of Samples': samples,
            'Mutation Samples': m_samples
        })

    counter += 1
    # Check if the counter is divisible by 10
    if counter % 10 == 0:
        print(f"Checked {counter} rsIDs...")

snp_info_df = pd.DataFrame(snp_info_list)
unique_rsIDs = snp_info_df['rsID'].unique()
vcf_data_filtered = vcf_data[vcf_data['ID'].isin(unique_rsIDs)]

merged_data = vcf_data_filtered.merge(snp_info_df, left_on='ID', right_on='rsID')
merged_data['INFO'] = merged_data.apply(lambda row: row['INFO'] + ';' + create_info_element(row), axis=1)
final_data = merged_data.drop(columns=['Nucleotide Changes', 'Number of Studies', 'Number of Samples', 'Mutation Samples', 'rsID'])
vcf_data_filtered.update(final_data[['ID', 'INFO']])

with open(output_filtered_VCF, 'w') as f:
    f.write("##fileformat=VCFv4.3\n##fileDate=20090805\n##source=myImputationProgramV3.1\n##reference=file:///seq/references/1000GenomesPilot-NCBI36.fastas\n#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA19909\n")
# Write vcf_data_filtered to a file without the index and header
vcf_data_filtered.to_csv(output_filtered_VCF, sep='\t', index=False, header=False, mode='a')


snp_info_df.to_csv(output_SNP_info, sep='\t', index=False)


# After the loop, you can print the total number of rsIDs checked
print(f"Total rsIDs checked: {counter}")

# Convert the list of dicts to a DataFrame
snp_info_df = pd.DataFrame(snp_info_list)

# Save the DataFrame to a CSV file
#snp_info_df.to_csv('snp_info.csv', index=False)

# Optionally, inspect the DataFrame
print(snp_info_df.head())

from Bio import Entrez
from Bio.SeqFeature import FeatureLocation, CompoundLocation
import pymysql
import pandas as pd
import warnings
from datetime import datetime
import time


def fetch_surrounding_genes2(chromosome, position, genome_version):
    # Connect to the UCSC database using PyMySQL
    conn = pymysql.connect(host='genome-mysql.cse.ucsc.edu', user='genome', database=genome_version)
    cursor = conn.cursor()

    # Define the query
    inside_query = f"SELECT name FROM knownGene WHERE chrom='{chromosome}' AND cdsStart <= {position} AND cdsEnd >= {position}"
    upstream_query = f"SELECT name FROM knownGene WHERE chrom='{chromosome}' AND cdsEnd < {position} ORDER BY cdsEnd DESC LIMIT 1"
    downstream_query = f"SELECT name FROM knownGene WHERE chrom='{chromosome}' AND cdsStart > {position} ORDER BY cdsStart ASC LIMIT 1"

    inside_query = f"""
    SELECT k.name AS transcript_id, x.geneSymbol AS ensg_id, k.cdsStart, k.cdsEnd, k.chrom
    FROM knownGene k
    LEFT JOIN kgXref x ON k.name = x.kgID
    WHERE chrom='{chromosome}' AND cdsStart <= {position} AND cdsEnd >= {position};
    """

    upstream_query = f"""
    SELECT k.name AS transcript_id, x.geneSymbol AS ensg_id, k.cdsStart, k.cdsEnd, k.chrom
    FROM knownGene k
    LEFT JOIN kgXref x ON k.name = x.kgID
    WHERE k.chrom='{chromosome}' AND k.cdsEnd < {position}
    ORDER BY k.cdsEnd DESC
    LIMIT 5;
    """

    downstream_query = f"""
    SELECT k.name AS transcript_id, x.geneSymbol AS ensg_id, k.cdsStart, k.cdsEnd, k.chrom
    FROM knownGene k
    LEFT JOIN kgXref x ON k.name = x.kgID
    WHERE k.chrom='{chromosome}' AND k.cdsStart > {position}
    ORDER BY k.cdsStart ASC
    LIMIT 5;
    """
    '''
    query = "SHOW tables;"
    cursor.execute(query)
    results = cursor.fetchall()
    '''

    # Execute the inside query and transform it into a DataFrame
    inside_df = pd.read_sql(inside_query, conn)
    inside_df['chrom'] = chromosome
    inside_df['position'] = 'inside'

    # Execute the upstream query and transform it into a DataFrame
    upstream_df = pd.read_sql(upstream_query, conn)
    upstream_df['chrom'] = chromosome
    upstream_df['position'] = 'upstream'

    # Execute the downstream query and transform it into a DataFrame
    downstream_df = pd.read_sql(downstream_query, conn)
    downstream_df['chrom'] = chromosome
    downstream_df['position'] = 'downstream'

    data = [downstream_df, upstream_df, inside_df]
    nearest_genes_df = pd.concat(data)

    '''
    # Execute the query
    cursor.execute(inside_query)
    results_inside = cursor.fetchall()

    cursor.execute(upstream_query)
    results_upstream = cursor.fetchall()

    cursor.execute(downstream_query)
    results_downstream = cursor.fetchall()


    nearest_genes = {'downstream_gene': results_downstream, 'upstream_gene': results_upstream, 'inside_gene':  results_inside}
    '''
    conn.close()

    # Return surrounding genes
    return nearest_genes_df
    #return results_inside


def fetch_surrounding_genes(chromosome, position, genome_version):
    # Initialize variables for retry logic
    connected = False
    attempts = 0
    max_attempts = 100
    sleep_time = 5  # time to wait before retrying in seconds

    # Retry connection loop
    while not connected and attempts < max_attempts:
        try:
            # Connect to the UCSC database using PyMySQL
            conn = pymysql.connect(host='genome-mysql.cse.ucsc.edu', user='genome', database=genome_version)
            connected = True
        except (pymysql.err.OperationalError, TimeoutError) as e:
            print(f"Connection attempt {attempts + 1} failed: {e}")
            time.sleep(sleep_time)  # wait before retrying
            attempts += 1

    if not connected:
        print("Failed to connect to the database after several attempts.")
        return pd.DataFrame()  # return an empty dataframe

    # The connection was successful, so we can proceed with the queries
    try:
        inside_query = f"""
        SELECT k.name AS transcript_id, x.geneSymbol AS ensg_id, k.cdsStart, k.cdsEnd, k.chrom
        FROM knownGene k
        LEFT JOIN kgXref x ON k.name = x.kgID
        WHERE chrom='{chromosome}' AND cdsStart <= {position} AND cdsEnd >= {position};
        """

        upstream_query = f"""
        SELECT k.name AS transcript_id, x.geneSymbol AS ensg_id, k.cdsStart, k.cdsEnd, k.chrom
        FROM knownGene k
        LEFT JOIN kgXref x ON k.name = x.kgID
        WHERE k.chrom='{chromosome}' AND k.cdsEnd < {position}
        ORDER BY k.cdsEnd DESC
        LIMIT 5;
        """

        downstream_query = f"""
        SELECT k.name AS transcript_id, x.geneSymbol AS ensg_id, k.cdsStart, k.cdsEnd, k.chrom
        FROM knownGene k
        LEFT JOIN kgXref x ON k.name = x.kgID
        WHERE k.chrom='{chromosome}' AND k.cdsStart > {position}
        ORDER BY k.cdsStart ASC
        LIMIT 5;
        """

        # Inside query execution
        inside_df = pd.read_sql(inside_query, conn)
        inside_df['chrom'] = chromosome
        inside_df['position'] = 'inside'

        # Upstream query execution
        upstream_df = pd.read_sql(upstream_query, conn)
        upstream_df['chrom'] = chromosome
        upstream_df['position'] = 'upstream'

        # Downstream query execution
        downstream_df = pd.read_sql(downstream_query, conn)
        downstream_df['chrom'] = chromosome
        downstream_df['position'] = 'downstream'

        # Combine dataframes
        nearest_genes_df = pd.concat([downstream_df, upstream_df, inside_df])

    except Exception as e:
        print(f"An error occurred while executing queries: {e}")
        nearest_genes_df = pd.DataFrame()  # return an empty dataframe in case of error
    finally:
        # Make sure the connection is closed even if an exception occurred
        conn.close()

    return nearest_genes_df

if __name__ == "__main__":

    # Suppress all warnings
    warnings.filterwarnings("ignore")

    output_file = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\GainOfFunction\\SNPs_in_highAcetylation_AND_immuneGenes.vcf"

    immune_genes_df = pd.read_csv("C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\found_SNPs\\immune_genes.tsv", sep="\t")

    allSNPs_that_survive_thisFilter = pd.DataFrame()

    allSNPs_that_survived_acetylationFilter = pd.read_csv("C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\GainOfFunction\\SNPs_in_highAcetylation.vcf", sep='\t', comment='#', header=None)

    genome_version = "hg38"
    for index, row in allSNPs_that_survived_acetylationFilter.iterrows():

        if index%10 == 0:
            current_time = datetime.now()
            print(str(index) + " @ " + str(current_time))

        chromosome = "chr" + row[0]
        position = row[1]

        nearest_genes = fetch_surrounding_genes(chromosome, position, genome_version)
        nearest_genes['transcript_id'] = nearest_genes['transcript_id'].str.replace(r'\.\d+', '', regex=True)

        filtered_df = nearest_genes.merge(immune_genes_df,
                                          left_on=['transcript_id', 'ensg_id'],
                                          right_on=['ensembl_transcript_id', 'external_gene_name'],
                                          how='inner')
        filtered_df.drop(['ensembl_transcript_id', 'ensembl_gene_id'], axis=1, inplace=True)

        if not filtered_df.empty:
            #unique_positions = filtered_df['external_gene_name'].unique()
            #unique_immuneGenes = filtered_df['position'].unique()

            already_added_genes = set()
            info_string = row[7]

            for index2, row2 in filtered_df.iterrows():
                gene = row2["external_gene_name"]
                position = row2["position"]
                go_number = row2["go_id"]
                go_term = row2[len(row2)-1]

                if gene not in already_added_genes:
                    info_string = info_string + ";(gene='" + gene + "',position='"+position+"',GeneOntology='"+go_term+"["+go_number+"]')"
                    already_added_genes.add(gene)

            df_new = pd.DataFrame({
                'CHROM': [row[0]],
                'POS': [row[1]],
                'ID': [row[2]],
                'REF': [row[3]],
                'ALT': [row[4]],
                'QUAL': [row[5]],
                'FILTER': [row[6]],
                'INFO': [info_string],
                'FORMAT': [row[8]],
            })

            data = [allSNPs_that_survive_thisFilter, df_new]
            allSNPs_that_survive_thisFilter = pd.concat(data)

    with open(output_file, 'w') as f:
        f.write(
            "##fileformat=VCFv4.3\n##fileDate=20090805\n##source=myImputationProgramV3.1\n##reference=file:///seq/references/1000GenomesPilot-NCBI36.fastas\n#")

    allSNPs_that_survive_thisFilter.to_csv(output_file, index=False, sep='\t', index_label=False, mode="a")

    print("X")
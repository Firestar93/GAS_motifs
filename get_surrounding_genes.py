from Bio import Entrez
from Bio.SeqFeature import FeatureLocation, CompoundLocation
import mysql.connector


def fetch_surrounding_genes(chromosome, position, genome_version):
    # Connect to the UCSC database
    conn = mysql.connector.connect(host="genome-mysql.cse.ucsc.edu", user="genome", database=genome_version)
    cursor = conn.cursor()

    # Define the query
    query = f"SELECT name, chromStart, chromEnd FROM knownGene WHERE chrom='{chromosome}' AND chromStart <= {position} AND chromEnd >= {position}"

    # Execute the query
    cursor.execute(query)
    results = cursor.fetchall()

    # Close the connection
    conn.close()

    # Return surrounding genes
    return results


if __name__ == "__main__":
    chromosome = "chr1"
    position = 1000000
    genome_version = "hg38"

    genes = fetch_surrounding_genes(chromosome, position, genome_version)
    print("Surrounding genes are:")
    for gene in genes:
        print(gene)

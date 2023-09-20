import pandas as pd

input_bed = "/home/markus/Dropbox/UNI/NIH_Projects/GAS_motifs_SNPs/found_SNPs/bed_files/combined_bed.bed"
input_SNPs = "/home/markus/Dropbox/UNI/NIH_Projects/GAS_motifs_SNPs/found_SNPs/VCF/potential_SNPs.vcf"

bed_file = pd.read_csv(input_bed, sep='\t', header=None)
snps = pd.read_csv(input_SNPs, sep ='\t', comment="#", header=None)

df_importantSNPs = pd.DataFrame(columns=['#CHROM','POS', 'ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','NA19909'])
counter = 0
current_chrom = "hello"

for bed_index, bed_row in bed_file.iterrows():
    chromosome_motif = bed_row[0]
    start_index = int(bed_row[1])
    end_index = int(bed_row[2])

    if current_chrom != chromosome_motif:
        print("CHROMOSOME CHANGE" + str(chromosome_motif))
        current_chrom = chromosome_motif
        chrom_snps = snps.loc[snps[0] == str(chromosome_motif)]
        chrom_snps = chrom_snps.reset_index(drop=True)


    length_motif = end_index - start_index
    if length_motif == 9:
        for snp_index, snp_row in chrom_snps.iterrows():
            chromosome_snp = snp_row[0]
            position_snp = int(snp_row[1])
            rsID = snp_row[2]

            if start_index <= position_snp <= end_index:
                #print("X: " + str(counter) + " " + str(chromosome_snp))
                counter = counter +1
                row_to_add = pd.DataFrame({'#CHROM': [snp_row[0]], 'POS': [snp_row[1]], 'ID': [snp_row[2]], 'REF' : [snp_row[3]], 'ALT' : [snp_row[4]], 'QUAL' : ['.'], 'FILTER' : ['.'], 'INFO': ['.'], 'FORMAT': ['.'], 'NA19909': ['.']})
                df_importantSNPs = pd.concat([df_importantSNPs, row_to_add]).reset_index(drop=True)

df_importantSNPs.to_csv("/home/markus/Dropbox/UNI/NIH_Projects/GAS_motifs_SNPs/found_SNPs/VCF/3gappers.vcf", index=False, sep='\t', index_label=False)
print("X")




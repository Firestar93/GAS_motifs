### REANALYSIS ###

1. AlterChromosomeFasta_for_FIMO.py (makes the chr.fa readable for FIMO)
2. Run FIMO
3. AdjustFIMO_results.py to adjust FIMO (0-based and 1-based bias)
4. IdentifySNPs_from_dbSNP.py to identify destroying SNPs
5. create_VCF.py to create VCF file
6. plot_3AND4gappers_all_AND_perChr.R
7. getGeneClassification.R to get gene list classifications of interests
8. get_surrounding_genes.py to get genes surrounding SNPs in high acetylation and filter SNPs with list from step before
9. 


####OLD CODE #####

Code to analyze GAS motifs after found!

main.py => detection of GAS motifs
motif.py => get bed files for the whole GAS motifs
generate_link.py => generate genome browser link from main.py

createVCF.py => generate VCF file from main.py
createVCF_AFTER_H3K27ac.py => generate VCF file from after using bedtools window
getGenomeBrowserLink_BEDTOOLSWINDOW.py => generate genome browser link from from bedtools window file
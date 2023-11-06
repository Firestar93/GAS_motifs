### REANALYSIS (GAIN OF FUNCTION) ###
1. Run FIMO (run_fimo_automatically.py)
2. merge_and_filter_FIMO.py to get all almost GAS motifs with only one change, name them, and put it into VCF
3. get_create_GAS_motif_SNPs.py get SNPs that create a GAS motif 
4. use genome browser with wiggle cutoff and window to get genes in high acetylation (higher than 20)
5. get_surrounding_genes.py to get genes surrounding SNPs in high acetylation and filter SNPs with list from step before
6. check_dbSNP.py checks dbSNP if the nucleotide to create a GAS motif was observed in at least 2 samples from either the same or two different studies
7. X




### REANALYSIS (LOSS OF FUNCTION AND 4 SPACERS) ###

1. AlterChromosomeFasta_for_FIMO.py (makes the chr.fa readable for FIMO)
2. Run FIMO
3. AdjustFIMO_results.py to adjust FIMO (0-based and 1-based bias)
4. IdentifySNPs_from_dbSNP.py to identify destroying SNPs
5. create_VCF.py to create VCF file
6. plot_3AND4gappers_all_AND_perChr.R
7. getGeneClassification.R to get gene list classifications of interests
8. get_surrounding_genes.py to get genes surrounding SNPs in high acetylation and filter SNPs with list from step before
9. get_highly_induced_targetGenes.R to get highly induced IFN and ILN genes from outside data (Jakub, HK, GEO) and produce normalized counts
10. filter_IFN_induced.py to filter the SNPs
11. plots_GAS_motifs_IFN_induced.R to plot expression
12. getOnly4Spacers.py will give you only 4 spacers to search for gain of function SNPs (i.e., makes a 3 spacer out of a 4 spacer)


####OLD CODE #####

Code to analyze GAS motifs after found!

main.py => detection of GAS motifs
motif.py => get bed files for the whole GAS motifs
generate_link.py => generate genome browser link from main.py

createVCF.py => generate VCF file from main.py
createVCF_AFTER_H3K27ac.py => generate VCF file from after using bedtools window
getGenomeBrowserLink_BEDTOOLSWINDOW.py => generate genome browser link from from bedtools window file
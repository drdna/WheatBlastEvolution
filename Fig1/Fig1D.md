## Calculation of pairwise nucleotide diversity
1. Lineage affilation information was added to the pairwise nucleotide divergences determined by iSNPcaller using the [Pairwise_distant_boxplot.pl](/Fig1/Pairwise_distance_boxplot.pl) script.
'''bash
perl Pairwise_distance_boxplot boxplot.strain.idfile AllSNPCountsJan2021.txt > Pairwise_distances.txt
'''
2. The pairwise distances.txt file was used as input to the [Fig1D_PairwiseDistances.R](/Fig1/Fig1D_PairwiseDistances.R) script. 
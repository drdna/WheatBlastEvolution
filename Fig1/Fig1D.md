## Calculation of pairwise nucleotide diversity
1. Lineage affilation information was added to the pairwise nucleotide divergences determined by iSNPcaller using the [Pairwise_distances_boxplot.pl](/Fig1/Pairwise_distances_boxplot.pl) script.
```bash
perl Pairwise_distances_boxplot.pl boxplot.strain.idfile AllSNPCountsJan2021.txt > Pairwise_distances.txt
```
2. The resulting [Pairwise distances.txt](/Fig1/Pairwise distances.txt) file was used as input to the [Fig1D_PairwiseDistances.R](/Fig1/Fig1D_PairwiseDistances.R) script. 

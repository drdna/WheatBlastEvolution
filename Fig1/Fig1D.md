## Calculation of pairwise nucleotide diversity
1. Lineage affilation information was added to the pairwise nucleotide divergences determined by iSNPcaller using the [Pairwise_distances_boxplot.pl](/Fig1/Pairwise_distances_boxplot.pl) script.
```bash
perl Pairwise_distances_boxplot.pl boxplot.strain.idfile AllSNPCountsJan2021.txt > PairwiseDistances.txt
```
2. The resulting [PairwiseDistances.txt](/Fig1/PairwiseDistances.txt) file was used as input to the [Fig1D_PairwiseDistances.R](/Fig1/Fig1D_PairwiseDistances.R) script to create the output plot: 
 
![/Fig1/Fig1D_PairwiseDistances.png](/Fig1/Fig1D_PairwiseDistances.png). 

## Fig 1C. Neighbor-joining tree
1. OPen terminal wondow #1 and start a new SNP calling project:
```bash
perl iSNPcaller_MT.pl WheatBlast
```
2. Open new terminal window (#2) and copy genomes into GENOMES directory:
```bash
cp RAW_GENOMES/*fasta WheatBlast/GENOMES/
```
3. Type "run" in terminal window # 1 to start genome masking, pairwise alignment and SNP calling pipeline.
4. Generate pairwise distance matrix in MEGA format:
```bash
perl Pairwise_matrix_MEGA.pl WheatBLAST/SNP_COUNTS/SNP_counts_*txt > AllSNPCountsJan2021.txt
```
5. Build neighbor-joining tree using MEGA X.
   
## Fig 1D. Calculation of pairwise nucleotide diversity
1. Lineage affilation information was added to the pairwise nucleotide divergences determined by iSNPcaller using the [Pairwise_distances_boxplot.pl](/Fig1/Pairwise_distances_boxplot.pl) script.
```bash
perl Pairwise_distances_boxplot.pl boxplot.strain.idfile AllSNPCountsJan2021.txt > PairwiseDistances.txt
```
2. The resulting [PairwiseDistances.txt](/Fig1/PairwiseDistances.txt) file was used as input to the [Fig1D_PairwiseDistances.R](/Fig1/Fig1D_PairwiseDistances.R) script to create the output plot: 
 
![/Fig1/Fig1D_PairwiseDistances.png](/Fig1/Fig1D_PairwiseDistances.png). 

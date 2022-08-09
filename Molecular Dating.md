## Measurement of Phylogenetic Signal

The Chr1Chr2Chr5_dated.fasta file comprising all SNPs on the co-inherited regions of chromosomes 1, 2 and 5 was imported into MEGAX and used to build a distance tree with 1,000 bootstrap replications. The resulting tree [Chr1Chr2Chr5_dated.tre](/data/Chr1Chr2Chr5_data.tre) was then used as an input to TempEst v1.5.3 (http://tree.bio.ed.ac.uk/software/tempest/) and the root-to-tip distances were calculated. The data were output as a table, and lineage information was added using the [Root-to-tip_addLineages.pl](/scripts/Root-to-tip_addLineages.pl) script. Correlations between sampling dates and phylogenetic distance were then calculated using a custom R script [Root2Tip.R](/scripts/Root2Tip.R).
![Distances vs Dates](/data/Root2Tip.png)
## Molecular Dating

A carefully curated SNP dataset see [filtering SNP calls](PhylogeneticAnalyses.md#filtering-to-remove-false-snp-calls)

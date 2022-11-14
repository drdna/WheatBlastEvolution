## Measurement of Phylogenetic Signal
The [Chr1Chr2Chr5_dated.fasta](/Fig5/Chr1Chr2Chr5_dated.fasta) file comprising all SNPs on the co-inherited regions of chromosomes 1, 2 and 5 was imported into MEGAX and used to build a distance tree with 1,000 bootstrap replications. The resulting tree [Chr1Chr2Chr5_dated.tre](/Fig5/Chr1Chr2Chr5_data.tre) was then used as an input to TempEst v1.5.3 (http://tree.bio.ed.ac.uk/software/tempest/) and the root-to-tip distances were calculated. The data were output as a table, and lineage information was added using the [Root-to-tip_addLineages.pl](/scripts/Root-to-tip_addLineages.pl) script. Correlations between sampling dates and phylogenetic distance were then calculated using a custom R script [Root2Tip.R](/scripts/Root2Tip.R).
![Distances_vs_Dates](/Fig5/Root2Tip.png)

Figure 1. Root-to-tip distances versus sampling dates.
## Choosing a Nucleotide Substitution Model
A carefully curated SNP dataset was generated by [filtering out false SNP calls](PhylogeneticAnalyses.md#filtering-to-remove-false-snp-calls) and a custom script [ConstantSites.pl](/scripts/ConstantSites.pl) was then used to create a dataset comprising all of the invariant nucleotides from the regions of the B71 reference genome used for SNP calling. Sites were included only if they were present in every samples and did NOT reside in repetitive regions of the genome. The resulting data were exported in phylip format and used as input to [PartitionFinder2](https://github.com/brettc/partitionfinder). The configuration file used for the run was [partition_finder.cfg](/Fig5/partition_finder.cfg) and the model that was selected was HKY + G + I ([best scheme](/Fig5/best_scheme.txt)).
## Molecular Dating
The [Chr1Chr2Chr5_dated.fasta](/Fig5/Chr1Chr2Chr5_dated.fasta) file was imported into BEAUTi and the HKY + G model was selected with an invariant (+I) proportion of 0.985. Runs were performed using both the strict and relaxed clock models (to account for the different substitution rates inferred from the TempEst analysis, Figure 1), with an lower substitution rate limit of 5E-8 and an upper limit of 5E-7. After the xml file was generated upon saving, we manually edited the file to add a constantSiteWeights parameter (A:840276 C: 914905 G: 910177 T: 839104). These values were determined using the [ConstantSites.pl](/scripts/ConstantSites.pl) script. A total 100 million iterations with 10 million burn-in samples were performed for each condition, after which the operator switches were tuned to the settings recommended at runtime completion. This step was repeated until no further recommendations were given, or until the recommended setting simply "flip-flopped" backwards and forwards after each adjustment. Once the optimal paramters had been established, we performed 20 runs for both the strict and relaxed clock settings, with each run comprising 100 million iterations, after 10 million burn-in samplings. Upon completion we combined the logs and the trees for the 5 best runs for both the strict and relaxed clock models (based on the highest ESS values for clock rate and tree height). We then used TreeAnnotator to select the maximum clade credibility tree. The final trees were drawn in FigTree v1.4.4 (http://tree.bio.ed.ac.uk/software/figtree/) and exported in .svg format and then manually edited in Inkscape (https://inkscape.org).

![BEASTtree.tiff](/Fig5/BEASTtree.tiff)
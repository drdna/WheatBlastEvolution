## ChromoPainting of Wheat Blast Genomes
1. Read the iSNPcaller SNP calls to generate a haplotypes file (Po_haplotypes.txt):
```bash
perl Generate_haplotypes.pl B71_SNPs B71_BLASTs Po_haplotypes.txt B71_reference.fasta Chr
```
2. Convert haplotype information into the format required by ChromoPainter - and remove sites with missing data in any strain (note: poor quality genomes were not included in this analysis owing to too many missing datapoints).
```bash
perl Make_chromopainter_Nd_infiles.pl donor_list.txt B71 Poryzae_CP Po_haplotypes.txt none   # can include a comma-separated list of strains to exclude as last argument)
Cull_missing_sites.pl Poryzae_CP
```

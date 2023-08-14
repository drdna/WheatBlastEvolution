# Analysis of secondary contributions from certain donor populations

## 1. Print haplotypes of representative isolates from the PoE1/PoSt populations
1. Extract suspected 2PoE1 & 2PoSt regions from respective ChromoPainter (.cp) file and export as .csv:
```bash
perl Sort_by_populationv3.pl CPAug3121 chr1 36450 1600 | perl 2PoE_2PoSt_haplotype_grps.pl - >2PoE1_Chr1_36450_110.csv
perl Sort_by_populationv3.pl CPAug3121 chr7 1790 2000 | perl 2PoE_2PoSt_haplotype_grps.pl - > 2PoSt_Chr7_1790_2000.csv
```
2. Import csv files into ExtDataFig6.R script.

## 2. Print haplotypes of representative PoX population members
1. Identify the approximate nucleotide position where the PoX contribution on the right arm of chromosome 2 starts:

2. Determine the corresponding column in the ChromoPainter (.cp) file
```bash
awk 'NR==3 {i=0; j=0; while($i > 6250000) i++; print $i; j++; if(j==1) {exit}}' CPAug3121/B71.chr2.V2.complete.cp
```
3. Get an arbitrary end point beyond the presumptive crossover between the 2PoX and PoX contributions:
```bash
awk 'NR==3 {i=0; j=0; while($i > 6400000) i++; print $i; j++; if(j==1) {exit}}' CPAug3121/B71.chr2.V2.complete.cp
```
4. Use these coordinates with the 2PoX_haplotypes.pl script to extract haplotypes (in .fasta format) for this region from select strains:
```bash
perl 2PoX_haplotypes.pl > 2PoX_haplotypes.fasta
```
5. Represent haplotypes in column format for reading into R:
```bash
cat 2PoX_even_more_extended_haplotypes.txt | perl -ne 'chomp($_); ($s = $_) =~ s/>(.+)/$1/ if $_ =~ /^>/; @L = split(//, $_) if $_ !~ /^>/; print "$s\t@L\n" if $_ !~ /^>/' > df_for_2PoX.txt
```
6. Extract corresponding nucleotide positions from ChromoPainter haplotype file:
```bash
awk 'NR == 3 {for(i=66000;i<=71285; i++) print $i}' CPAug3121/B71.chr2.V2.complete.cp > 2PoX_varsites.txt
```
## 3. Convert phylogenetic distances to differences/total number of genomic sites queried (not just variant sites).
This was accomplished by counting non-repetitive DNA positions in the B71v2sh self-aligned alignment string and scaling the tree distances accordingly:
```bash
awk '$1 ~ /Chr1/ {print substr($2, 4738666, 4874894-4738666+1)}' B71v2sh/B71v2sh.B71v2sh_alignments | grep 1 -o | wc -l
```
answer = 136147
```bash
awk '$1 ~ /Chr7/ {print substr($2, 171533, 232826-171533+1)}' B71v2sh/B71v2sh.B71v2sh_alignments | grep 1 -o | wc -l
```
answer = 60716

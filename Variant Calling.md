# Variant Calling
### iSNPcaller
Standard SNP callers fail to detect repetitive portions of fungal genomes which results in the reporting of false variants between non-allelic sequences. Repeat-masking programs also miss large numbers of repeated sequences. To address these problem we built a BLAST-based SNP caller that pre-masks all repeats in the reference and query genomes before alignment, and then performs a second round of masking after alignment to screen out "cryptic" repeats that are not detected using genome self-comparison. Our program, [iSNPcaller](https://github.com/drdna/iSNPcaller) works in an incremental fashion where newly added genomes are first compared with one another, and then with those that have been previously analyzed.
iSNPcaller performs the following operations:
a) Rewrites sequence header lines using a standard format: >genomeID_contig1, >genomeID_contig2, etc.)
b) Creates a repeat-masked version of every genome assembly
c) Blasts each genome against all others in pairwise fashion
d) Determines # of uniquely aligned nucleotide positions for each genome x genome comparision
e) Performs SNP calling and reports: i) total # of SNPs; ii) total number of uniquely aligned nucleotide positions; and iii) SNPs/Mb uniquely aligned sequence
f) Moves analyzed sequences/results into a "PROCESSED" directory, allowing new genomes to be analyzed in an incremental fashion
## SNP calling based on pairwise alignments between all genomes
1. iSNPcaller (multi-threaded version) was used to create a project directory, with subfolders for holding intermadiate analyses and final outputs:
```bash
perl iSNPcaller_MT.pl WheatBlast
```
2. Genomes were copied into the newly-created GENOMES directory:
```bash
cp RAW_GENOMES/*fasta WheatBlast/GENOMES/
```
3. iSNPcaller was then run in multi-threaded mode on a High Performance Computing Cluster using the [iSNPcaller_MT.sh](/scripts/pairwiseVariantCalling/iSNPcaller_MT.sh) SLURM script:
```bash
sbatch $scripts/iSNPcaller_MT.sh WheatBlast
```
## SNPcalling against the B71 reference genome:
1. Each genome assembly was run through a custom script that masks all nucleotide positions that occur in multiple alignments when the genome is BLASTed against itself:
```bash
mkdir MASKED_GENOMEs
cp RAW_GENOMEs/*.fasta MASKED_GENOMEs
perl RMSA_MT.pl MASKED_GENOMEs
```
3. The B71 reference genome was then BLASTed against each of the masked genome assemblies:
```bash
mkdir B71v5_BLAST
cd MASKED_GENOMEs
for f in `ls *masked.fasta`; do blastn -query B71v5_nh_masked.fasta -subject $f -evalue 1e-20 -max_target_seqs 2000 -outfmt '6 qseqid sseqid qstart qend sstart send btop' > ../B71v5_BLAST/B71v5.$f.BLAST; done
```
4. SNPs were then called using the SNPcalling module of iSNPcaller:
```bash
cd ..
perl Run_SU4.pl B71v5_BLAST B71v5_SNPs
```
## SNPcalling against the B71 reference genome using GATK:
SNPs were called using a standard Bowtie2/GATK pipeline using the [BWT2-GATK.sh](/scripts/bowtieGATK/BWT2-GATK.sh) SLURM script. The variant call format file was then filtered to remove: i) sites that occurred in repeat regions of the reference genome (to ensure that all calls were between allelic loci); ii) heterozygous calls (alt:ref ratio < 20; to avoid calling variants between non allelic loci, due to repeat regions in the query genome); and iii) variant calls with low coverage (DP < 10; usually false calls caused by poor sequence quality in homopolymer tracts).

1. The B71 reference genome was indexed using bowtie2-build:
```bash
bowtie2-build B71.fasta B71_index/B71
```
2. Sequence reads were aligned using bowtie2 and genotyping was performed using GATK version 4.1.4.1 using the [BWT2-GATK.sh](/scripts/bowtieGATK/BWT2-GATK.sh) SLURM script.
```bash
for f in `ls FASTQ_DIRECTORY/*_1.fastq.gz | awk -F '/|_' '{print $3}`; do sbatch BWT2-GATK.sh B71.fasta FASTQ_DIRECTORY $f; done
```
## Filtering to remove false SNP calls
3. The "snps-only" VCF files were copied into a new directory and illegal SNP calls were then filtered out using the [SmartSNPsV2.pl](/scripts/bowtieGATK/SmartSNPsV2.pl) script:
```bash
for f in `ls VCF_FILES/*vcf`; do SmartSNPs.pl B71_ALIGN_STRINGs/B71.B71_alignments $f 20 10; done   # alt:ref ratio >= 20; read coverage >= 10
```

## Manual filtering based on comparison between iSNPcaller and BWT/GATK variant datasets
1. The SmartSNPs-filtered data were summarized using the [Summarize_SNPs.pl](/scripts/bowtieGATK/Summarize_SNPs.pl) script. This produces a convenient output format allowing manual inspection of the data to identify possible problems (especially calls in repeated regions that escape repeat detection):
```bash
perl Summarize_SNPs.pl CHR1CHR2CHR5_FINAL > Chr1Ch2Chr5_sites.txt
```
(note: inspection of the resulting output file revealed no obvious problems with the dataset)
3. Next, we used the [GATKviSNPcaller.pl](/scripts/bowtieGATK/GATKviSNPcaller.pl) script to compare variant calls made using the  "genome assembly x reference genome" strategy versus the "reads x reference genome" approach.
```bash
perl GATKviSNPcaller.pl samples.txt Chr1Chr2Chr5_sites.txt B71v5_SNPs > Chr1Chr2Chr5_GATKviSNPs.txt
```
3.   Any differences (either in variant positions and/or which samples possessed a given variant) were investigated to identify the reason for the discrepancy. Confirmed "problem" sites were recorded in a "disallowed-sites" file (for false calls), or in a "add-back" file (for legitimate calls filtered out by the SmartSNPs script). The SNP call dataset was then updated using the [Summarize_SNPs_no_dodgy.pl](/scripts/bowtieGATK/Summarize_SNPs_no_dodgy.pl) script:
```bash
perl Summarize_SNPs_no_dodgy.pl CHR1CHR2CHR5_VCFs > Chr1Chr2Chr5_final.txt
```

# Variant Calling
Standard SNP callers fail to detect repetitive portions of fungal genomes which results in the reporting of false variants between non-allelic sequences. The common repeat-masking programs also miss large numbers of repeated sequences. To address this problem we built a BLAST-based SNP caller that pre-masks  all repeats in the reference and query genomes before alignment, and then performs a second round of masking after alignment to screen out "cryptic" repeats that are not detected using genome self-comparision. Our program, [iSNPcaller_MT.sh](/scripts/iSNPcaller_MT.sh) works in an incremental fashion where newly added genomes are first compared with one another, and are with those that have been previously analyzed.

## SNP calling between all genomes using pairwise alignments
1. The correct directory structure for iSNPcaller was created by running the iSNPdirectory script:
```bash
iSNPdirectory <WheatBlast>
```
2. Fasta files for each genome sequence were then copied into the *GENOMES* directory
3. iSNPcaller was then run in multi-threaded mode using the [iSNPcaller_MT.sh](/scripts/iSNPcaller_MT.sh) SLURM script.

## SNPcalling against the B71 reference genome using iSNPcaller:

1. The fasta headers in each genome assembly were converted to a standard format:
```bash
perl SimpleFastaHeaders_SB.pl RAW_GENOMEs
```
2. Each genome assembly was run through a custom script that masks all nucleotide positions that occur in multiple alignments when the genome is BLASTed against itself:
```bash
mkdir MASKED_GENOMEs
mv RAW_GENOMEs/*.fasta MASKED_GENOMEs
perl RMSA_MT.pl MASKED_GENOMEs
```
3. The B71 reference genome was then BLASTed against each of the masked genome assemblies:
```bash
mkdir B71v5_BLAST
cd MASKED_GENOMEs
for f in `ls MASKED_GENOMEs/*masked.fasta`; do blastn -query B71v5_nh_masked.fasta -subject $f -evalue 1e-20 -max_target_seqs 2000 -outfmt '6 qseqid sseqid qstart qend sstart send btop' > ../B71v5_BLAST/B71v5.$f.BLAST; done
```
4. SNPs were then called using the SNPcalling module of iSNPcaller:
```bash
cd ..
perl Run_SU4.pl B71v5_BLAST B71v5_SNPs
```
## SNPcalling against the B71 reference genome using GATK:
SNPs were called using a standard GATK pipeline. The variant call format file was then filtered to remove: i) sites that occurred in repeat regions of the reference genome (to ensure that all calls were between allelic loci); ii) heterozygous calls (alt:ref ratio < 20; to avoid calling variants between non allelic loci, due to repeat regions in the query genome); and iii) variant calls with low coverage (DP < 10; usually false calls caused by poor sequence quality in homopolymer tracts).

1. The B71 reference genome was indexed using bowtie2-build:
```bash
bowtie2-build B71.fasta B71_index/B71
```
2. Sequence reads were aligned using bowtie2 and genotyping was performed using GATK using the [BWT2-GATK.sh](/scripts/BWT2-GATK.sh) SLURM script.
```bash
for f in `ls FASTQ_DIRECTORY/*_1.fastq.gz | awk -F '/|_' '{print $3}`; do sbatch BWT2-GATK.sh B71.fasta FASTQ_DIRECTORY $f; done
```
## Filtering to remove false SNP calls
3. The "snps-only" VCF files were copied into a new directory and illegal SNP calls were then filtered out using the SmartSNPs.pl script:
```bash
for f in `ls VCF_FILES/*vcf`; do SmartSNPs.pl B71_ALIGN_STRINGs/B71.B71_alignments $f 20 10; done   # alt:ref ratio >= 20; read coverage >= 10
```

## Manual filtering based on comparison between iSNPcaller and BWT/GATK variant datasets
1. The SmartSNPs-filtered data were summarized using the Summarize_SNPs.pl script. This produces a convenient output format allowing manual inspection of the data to identify possible problems (especially calls in repeated regions that escape repeat detection):
```bash
perl Summarize_SNPs.pl CHR1CHR2CHR5_FINAL > Chr1Ch2Chr5_sites.txt
```
(note: inspection of the resulting output file revealed no obvious problems with the dataset)
3. Next, we used the GATKviSNPcaller.pl script to compare variant calls made using the  "genome assembly x reference genome" strategy versus the "reads x reference genome" approach.
```bash
perl GATKviSNPcaller.pl samples.list Chr1Chr2Chr5_sites.txt B71v5_SNPs > Chr1Chr2Chr5_GATKviSNPs.txt
```
3.   Any differences (either in variant positions and/or which samples possessed a given variant) were investigated to identify the reason for the discrepancy. Confirmed "problem" sites were recorded in a "disallowed-sites" file (for false calls), or in a "add-back" file (for legitimate calls filtered out by the SmartSNPs script). The SNP call dataset was then updated using the Summarize_SNPs_no_dodgy.pl script:
```bash
perl Summarize_SNPs_no_dodgy.pl CHR1CHR2CHR5_VCFs > Chr1Chr2Chr5_final.txt
```
4. Finally, the Create_nexus_alignment.pl script was used to generate a sequence alignment file (in nexus format), taking into account the reference base at each position, as well as alignment information (yes/no) across the site in question.
```bash
perl Create_fasta_alignment.pl Chr1Chr2Chr5_final.txt ALIGNSTRINGS
```
This produced the alignment file: Chr1Chr2Chr5_final.fasta.
5. Lastly, date information was appended to each sample identifier using the AddDates.pl script;
```bash
perl AddDates.pl samples.txt Chr1Chr2Chr5_final.fasta > Chr1Chr2Chr5_dated.fasta
```
The resulting file was imported into BEAUTI to gerenate the BEAST input file.

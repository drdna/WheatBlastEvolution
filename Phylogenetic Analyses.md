# Phylogenetic Analysis based on pairwise distance data
Standard SNP callers fail to detect repetitive portions of fungal genomes which results in the reporting of false variants between non-allelic sequences. The common repeat-masking programs also miss large numbers of repeated sequences. To address this problem we built a BLAST-based SNP caller that pre-masks  all repeats in the reference and query genomes before alignment, and then performs a second round of masking after alignment to screen out "cryptic" repeats that are not detected using genome self-comparision. Our program, [iSNPcaller_MT.sh](/scripts/iSNPcaller_MT.sh) works in an incremental fashion where newly added genomes are first compared with one another, and are with those that have been previously analyzed.

## SNPcalling between all genomes using pairwise alignments
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
3. The "snps-only" VCF files were copied into a new directory and illegal SNP calls were then filtered out using the SmartSNPs.pl script:
```bash
for f in `ls VCF_FILES/*vcf`; do SmartSNPs.pl B71_ALIGN_STRINGs/B71.B71_alignments $f 20 10; done   # alt:ref ratio >= 20; read coverage >= 10
```
4. The resulting filtered files were cross-referenced against the iSNPcaller calls using the [BWT2-GATK.pl](/scripts/.pl) script and any sites that were not in perfect correspondence (i.e. same SNP called in same set of strains) were further investigated by examining the raw read data to determine the reason for incongruencies, and the SNP call dataset was corrected as indicated (site rejected, or missed calls added).


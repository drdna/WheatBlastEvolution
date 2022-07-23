# SNP calling pipelines

## Incremental SNP calling using iSNPcaller multi-threaded
This script masks all repeats in the reference and query genomes before alignment and then re-masks "cryptic" repeats that were not detected using genome self-comparision before doing the SNP calling.

1. The correct directory structure for iSNPcaller was created by running the iSNPdirectory script:
```bash
iSNPdirectory <WheatBlast>
```
2. Fasta files for each genome sequence were then copied into the *GENOMES* directory
3. iSNPcaller was then run in multi-threaded mode using the [iSNPcaller_MT.sh](/scripts/iSNPcaller_MT.sh) SLURM script.

## SNPcalling against the B71 reference genome

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
for f in `ls VCF_FILES/*vcf`; do SmartSNPs.pl B71_ALIGN_STRINGs/B71.B71_alignments $f 20 10; done
```
4. The resulting filtered files were cross-referenced against the iSNPcaller calls using the [BWT2-GATK.pl](/scripts/.pl) script and any sites that were not in perfect correspondence (i.e. same SNP called in same set of strains) were further investigated by examining the raw read data to determine the reason for incongruencies, and the SNP calll dataset was corrected as indicated (site rejected, or missed calls added).


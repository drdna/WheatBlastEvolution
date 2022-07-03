# SNP calling pipelines

## Incremental SNP calling using iSNPcaller multi-threaded
This script masks all repeats in the reference and query genomes before alignment and then re-masks "cryptic" repeats that were not detected using genome self-comparision before doing the SNP calling.

1. The correct directory structure for iSNPcaller was created by running the iSNPdirectory script:
```bash
iSNPdirectory <WheatBlast>
```
2. Fasta files for each genome sequence were then copied into the *GENOMES* directory
3. The iSNPcaller was then run using a [SLURM script](/scripts/iSNPcaller_MT.sh).

## SNPcalling against the B71 reference genome

1. The B71 reference genome was indexed using bowtie2-build:
```bash
bowtie2-build B71.fasta B71_index/B71
```
2. Sequence reads were aligned using bowtie2 and genotyping was perforemd using GATK using the BWT2-GATK.sh script


# SNP calling pipelines

## Incremental SNP calling using iSNPcaller multi-threaded

1. The correct directory structure for iSNPcaller was created by running the iSNPdirectory script:
```bash
iSNPdirectory <WheatBlast>
```
2. Fasta files for each genome sequence were then copied into the *GENOMES* directory
3. The iSNPcaller was then run using a [SLURM script] (/scripts)


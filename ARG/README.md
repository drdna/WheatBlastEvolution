# Analysis of Recombinational History using ARGweaver

1. An appropriate dataset was assembled for chromosome 2 using the custom script Generate_ARGsites.pl. ARGweaver struggles with large datasets (even on the supercomputer) and therefore we included a modest set of strains in the analysis, each representing a different chromosomal haplotype: 

### PoL1 members
ATCC64557

BTBa-B1

P28

Po221

Pg1213-22

PtKY18-1

U234

### PoT members:
Br126.1

Br127.11

T1-1

T2-1

T47-3

### Candidate donors:
Bm88324 (PoU1)

Br35 (PoSt)

CD156 (PoE1)

Cd88215 (PoC1)

U168 (PoLu)

U169 (PoE1)

U232 (PoX)

U75 (PoX)

### Non-donors:
EiJA178 (PoE3)

MrJA49 (PoM)

2. ARGweaver was run for 200 iterations using 100 time intervals using a mutation rate of 2E-7 and recombination rate of 1.5E-9.
```bash
## ARGUMENTS

sitesfile=$1

mutrate=$2

recombrate=$3

region=$4

outprefix=${sitesfile/\.sites/}_${mutrate}_${recombrate}_${region}/${sitesfile/\.sites/}_${region}

## RUN

source /project/farman_uksr/miniconda/etc/profile.d/conda.sh

conda activate bioinfo

arg-sample -s $sitesfile --ntimes 100  -n 200 --sample-step 1 -m $mutrate -r $recombrate -o $outprefix --region $region --overwrite

conda deactivate
```

3. An ancestral recombination graph was then build using smc2arg.py script:
```bash
for f in `ls SMC_files/*gz`; do python2 smc2arg $f ${f/gz/arg}; done
```
4. A custom script was then used to iterate through the result .arg files to determine for each candidate donor isolate, the most recent inferred convergence date (time to most recert recombination event, TMRRE) between that isolate and any member of the PoL1/PoT populations
```bash
for f in `ls SMC_files/*arg`; do perl ARGiterator.pl $f >> TMRREs.txt; done
```
5. Plot TMRREs using a custom R script ([TMRREs.R](/ARG/TMRREs.R))
![TMRREs.png](/ARG/TMRREs.png)

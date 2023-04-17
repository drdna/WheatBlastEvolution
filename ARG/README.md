# Analysis of Recombinational History using ARGweaver

## 1. Assemble dataset:
An appropriate dataset was assembled for chromosome 2 using the custom script [Generate_ARGsites.pl](/ARG/Generate_ARGsites.pl). ARGweaver struggles with large datasets (even on the supercomputer) and therefore we included a modest set of strains in the analysis. For PoL1/PoT we included a single representive of each haplotype for this chromosome. We also included a at least one member of each candidate donor population and three putative non-donors.

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

U168 (PoLu)

U169 (PoE1)

U232 (PoX)

U75 (PoX)

### Non-donors:
Cd88215 (PoC1)

EiJA178 (PoE3)

MrJA49 (PoM)

## 2. Run ARGweaver:
ARGweaver was run for 200 iterations using 100 time intervals using a mutation rate of 2E-7 and recombination rate of 1.5E-9.
```bash
## ARGweaver.sh

# Script for running ARGweaver:

# Usage: sbatch ARGweaver.sh <sites-file> <mutation-rate> <recombination-rate> <region-to-analyze>

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

## 3. Build ancestral recombination graphs:
Ancestral recombination graphs were then built using smc2arg.py script:
```bash
for f in `ls SMC_files/*gz`; do python2 smc2arg $f ${f/gz/arg}; done
```
## 4. Determine times to most recent recombination events
A custom script was then used to iterate through the result .arg files to determine for each candidate donor isolate, the most recent inferred convergence date (time to most recent recombination event, TMRRE) between that isolate and any member of the PoL1/PoT population.
```bash
for f in `ls SMC_files/*arg`; do perl ARGiterator.pl $f >> TMRREs.txt; done
```
## 5. Plot TMRREs:
The custom R script ([TMRREs.R](/ARG/TMRREs.R)) was used to generate a plot showing the TMRRE distributions for each candidate donor. Note how the isolates from the main candidate donor populations Bm88324 (PoU1), Br35 (PoSt), U168 (PoLu), U168 (PoE1), U75 and U232 (PoX) all have distributions heavily weighted toward T0, while the others are distributed over a considerable timeframe. Here it should be noted that most isolates show a small number of TMRREs near zero because the small size of non-recombinant chromosome blocks results in several comparision "windows" containing zero SNPs. Additionally, some TMRREs are overestimated because ARGweaver only uses time to convergence to estimate recombination date, and this will be highly dependent on the time to most recent common ancestor between the isolate selected as the candidate donor and the actual donor.
![TMRREs.png](/ARG/TMRREs.png)

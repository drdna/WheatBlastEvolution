# Analysis of Recombinational History using ARGweaver

## 1. Assemble dataset:
The custom script [Generate_ARGsites.pl](/ARG/scripts/Generate_ARGsites.pl) was used to construct an appropriate dataset for chromosome 2 which contains intrigressions from all of the major swarm donors (note: ARGweaver struggles with large datasets - even on the supercomputer - and therefore we included a modest set of strains in the analysis). For PoL1/PoT we included a single representive of each haplotype for chromosome 2. We also included at least one member of each candidate donor population and three putative non-donors.

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
Ancestral recombination graphs were then built using the smc2arg.py script (included in ARGweaver distribution):
```bash
for f in `ls SMC_files/*gz`; do python2 smc2arg $f ${f/gz/arg}; done
```
## 4. Construct tree sequence for select regions of chromosome 2:
a. Use the [ARGweaver.py](/ARG/scripts/ARGweaver.py) script to build maximum clade credibility trees for chromosome positions 0.5, 1, 2, 3, 4, 5, 6, and 7 Mb:
```bash
python ARGweaver.py
```
b. Use [PlotTanglegrams.py](/ARG/scripts/PlotTanglegrams.py) script to plot tree sequnce as a tanglegram:
```bash
python plot-tanglegrams.py sim-trees/
```
c. Use Fig6_ARG.R script to build chromopaintings of chromosome 2 for representative PoL1/PoT members and hypothetical donor isolates

d. Manually merge and edit the resulting pdfs in Illustrator:

![TreeSequence.png](/ARG/tanglegram-ML-trees.png)

## 5. Determine times to most recent recombination events
A custom script ([ARGiteratorNonTL.pl](/ARG/scripts/ARGiteratorNonTL.pl))was then used to iterate through the .arg files to determine for each candidate donor isolate, the most recent inferred convergence date (time to most recent recombination event, TMRRE) between that isolate and any member of the PoL1/PoT population.
```bash
for f in `ls SMC_files/*arg`; do perl ARGiterator.pl $f >> TMRREs.txt; done
```
## 6. Plot TMRREs:
The custom R script ([TMRREs.R](/ARG/scripts/TMRREs.R)) was used to generate a plot showing the TMRRE distributions for each candidate donor. Note how the isolates from the main candidate donor populations Bm88324 (PoU1), Br35 (PoSt), U168 (PoLu), U168 (PoE1), U75 and U232 (PoX) all have distributions heavily weighted toward T0, while the others are distributed over a considerable timeframe. Here it should be noted that most isolates show a small number of TMRREs near zero because the small size of non-recombinant chromosome blocks results in several comparision "windows" containing zero SNPs. Additionally, some TMRREs are overestimated because ARGweaver only uses time to convergence to estimate recombination date, and this will be highly dependent on the time to most recent common ancestor between the isolate selected as the candidate donor and the actual donor.

![TMRREs.png](/ARG/TMRREs.png)

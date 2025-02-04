# Analysis of Recombinational History using ARGweaver

## 1. Assemble dataset:
The custom script [Generate_ARGsites.pl](/ARG/scripts/Generate_ARGsites.pl) was used to construct a downsampled dataset for chromosome 2 ([Chr2ARGvalidSkip0.sites](/ARG/data/Chr2ARGvalidSkip0.sites)) which contains introgressions from all of the major swarm donors (note: ARGweaver struggles with large datasets - even on the supercomputer - and therefore we included a modest set of strains in the analysis). For PoL1/PoT we included a single representive of each haplotype for chromosome 2. We also included at least one member of each candidate donor population and three putative non-donors.

### PoL1 members
ATCC64557, BTBa-B1, P28, Po221, Pg1213-22, PtKY18-1, U234

### PoT members:
Br126.1, Br127.11, T1-1, T2-1, T47-3

### Candidate donors:
Bm88324 (PoU1), Br35 (PoSt), CD156 (PoE1), U168 (PoLu), U169 (PoE1), U232 (PoX), U75 (PoX)

### Non-donors:
Cd88215 (PoC1), EiJA178 (PoE3), MrJA49 (PoM)

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
## 3. Construct tree sequence for select regions of chromosome 2:
a. Download .smc.gz files from MCC supercomputer and use the [ARGweaver.py](/ARG/scripts/ARGweaver.py) script to build maximum clade credibility trees for chromosome positions 0.5, 1, 2, 3, 4, 5, 6, and 7 Mb:
```bash
mkdir smc_trees
scp ${dtnmcc}Chr2ARGsitesSkip4/*smc.gz smc_trees
gunzip smc_trees/*gz
python ARGweaver.py
mkdir sim-trees
mv MCCT*tre sim_trees
```
b. Use [PlotTanglegrams.py](/ARG/scripts/PlotTanglegrams.py) script to plot tree sequence as a tanglegram:
```bash
python PlotTanglegrams.py sim-trees/
```
c. Use [CP4ARG.R](/ARG/scripts/CP4ARG.R) script to build chromopaintings of chromosome 2 for the PoL1/PoT haplotypes represented in the Ancestral Recombination Graph as well as inferred chromopaintings for the donor isolates.

d. Manually merge and edit the resulting pdfs in Illustrator:

![TreeSequence.png](/ARG/tanglegram-ML-trees.png)

## 4. Build ancestral recombination graphs:
Ancestral recombination graphs were then built using the smc2arg.py script (included in ARGweaver distribution):
```bash
for f in `ls SMC_trees/*gz`; do python2 smc2arg $f ${f/gz/arg}; done
```

## 5. Determine times to most recent recombination events (TMRREs)
A custom script ([ARGiteratorNonTL.pl](/ARG/scripts/ARGiteratorNonTL.pl))was then used to iterate through the .arg files to determine for each candidate donor isolate, the most recent inferred convergence date (time to most recent recombination event, TMRRE) between that isolate and any member of the PoL1/PoT population.
```bash
for f in `ls SMC_files/*arg`; do perl ARGiterator.pl $f >> TMRREs.txt; done
```
## 6. Plot TMRREs:
The custom R script ([TMRREs.R](/ARG/scripts/TMRREs.R)) was used to generate a plot showing the TMRRE distributions for each candidate donor, over the 200 MCMC iterations. Note how the isolates from the main candidate donor populations Bm88324 (PoU1), Br35 (PoSt), U168 (PoLu), U168 (PoE1), U75 and U232 (PoX) all have TMRREs of zero, or  heavily weighted towards zero. In contrast non-donors have TMRREs that are distributed over a considerable timeframe. Here it should be noted that all isolates show a small number of TMRREs near zero because the small size of non-recombinant chromosome blocks results in several comparison "windows" containing zero SNPs. Additionally, some TMRREs are larger than expected (e.g. for PoX-related isolates U75 and U232). This is because ARGweaver only uses time to convergence to estimate recombination dates, and this will be highly dependent on the time to most recent common ancestor between the isolate selected as the sampled candidate donor and the actual donor (unsampled).

![TMRREs.png](/ARG/TMRREs.png)

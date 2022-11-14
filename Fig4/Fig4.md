## Fig4
1. Use the [ChromoPaint_to_R_compressed.pl](/Fig4/ChromoPaint_to_R_compressed.pl) script to compress the B71 copyprobsperlocus datasets by removing redundant datapoints. This is necessary to generate a ChromoPainting output file with a manageable size:
```bash
mkdir CPtest
cp CP/B71*copyprobs* CPtest/
for f in `ls CPtest/*out.copyprobs*`; do perl ChromoPaint_to_R_compressed.pl $f; done
```
2. Modify the [Fig4_ChromoPaint.R](/Fig4/Fig4_ChromoPaint.R) script to:
i) point to the HaplotypesTable.xlsx spreadsheet that lists haplotype designations for each strain:
```bash
strains <- read_excel("~/HaplotypesTable.xlsx", sheet = "Sheet1")
```
ii) point to the directory containing the down-sampled ChromoPainter output files:
```bash
file.dir <- "~/CPtest/"
```
3. Run the R script to generate the figure:

![Fig4_ChromoPaintings.tiff](/Fig4/Fig4_ChromoPaintings.tiff)

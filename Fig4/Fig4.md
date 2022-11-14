## Fig4
1. Use the 
2. Modify the [Fig4_ChromoPaint.R](/Fig4/Fig4_ChromoPaint.R) script to:
i) point to the HaplotypesTable.xlsx spreadsheet that lists haplotype designations for each strain:
```bash
strains <- read_excel("~/HaplotypesTable.xlsx", sheet = "Sheet1")
```
ii) point to the directory containing the down-sampled ChromoPainter output files:
```bash
file.dir <- "~/CPtest/"
```
2. Run the script to generate the figure:

![Fig4_ChromoPaintings.tiff](/Fig4/Fig4_ChromoPaintings.tiff)

## Fig3
1. Modify the [Fig3_B71vPoTPoL.R](/Fig3/Fig3_B71vPoTPoL.R) script to: i) point to the strain.idfile in the ChromoPainter input directory:
```bash
cladeName <- read.table("~/strain.idfile", header = FALSE)
```
ii) point to the ShinyHaplotypes output directory (in two places):
```bash
vec <- list.files("~/NEE_SHINY/")
```
```bash
df <- read.table(paste("~/NEE_SHINY/", "Chr1", ".", "ATCC64557", ".diffs", sep = ""), header = TRUE, row.names = 1, check.names = FALSE)3. Run the script ```
2. Run the code to generate the figure:

![B71vPoTPoL.tiff](/Fig3/B71vPoTPoL.tiff)

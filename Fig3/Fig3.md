## Fig3
1. Modify the [Fig3_B71vPoTPoL.R](/Fig3/Fig3_B71vPoTPoL.R) script to:
i) point to the strain.idfile in the ChromoPainter input directory:
```bash
cladeName <- read.table("~/strain.idfile", header = FALSE)
```
ii) point to target files in the ShinyHaplotypes output directory
```bash
df <- read.table(paste("~/NEE_SHINY/",  "Chr", f, ".", "B71", ".diffs", sep = ""), header = TRUE, row.names = 1, check.names = FALSE)
```
2. Run the script to generate the figure:

![B71vPoTPoL.tiff](/Fig3/B71vPoTPoL.tiff)

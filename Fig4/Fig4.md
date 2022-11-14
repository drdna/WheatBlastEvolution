## Fig4
1. Modify the [Fig4_ChromoPaint.R](/Fig4/Fig4_ChromoPaint.R) script to:
i) point to the strain.idfile in the ChromoPainter input directory:
```bash
cladeName <- read.table("~/strain.idfile", header = FALSE)
```
ii) point to target files in the ShinyHaplotypes output directory
```bash
df <- read.table(paste("~/NEE_SHINY/",  "Chr", f, ".", "B71", ".diffs", sep = ""), header = TRUE, row.names = 1, check.names = FALSE)
```
2. Run the script to generate the figure:

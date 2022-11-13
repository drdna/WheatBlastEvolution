## Generation of Figure 2:
1. Use the shell script, Run_CP.sh, to run ChromoPainter, followed by the [ChromoPaint_to_R.pl](/Fig2/ChromoPaint_to_R.pl) script to generate suitable inputs for R. Corresponding "widths" files are also generated that include distances between adjacent SNP sites, so that stacked bars to be plotted without gaps between them:
```bash
./Run_CP.sh
```
2. Fig2A was generated using the [Fig2A_B71_ChromoPaint.R](/Fig2/Fig2A_B71_ChromoPaint.R) script to point to the ChromoPainter "forR" and "widths" output files, [Out_3](/Fig2/Out_3.tar.gz):

![Fig2A.pdf](/Fig2/Fig2A.pdf)

```bash
for f in `ls perl 

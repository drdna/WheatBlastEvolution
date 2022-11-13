## Generation of Figure 2:
1. Use the shell script, Run_CP.sh, to run ChromoPainter, followed by the [ChromoPainter_to_R.pl](/Fig2/ChromoPainter_to_R.pl) script to generate suitable inputs for R. Corresponding "widths" files are also generated that include distances between adjacent SNP sites, so that stacked bars to be plotted without gaps between them:
```bash
./Run_CP.sh
```
2. Fig2A was generated using the [Fig2A_B71_ChromoPaint.R](/Fig2/Fig2A_B71_ChromoPaint.R) script to point to the directory [Out_3](/Fig2/Out_3.tar.gz) that contiain the output files the The copyprobs output file from ChromoPainter was used to generate a suitable input for R, along with a corresponding "widths" file that includes distances between adjacent SNP sites to allow stacked bars to be plotted without gaps between them:
```bash
for f in `ls perl 

## Generation of Figure 2A:
1. Use the shell script, Run_CP.sh, to run ChromoPainter, followed by the [ChromoPaint_to_R.pl](/Fig2/ChromoPaint_to_R.pl) script to generate suitable inputs for R. Corresponding "widths" files are also generated that include distances between adjacent SNP sites, so that stacked bars to be plotted without gaps between them.
```bash
cd Out_3/
./Run_CP.sh
```
2. Fig2A was generated using the [Fig2A_B71_ChromoPaint.R](/Fig2/Fig2A_B71_ChromoPaint.R) script which was run from within directory [Out_3](/Fig2/Out_3.tar.gz) containing the ChromoPainter "forR" and "widths" output files:

![Fig2A.tiff](/Fig2/Fig2A.tiff)

## Generation of Figure 2B:
1. Use the [ShinyHaplotypes.pl](/ShinyHaplotypes/ShinyHaplotypes.pl) script to generate haplotype divergence datasets (use the strain.idfile that was created inside the CP directory). Useage: perl ShinyHaplotypes.pl [strain.idfile] [CPdir] [window_size] [step_size>].
```bash
perl ShinyHaplotypes.pl strain.idfile CP 200 40
```
2. Fig2B was generated by using the [Fig2B_B71vRepresentativeDonors.R](/Fig2/ig2B_B71vRepresentativeDonors.R) script to point to the CP_SHINY directory generated by ShinyHaplotypes.pl (two places in script). The script also needs to read the original strain.idfile to gather lineage affiliation information.
  
![B71vRepresentativeDonors.tiff](/Fig2/B71vRepresentativeDonors.tiff)
  
## Calculate proportion of windows showing zero divergence between B71 and nearest candidate donor isolate
1. Import ShinyHaplotypes .diff files into the [ZeroDivergenceWindows.R](/Fig2/ZeroDivergenceWindows.R) script:
```bash
[1] "Chr1: windows exhibiting zero divergence = 0.83%"
[1] "Chr2: windows exhibiting zero divergence = 0.92%"
[1] "Chr3: windows exhibiting zero divergence = 0.91%"
[1] "Chr4: windows exhibiting zero divergence = 0.52%"
[1] "Chr5: windows exhibiting zero divergence = 0.53%"
[1] "Chr6: windows exhibiting zero divergence = 0.91%"
[1] "Chr7: windows exhibiting zero divergence = 0.89%"
```

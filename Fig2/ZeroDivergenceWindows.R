# Iterate through ShinyHaplotypes diff files and determine proportion of windows 
# showing zero divergence between B71 reference isolate and the best donor(s)

for (f in 1:7) {
  
  print(paste0("Chr", f, sep = ""))
  
  x <- as.data.frame(read.table(paste("~/NEE_SHINY/Chr", f, ".B71.diffs", sep = ""), header = T, row.names = 1))

  y <- x[x$clade != "L",]
  
  z <- y[y$clade != "T",]
  
  windows <- as.numeric(apply(z,1,min))

  nonzero <- windows[as.numeric(apply(z,1,min)) > 0]

  proportion <- 1 - (length(nonzero)/length(windows))

  print(paste0("Non-zero windows = ", proportion, sep = ""))

}

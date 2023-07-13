# Iterate through ShinyHaplotypes diff files and determine proportion of windows 
# showing zero divergence between B71 reference isolate and the best candidate donor(s)

for (f in 1:7) {
  
  x <- as.data.frame(read.table(paste("~/NEE_SHINY/Chr", f, ".B71.diffs", sep = ""), header = T, row.names = 1))

  y <- x[x$clade != "L",]
  
  z <- y[y$clade != "T",]

  # determine minimum value in each column
  windows <- as.numeric(apply(z,1,min))

  # extract columns with non-zero values
  nonzero <- windows[as.numeric(apply(z,1,min)) > 0]

  # calculate proportion of columns with zero as lowest value
  proportionZero <- 1 - (length(nonzero)/length(windows))

  print(paste0("Chr", f, ": windows exhibiting zero divergence = ", round(proportionZero, 2), "%", sep = ""))

}

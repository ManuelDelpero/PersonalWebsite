
getLODs <- function(anovas, phenotype, covariate = "SNP") {
  covars <- unlist(unique(lapply(lapply(lapply(anovas, "[", phenotype), function(x) { return(x[[1]]) }), rownames)))
  cIdx <- which(covars == covariate)
  x <- lapply(lapply(lapply(anovas, "[", phenotype), 
                          function(x) { return(as.numeric(x[[1]][,"Pr(>F)"])) }),"[", cIdx)
  snps <- names(x)
  lodscores <- -log10(as.numeric(x))
  names(lodscores) <- snps
  return(lodscores)
}

plotQTL <- function(anovas.male, anovas.female, phenotype = "Longevity", covariate = "SNP") {
  LODS.males <- getLODs(anovas.male, phenotype)
  LODS.females <- getLODs(anovas.female, phenotype)
  plot(c(0, max(map[, "cumPos"])), c(0, 6), t = 'n', main = phenotype, ylab = "LOD", xlab = "Chromosome", xaxt='n',xaxs='i', yaxs='i')
  for(chr in c(1:19,"X")) {
    onChr <- which(map[, "Chr"] == chr)
    points(map[onChr, "cumPos"], LODS.males[onChr], t = 'l', col = "dodgerblue")
    points(map[onChr, "cumPos"], LODS.females[onChr], t = 'l', col = "darksalmon")
  }
  abline(v = map[which(!diff(as.numeric(map[, "Chr"]) %% 2) == 0), "cumPos"]+10, lty=2)
  mids <- apply(cbind(c(0, map[which(!diff(as.numeric(map[, "Chr"]) %% 2) == 0), "cumPos"]+10), c(map[which(!diff(as.numeric(map[, "Chr"]) %% 2) == 0), "cumPos"]+10, max(as.numeric(map[, "cumPos"])))),1,mean)
  axis(1, at = mids,c(1:19,"X"))
  abline(h = -log10(c(0.1, 0.05, 0.01) / (5*length(length(LODS.males)))) , col=c("red", "orange", "green"))
}

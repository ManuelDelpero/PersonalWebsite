

plotQTL <- function(phenotype = "Longevity", mSNP = 3, fSNP = 1){
  LODS.males <- lapply(lapply(lapply(anovas.male, "[", phenotype), 
               function(x) { return(as.numeric(x[[1]][,"Pr(>F)"])) }),"[", mSNP)
  LODS.females <- lapply(lapply(lapply(anovas.female, "[", phenotype), 
               function(x) { return(as.numeric(x[[1]][,"Pr(>F)"])) }),"[", fSNP)
  plot(c(0, max(map[, "cumPos"])), c(0, 8), t = 'n', main = phenotype, ylab = "LOD", xlab = "Chromosome", xaxt='n',xaxs='i', yaxs='i')
  for(chr in c(1:19,"X")) {
    onChr <- which(map[, "Chr"] == chr)
    points(map[onChr, "cumPos"], -log10(unlist(LODS.males)[onChr]), t = 'l', col = "dodgerblue")
    points(map[onChr, "cumPos"], -log10(unlist(LODS.females)[onChr]), t = 'l', col = "darksalmon")
  }
  abline(v = map[which(!diff(as.numeric(map[, "Chr"]) %% 2) == 0), "cumPos"]+10, lty=2)
  mids <- apply(cbind(c(0, map[which(!diff(as.numeric(map[, "Chr"]) %% 2) == 0), "cumPos"]+10), c(map[which(!diff(as.numeric(map[, "Chr"]) %% 2) == 0), "cumPos"]+10, max(as.numeric(map[, "cumPos"])))),1,mean)
  axis(1, at = mids,c(1:19,"X"))
  abline(h = -log10(c(0.1, 0.05, 0.01) / (5*length(length(LODS.males)))) , col=c("red", "orange", "green"))
}

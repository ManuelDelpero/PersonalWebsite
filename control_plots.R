#
# ITP overview control plot for phenotype data
# (c) Danny Arends, HU-Berlin
#

op <- par(mfrow = c(2,3))
for(phenotype in phenotypes) {
  plot(c(0, length(treatments)), c(0, max(pheno[, phenotype], na.rm = TRUE)), t = 'n', main = phenotype, ylab= "measurement", xlab="Center", xaxt='n', las=2)
  for( sex in as.numeric(names(table(pheno[, "sex"])))){
    i = 1;
    for(treatment in treatments) {
      sidx <- which(pheno[, "Drug_Treatment"] == treatment & pheno[, "sex"] == sex)
      Y = pheno[sidx, phenotype]
      boxplot(Y, at = i + (0.5 * sex) - 1, add = TRUE, notch=TRUE, yaxt='n', col= c("darksalmon", "dodgerblue")[sex+1])
      i <- i + 1
    }
  }
  axis(1, at = -0.5 + 1: length(treatments), treatments, las = 2)
}

plot(c(-0.25, 11.75), c(0, 80), t = 'n', main = "Bodyweight by Center", ylab= "measurement", xlab="Treatment", xaxt='n', las=2, xaxs='i')
i <- 1
for(center in centers) {
  for(phenotype in phenotypes[-5]) {
    for( sex in as.numeric(names(table(pheno[, "sex"])))){
      sidx <- which(pheno[, "Center"] == center & pheno[, "sex"] == sex)
      Y = pheno[sidx, phenotype]
      boxplot(Y, at = i + (0.5 * sex) - 1, add = TRUE, notch=TRUE, yaxt='n', col= c("darksalmon", "dodgerblue")[sex+1])
    }
    i <- i + 1
  }
  abline(v = i - 1.25)
}
axis(1, at = c(2,6,10), centers)

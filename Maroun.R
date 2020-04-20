# Maroun data
library(qtl)
setwd("/home/danny/Github/Maroun")
load("All_Samples_Unique.RData")

itpcross <- calc.genoprob(calluniqueFilled)

setwd("/home/danny/Github/Maroun")

# Format: Write out the phenotypes
pheno <- pull.pheno(itpcross)
row.names(pheno) <- pheno[,"ITP_ID"]
write.table(pheno, "phenotypes.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Format: Write out genotypes
gts <- t(pull.geno(itpcross))
colnames(gts) <- rownames(pheno)
write.table(gts, "genotypes.SNProw.geno", sep = "\t", row.names=TRUE, quote=FALSE)
write.table(t(gts), "genotypes.INDrow.geno", sep = "\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

# Format: Create our map, and write it out
map <- cbind("Chr" = rep(NA, nrow(gts)), 
             Locus = rownames(gts), 
             "cM[F]" = rep(NA, nrow(gts)), 
             "cM[M]" = rep(NA, nrow(gts)), 
             Mb = rep(NA, nrow(gts)))

## - physical map
for(x in 1:nrow(gts)){
  r <- find.markerpos(itpcross, marker=rownames(gts)[x])
  map[x, "Chr"] = r[1, "chr"]
  map[x, "Mb"] = r[1, "pos.female"]
}

## - sex specific cM maps (Haldane)
ssmap <- est.map(itpcross, verbose=TRUE, omit.noninformative=FALSE, offset = 0)

for(x in 1:nrow(gts)) {
  marker <- rownames(gts)[x] # Marker name
  lidx <- which(unlist(lapply(lapply(ssmap,function(x){colnames(x) == marker}), any))) # list index
  map[x, "cM[F]"] = round(as.numeric(ssmap[[lidx]][1, marker]), 2)
  map[x, "cM[M]"] = round(as.numeric(ssmap[[lidx]][2, marker]), 2)
}
write.table(map, "map.txt", sep = "\t", quote = FALSE, row.names=FALSE)









setwd("D:/Ddrive/Collegues/Rob Williams")
write.table(gts, "maroun.num.SNProw.geno", sep = "\t", row.names=FALSE, quote=FALSE)
write.table(t(gts), "maroun.num.INDrow.geno", sep = "\t", row.names=TRUE, col.names=FALSE, quote=FALSE)



cat("@type:4-way\n", file = "maroun.num.geno")
cat("@name:maroun\n", file = "maroun.num.geno", append=TRUE)
cat("@mat:AB\n", file = "maroun.num.geno", append=TRUE)
cat("@pat:CD\n", file = "maroun.num.geno", append=TRUE)
cat("@unk:NA\n", file = "maroun.num.geno", append=TRUE)

write.table(gts, "maroun.num.geno", sep = "\t", row.names=FALSE, quote=FALSE, append = TRUE)

gts <- t(pull.geno(cuniqueFilled))
colnames(gts) <- paste0("I", 1:ncol(gts))

gts[gts == 1] <- "AC"
gts[gts == 2] <- "BC"
gts[gts == 3] <- "AD"
gts[gts == 4] <- "BD"
gts[gts == 5] <- "A"
gts[gts == 6] <- "B"
gts[gts == 7] <- "C"
gts[gts == 8] <- "D"

colnames(gts) <- paste0("I", 1:ncol(gts))

gts <- cbind("Chr" = NA, Locus = rownames(gts), cM = NA, Mb = "", gts)

for(x in 1:nrow(gts)){
  r <- find.markerpos(cuniqueFilled, marker=rownames(gts)[x])
  gts[x, "Chr"] = r[1, "chr"]
  gts[x, "cM"] = r[1, "pos.female"]
}

cat("@type:4-way\n", file = "maroun.txt.geno")
cat("@name:maroun\n", file = "maroun.txt.geno", append=TRUE)
cat("@mat:AB\n", file = "maroun.txt.geno", append=TRUE)
cat("@pat:CD\n", file = "maroun.txt.geno", append=TRUE)
cat("@unk:NA\n", file = "maroun.txt.geno", append=TRUE)

write.table(gts, "maroun.txt.geno", sep = "\t", row.names=FALSE, quote=FALSE, append = TRUE)


pheno <- pull.pheno(cuniqueFilled)
rownames(pheno) <-  paste0("I", 1:nrow(pheno))
cat("ID\t", file = "maroun.pheno")
write.table(pheno, "maroun.pheno", sep = "\t", quote=FALSE,append=TRUE)


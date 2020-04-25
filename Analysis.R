setwd("/home/danny/Github/Maroun")

# - Load in the formatted data
pheno <- read.table("phenotypes.txt", sep = "\t", header = TRUE, row.names = 1)
geno <- read.table("genotypes.INDrow.geno", sep = "\t", header = TRUE)
map <- read.table("map.txt", sep = "\t", header = TRUE, row.names = 2)

## Which animals are female/ male
f <- rownames(pheno)[which(pheno[,"sex"] == 0)]
m <- rownames(pheno)[which(pheno[,"sex"] == 1)]

## E.G. Genotypes of first 10 females or males on X
geno[f[1:10], which(map[,1] == "X")]
geno[m[1:10], which(map[,1] == "X")]

table(unlist(apply(geno[m, which(map[,1] == "X")],1,unlist)))
table(unlist(apply(geno[f, which(map[,1] == "X")],1,unlist)))

# Fill the autosomes
autosomes <- geno[, map[,1] != "X"]
autosomeso <- autosomes
dim(autosomes)

chrX <- geno[, map[,1] == "X"]
chrXo <- chrX
dim(chrX)

autosomes[autosomes == 1] <- "AC"
autosomes[autosomes == 2] <- "BC"
autosomes[autosomes == 3] <- "AD"
autosomes[autosomes == 4] <- "BD"
autosomes[autosomes == 5] <- "A"
autosomes[autosomes == 6] <- "B"
autosomes[autosomes == 7] <- "C"
autosomes[autosomes == 8] <- "D"

#table(unlist(apply(autosomes, 2, unlist)))

chrX[chrX == 3] <- "-"
chrX[chrX == 4] <- "-"
for (x in 1:nrow(chrX)) {
  if (rownames(chrX)[x] %in% m) { # Male
    chrX[x, chrX[x,] == 1] <- "A-"
    chrX[x, chrX[x,] == 2] <- "B-"
  } else {
    chrX[x, chrX[x,] == 1] <- "AC"
    chrX[x, chrX[x,] == 2] <- "BC"
  }
}

table(unlist(apply(chrX, 1, unlist)))

codedgeno <- cbind(autosomes, chrX)[rownames(geno), colnames(geno)]
write.table(codedgeno, file="genotypes.INDrow.code", sep = "\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

genotypes <- read.table("genotypes.INDrow.code", sep = "\t", header = TRUE,na.strings=c("NA", "-"))

#We now can start mapping QTLs

library(lme4)

pheno[, "Batch_Plate"] <- as.character(pheno[, "Batch_Plate"])
pheno[which(is.na(pheno[, "Batch_Plate"])), "Batch_Plate"] <- "Sequenom"

# Figure out the best covariates for each phenotype
lmdata <- data.frame(
            "BW6" = pheno[, "Body_Weight_6m"],
            "BW12" = pheno[, "Body_Weight_12m"],
            "BW18" = pheno[, "Body_Weight_18m"],
            "BW24" = pheno[, "Body_Weight_24m"],
            "Longevity" = pheno[, "Longevity"],
            "Center" = as.factor(pheno[, "Center"]),
            "Sex" = as.factor(pheno[, "sex"]),
            "Cohort_Year" = as.factor(pheno[, "Cohort_Year"]),
            "Birth_Month" = as.factor(substr(pheno[, "Birth_Date"], 6, 7)),
            "Drug_Treatment" = as.factor(pheno[, "Drug_Treatment"]),
            "Batch_Plate" = as.factor(pheno[, "Batch_Plate"])
)

lmdata <- lmdata[-which(apply(lmdata[, 6:11],1, function(x){any(is.na(x))})),]

anova(lm(BW6 ~ Sex + Center + Cohort_Year + Birth_Month  + Drug_Treatment + Batch_Plate, data = lmdata))
anova(lm(BW12 ~ Sex + Center + Cohort_Year + Birth_Month  + Drug_Treatment + Batch_Plate, data = lmdata))
anova(lm(BW18 ~ Sex + Center + Cohort_Year + Birth_Month  + Drug_Treatment + Batch_Plate, data = lmdata))
anova(lm(BW24 ~ Sex + Center + Cohort_Year + Birth_Month  + Drug_Treatment + Batch_Plate, data = lmdata))
anova(lm(Longevity ~ Sex + Center + Cohort_Year + Birth_Month  + Drug_Treatment + Batch_Plate, data = lmdata))
anova(lm(Longevity ~ Sex + Center + Cohort_Year + Drug_Treatment + Batch_Plate + Sex:Center + Sex:Cohort_Year, data = lmdata))

model1 <- lm(Longevity ~ Sex + Center + Cohort_Year + Batch_Plate + Sex:Center + Drug_Treatment, data = lmdata)
anova(lm(Longevity ~ Sex + Batch_Plate + Cohort_Year + Drug_Treatment, data = lmdata))

m1 <- lm(Longevity ~ Sex + Center + Cohort_Year + Birth_Month  + Drug_Treatment + Batch_Plate, data = lmdata)
m2 <- lm(Longevity ~ Sex + Center + Batch_Plate + Cohort_Year + Drug_Treatment, data = lmdata)

AIC(m1,m2)

plot(Longevity ~ Batch_Plate, data = lmdata)

anovas <- apply(genotypes, 2, function(marker){
  lmdata <- cbind(lmdata, SNP = as.factor(marker))
  BW6 <- anova(lm(BW6 ~ Center + Cohort_Year + Birth_Month  + Drug_Treatment + SNP, data = lmdata))
  BW12 <- anova(lm(BW12 ~ Center + Cohort_Year + Birth_Month  + Drug_Treatment + SNP, data = lmdata))
  BW18 <- anova(lm(BW18 ~ Center + Cohort_Year + Birth_Month  + Drug_Treatment + SNP, data = lmdata))
  BW24 <- anova(lm(BW24 ~ Center + Cohort_Year + Birth_Month  + Drug_Treatment + SNP, data = lmdata))
  Longevity <- anova(lm(Longevity ~ Center + Cohort_Year + Birth_Month + Drug_Treatment + SNP, data = lmdata))
  return(list(BW6 = BW6, BW12 = BW12, BW18 = BW18, BW24 = BW24, Longevity = Longevity))
})

BW6 <- lapply(
             lapply(anovas, "[", "BW6"), 
             function(x) { return(as.numeric(x[[1]][1:5, "Pr(>F)"])) })
BW12 <- lapply(
             lapply(anovas, "[", "BW12"), 
             function(x) { return(as.numeric(x[[1]][1:5, "Pr(>F)"])) })
BW18 <- lapply(
             lapply(anovas, "[", "BW18"), 
             function(x) { return(as.numeric(x[[1]][1:5, "Pr(>F)"])) })
BW24 <- lapply(
             lapply(anovas, "[", "BW24"), 
             function(x) { return(as.numeric(x[[1]][1:5, "Pr(>F)"])) })
Longevity <- lapply(
             lapply(anovas, "[", "Longevity"), 
             function(x) { return(as.numeric(x[[1]][1:5, "Pr(>F)"])) })

dnames <- list(rownames(BW6), c("Center", "Cohort_Year","Birth_Month","Drug_Treatment","SNP"))
BW6 <- matrix(unlist(BW6), nrow(map),5, byrow=TRUE, dimnames = dnames)
BW12 <- matrix(unlist(BW12), nrow(map), 5, byrow=TRUE, dimnames = dnames)
BW18 <- matrix(unlist(BW18), nrow(map), 5, byrow=TRUE, dimnames = dnames)
BW24 <- matrix(unlist(BW24), nrow(map), 5, byrow=TRUE, dimnames = dnames)
Longevity <- matrix(unlist(Longevity), nrow(map), 5, byrow=TRUE, dimnames = dnames)

cdiff <- (as.numeric(map[,"Chr"]) %% 2) * 0.1
plot(c(0, nrow(map)), c(0, 27.5), t = 'n')
points(-log10(BW6[,"SNP"]), col = gray(0.2+cdiff) , t = 'b', pch=17)
points(-log10(BW12[,"SNP"]), col = gray(0.4+cdiff), t = 'b', pch=18)
points(-log10(BW18[,"SNP"]), col = gray(0.6+cdiff), t = 'b', pch=19)
points(-log10(BW24[,"SNP"]), col = gray(0.8), t = 'b', pch=20)
points(-log10(Longevity[,"SNP"]), col = 2 + ceiling(cdiff), t = 'b', pch=20)


#
# ITP initial data analysis
# (c) Danny Arends, HU-Berlin
#

setwd("/home/danny/Github/Maroun")

# - Load in the formatted data
pheno <- read.table("phenotypes.txt", sep = "\t", header = TRUE, row.names = 1)
geno <- read.table("genotypes.INDrow.geno", sep = "\t", header = TRUE)
map <- read.table("map.txt", sep = "\t", header = TRUE, row.names = 2)

map <- cbind(map, cumPos = NA)
pPos <- 0
for(chr in c(1:19,"X")){
  map[map[, "Chr"] == chr, "cumPos"] <- map[map[, "Chr"] == chr, "Mb"] + pPos
  pPos <- max(map[map[, "Chr"] == chr, "cumPos"]) + 20
}

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

# Convert the coding into genotypes
autosomes[autosomes == 1] <- "AC"
autosomes[autosomes == 2] <- "BC"
autosomes[autosomes == 3] <- "AD"
autosomes[autosomes == 4] <- "BD"
autosomes[autosomes == 5] <- "A"
autosomes[autosomes == 6] <- "B"
autosomes[autosomes == 7] <- "C"
autosomes[autosomes == 8] <- "D"

sum(autosomes == "A")
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

toSeason <- function(months){
  seasons <- c(rep("Winter", 2), rep("Spring",3), rep("Summer", 3), rep("Autumn", 3), rep("Winter", 1))
  return(seasons[as.numeric(months)])
}

# Perform some phenotype checks
pheno[, "Batch_Plate"] <- as.character(pheno[, "Batch_Plate"])
pheno[, "sex"] <- as.character(pheno[, "sex"])
pheno[pheno[, "sex"] == 0, "sex"] <- "F"
pheno[pheno[, "sex"] == 1, "sex"] <- "M"
pheno[, "sex"] <- as.factor(pheno[, "sex"])
pheno[which(is.na(pheno[, "Batch_Plate"])), "Batch_Plate"] <- "Sequenom"
pheno[, "Center"] <- as.character(pheno[, "Center"])

pheno <- cbind(pheno, "Year" = as.factor(substr(pheno[, "Birth_Date"], 0, 4)))
pheno <- cbind(pheno, "Month" = as.factor(substr(pheno[, "Birth_Date"], 6, 7)))
pheno <- cbind(pheno, "Season" = as.factor(toSeason(substr(pheno[, "Birth_Date"], 6, 7))))
colnames(pheno)[8] <- "Sex"
colnames(pheno)[9] <- "BW6"
colnames(pheno)[10] <- "BW12"
colnames(pheno)[11] <- "BW18"
colnames(pheno)[12] <- "BW24"
pheno <- pheno[, -c(3, 15, 18, 19, 20)]
pheno[1:5,]

phenotypes <- c("BW6","BW12", "BW18", "BW24", "Longevity")
covariates <- c("Sex", "Center", "Year", "Month", "Season", "Drug_Treatment", "Batch_Plate")
treatments <- table(pheno[,"Drug_Treatment"])
treatments <- c("Control", names(treatments[treatments > 150 & names(treatments) != "Control"]))
centers <- names(table(pheno[,"Center"]))

# Fixes due to boxplots: Bodyweight values 0
pheno[which(pheno[, "BW6"] < 10), "BW6"] = NA
pheno[which(pheno[, "BW12"] == 0), "BW12"] = NA
pheno[which(pheno[, "BW18"] == 0), "BW18"] = NA
pheno[which(pheno[, "BW24"] == 0), "BW24"] = NA

source("control_plots.R")

# Figure out the best covariates for each phenotype
lmdata <- data.frame(
            "BW6" = pheno[, "BW6"],
            "BW12" = pheno[, "BW12"],
            "BW18" = pheno[, "BW18"],
            "BW24" = pheno[, "BW24"],
            "Longevity" = pheno[, "Longevity"],
            "Center" = as.factor(pheno[, "Center"]),
            "Sex" = as.factor(pheno[, "Sex"]),
            "Year" = as.factor(pheno[, "Year"]),
            "Month" = as.factor(pheno[, "Month"]),
            "Season" = as.factor(pheno[, "Season"]),
            "Drug_Treatment" = as.factor(pheno[, "Drug_Treatment"]),
            "Batch_Plate" = as.factor(pheno[, "Batch_Plate"])
)
rownames(lmdata) <- rownames(pheno)

lmdata <- lmdata[-which(apply(lmdata[, 6:12], 1, function(x){any(is.na(x))})),]
males.ctrl <- lmdata[which(lmdata[, "Drug_Treatment"] == "Control" & lmdata[, "Sex"] == "M"),]
females.ctrl <- lmdata[which(lmdata[, "Drug_Treatment"] == "Control" & lmdata[, "Sex"] == "F"),]

# Minimal models for male control samples
anova(lm(BW6 ~ Center + Year, data = males.ctrl))
anova(lm(BW12 ~ Center + Year, data = males.ctrl))
anova(lm(BW18 ~ Center + Year, data = males.ctrl))
anova(lm(BW24 ~ Center, data = males.ctrl))
anova(lm(Longevity ~ Center + Year, data = males.ctrl))

# Minimal models for female control samples
anova(lm(BW6 ~ Center + Year + Season, data = females.ctrl))
anova(lm(BW12 ~ Center + Year + Season,, data = females.ctrl))
anova(lm(BW18 ~ Center + Year + Season, data = females.ctrl))
anova(lm(BW24 ~ Center + Year + Season, data = females.ctrl))
anova(lm(Longevity ~ 1, data = females.ctrl))

# We now can start mapping QTLs in control males and females
library(lme4)
anovas.male <- apply(genotypes[rownames(males.ctrl), ], 2, function(marker) {
  mdata <- cbind(males.ctrl, SNP = as.factor(marker))
  BW6 <- anova(lm(BW6 ~ Center + Year + SNP, data = mdata))
  BW12 <- anova(lm(BW12 ~ Center + Year + SNP, data = mdata))
  BW18 <- anova(lm(BW18 ~ Center + Year + SNP, data = mdata))
  BW24 <- anova(lm(BW24 ~ Center + SNP, data = mdata))
  Longevity <- anova(lm(Longevity ~ Center + Year + SNP, data = mdata))
  return(list(BW6 = BW6, BW12 = BW12, BW18 = BW18, BW24 = BW24, Longevity = Longevity))
})

anovas.female <- apply(genotypes[rownames(females.ctrl), ], 2, function(marker) {
  mdata <- cbind(females.ctrl, SNP = as.factor(marker))
  BW6 <- anova(lm(BW6 ~ Center + Year + Season + SNP, data = mdata))
  BW12 <- anova(lm(BW12 ~ Center + Year + Season + SNP, data = mdata))
  BW18 <- anova(lm(BW18 ~ Center + Year + Season + SNP, data = mdata))
  BW24 <- anova(lm(BW24 ~ Center + Year + Season + SNP, data = mdata))
  Longevity <- anova(lm(Longevity ~ SNP, data = mdata))
  return(list(BW6 = BW6, BW12 = BW12, BW18 = BW18, BW24 = BW24, Longevity = Longevity))
})

source("plot.R") # Plot functions, depends on anovas.male & anovas.female being available
# Quick visualization of QTLs in males and females (different model length for males / females and phenotype)
plotQTL("BW6", 3, 4)
plotQTL("BW12", 3, 4)
plotQTL("BW18", 3, 4)
plotQTL("BW24", 2, 4)
plotQTL("Longevity", 3, 1)

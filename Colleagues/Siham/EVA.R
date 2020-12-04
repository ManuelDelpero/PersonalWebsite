# Data processing for submission to the European Variation Archive (EVA)
#
# copyright (c) - Manuel Delpero
# first written december, 2020
# 

setwd("C:/Users/Manuel/Desktop/Siham")

genotypesRAW <- read.csv("RawSnps.txt", header = TRUE, check.names = FALSE, sep="\t", colClasses="character", na.strings=c("","NA"))
genotypesNEW <- read.csv("NewSNPs.txt", header = TRUE, check.names = FALSE, sep="\t", colClasses="character", na.strings=c("","NA"))


genotypesNEW <- genotypesNEW[,c(1,2,3,4)]

combined <- genotypesRAW[which(genotypesRAW[,2] %in% genotypesNEW[,1]),]

genotypesNEW <- cbind(genotypesNEW, combined)

# Convert in EVA format
for(x in 1:nrow(genotypesNEW)){
  ref <- as.character(genotypesNEW[x,"REF"])
  alt <- as.character(genotypesNEW[x,"ALT"])
  samples <- genotypesNEW[1,c(9:length(colnames(genotypesNEW)))]
  HOMref <- paste0(ref,ref)
  HOMalt <- paste0(alt,alt)
  HET1 <- paste0(ref,alt)
  HET2 <- paste0(alt,ref)
  genotypesNEW[x,grep(HOMref, genotypesNEW[x,])] = "0/0"
  genotypesNEW[x,grep(HOMalt, genotypesNEW[x,])] = "1/1"
  genotypesNEW[x,grep(HET1, genotypesNEW[x,])] = "0/1"
  genotypesNEW[x,grep(HET2, genotypesNEW[x,])] = "0/1"
  genotypesNEW[x, is.na(genotypesNEW[x,])] = "./."
}

write.table(genotypesNEW, file = "EVA_NEW_SNPs.txt", row.names = FALSE, sep = "\t", quote = FALSE)



genotypesNEW_names <- read.csv("FilteredSNPs_Genotypesname.txt", header = TRUE, check.names = FALSE, sep="\t", colClasses="character", na.strings=c("","NA"))

genotypesNEW_names <- genotypesNEW_names[which(genotypesNEW_names[,"POS"] %in% genotypesNEW[,"POS"]),]

columns <- colnames(genotypesNEW_names)[grep("1", colnames(genotypesNEW_names))]
columns <- gsub("1", "", columns)

# Calculate total and for breed allele frequencies 
tots <- c()
AllFreq <- matrix(NA, nrow(genotypesNEW_names), length(columns), dimnames=list(genotypesNEW_names[,"POS"], columns))
for(r in 1:nrow(genotypesNEW_names)){
  alt <- as.character(genotypesNEW_names[r,"ALT"])
  ref <- as.character(genotypesNEW[r,"REF"])
  HOMref <- paste0(ref,ref)
  HOMalt <- paste0(alt,alt)
  HET1 <- paste0(ref,alt)
  HET2 <- paste0(alt,ref)
  for (c in columns){
    breed <- genotypesNEW_names[r,grep(c, colnames(genotypesNEW_names[r,]))] 
	ALTtotbreed <- 2*length(grep(HOMalt, breed)) + length(grep(HET1, breed)) + length(grep(HET2, breed))
	ALLtotbreed  <- 2*length(grep(HOMalt, breed)) + 2*length(grep(HET1, breed)) + 2*length(grep(HET2, breed)) + 2*length(grep(HOMref, breed))
	FRalt <- ALTtotbreed/ALLtotbreed
	AllFreq[r,c] <- FRalt
  }
  ALTtot <- 2*length(grep(HOMalt, genotypesNEW_names[r,])) + length(grep(HET1, genotypesNEW_names[r,])) + length(grep(HET2, genotypesNEW_names[r,]))
  ALLtot  <- 2*length(grep(HOMalt, genotypesNEW_names[r,])) + 2*length(grep(HET1, genotypesNEW_names[r,])) + 2*length(grep(HET2, genotypesNEW_names[r,])) + 2*length(grep(HOMref, genotypesNEW_names[r,]))
  tot <- ALTtot/ALLtot
  tots <- c(tots, tot)
}

write.table(AllFreq, file = "EVA_NEW_FREQ_SNPs.txt", row.names = TRUE, sep = "\t", quote = FALSE)

AllFreqe <- matrix(NA, nrow(genotypesNEW_names), 1, dimnames=list(genotypesNEW_names[,"POS"], "Frequencies"))
for(r in 1:nrow(genotypesNEW_names)){
  Freq <- paste0("AFtotal=", tots[r], ";AFSaanen=", AllFreq[r,1], ";AFNbian=", AllFreq[r,2], ";AFDesert=", AllFreq[r,3], ";AFNilotic=", AllFreq[r,4], ";AFTaggar=", AllFreq[r,5], ";AFNubianibex=", AllFreq[r,6], ";AFBezoaribex=", AllFreq[r,7], ";AFAlpineibex=", AllFreq[r,8])
  AllFreqe[r,1] <- Freq
}  

write.table(AllFreqe, file = "EVA_NEW_FREQ_SNPs_Format.txt", row.names = TRUE, sep = "\t", quote = FALSE)
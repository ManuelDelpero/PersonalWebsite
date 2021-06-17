setwd("C:/Users/Manuel/Desktop/Siham")

data <- read.table("DGAT1.txt", sep = '\t', header=TRUE, colClasses = "character")
head(data)

# Get all the SAA animals
saa <- data[,grep("SAA", colnames(data))]

# Contry
country <- c("FRCH", "CHCH", "TZCH")


# clean the dataset
for (row in 1:nrow(saa)){
  for (col in 1:ncol(saa)){
    saa[row, col] <- paste0(unlist(strsplit(as.character(saa[row, col]), ""))[c(1,3)], collapse = "")
	if (saa[row, col] == ".."){
	  saa[row, col] = NA
	}
  }
}
	
# Divide SAA by country and calculate allele freq for each country
fr <- saa[,grep("FRCH", colnames(saa))]
ch <- saa[,grep("CHCH", colnames(saa))]
tz <- saa[,grep("TZCH", colnames(saa))]

allFreqFR <- matrix(NA, nrow = length(rownames(data)), ncol = 3, dimnames = list(rownames(data),c("alt", "ref", "total_animals")))
for (row in 1:nrow(fr)){
  tot <- (length(fr[row,]) - length(which(is.na(fr[row,])))) *2
  if ("01" %in% rownames(apply(fr[row,],1,  table))){
    het <- apply(fr[row,],1,  table)["01",]
	alt <- apply(fr[row,],1,  table)["01",]/tot
  }else{ 
    het = 0
	alt = 0
	ref = 1
  }
  if ("11" %in% rownames(apply(fr[row,],1,  table))){
    alt <- ((apply(fr[row,],1,  table)["11",] * 2) + het)/tot
  }else{ alt = alt}
  if ("00" %in% rownames(apply(fr[row,],1,  table))){
    ref <- ((apply(fr[row,],1,  table)["00",] * 2) + het)/tot
  }else{ ref = het/tot}
  if (length(apply(fr[row,],1,  table)) == 1){
    ref = 1
  }
  tot <- tot/2
  allFreqFR[row, 1] <- alt
  allFreqFR[row, 2] <- ref
  allFreqFR[row, 3] <- tot
}

allFreqCH <- matrix(NA, nrow = length(rownames(data)), ncol = 3, dimnames = list(rownames(data),c("alt", "ref", "total_animals")))
for (row in 1:nrow(ch)){
  tot <- (length(ch[row,]) - length(which(is.na(ch[row,])))) *2
  if ("01" %in% rownames(apply(ch[row,],1,  table))){
    het <- apply(ch[row,],1,  table)["01",]
	alt <- apply(ch[row,],1,  table)["01",]/tot
  }else{ 
    het = 0
	alt = 0
	ref = 1
  }
  if ("11" %in% rownames(apply(ch[row,],1,  table))){
    alt <- ((apply(ch[row,],1,  table)["11",] * 2) + het)/tot
  }else{ alt = alt}
  if ("00" %in% rownames(apply(ch[row,],1,  table))){
    ref <- ((apply(ch[row,],1,  table)["00",] * 2) + het)/tot
  }else{ ref = het/tot}
  if (length(apply(ch[row,],1,  table)) == 1){
    ref = 1
  }
  tot <- tot/2
  allFreqCH[row, 1] <- alt
  allFreqCH[row, 2] <- ref
  allFreqCH[row, 3] <- tot
}

allFreqTZ <- matrix(NA, nrow = length(rownames(data)), ncol = 3, dimnames = list(rownames(data),c("alt", "ref", "total_animals")))
for (row in 1:nrow(tz)){
  tot <- (length(tz[row,]) - length(which(is.na(tz[row,])))) *2
  if ("01" %in% rownames(apply(tz[row,],1,  table))){
    het <- apply(tz[row,],1,  table)["01",]
	alt <- apply(tz[row,],1,  table)["01",]/tot
  }else{ 
    het = 0
	alt = 0
	ref = 1
  }
  if ("11" %in% rownames(apply(tz[row,],1,  table))){
    alt <- ((apply(tz[row,],1,  table)["11",] * 2) + het)/tot
  }else{ alt = alt}
  if ("00" %in% rownames(apply(tz[row,],1,  table))){
    ref <- ((apply(tz[row,],1,  table)["00",] * 2) + het)/tot
  }else{ ref = het/tot}
  if (length(apply(tz[row,],1,  table)) == 1){
    ref = 1
  }
  tot <- tot/2
  allFreqTZ[row, 1] <- alt
  allFreqTZ[row, 2] <- ref
  allFreqTZ[row, 3] <- tot
}

allFreqFR <- cbind(data[,c(1:2)], allFreqFR)
allFreqCH <- cbind(data[,c(1:2)], allFreqCH)
allFreqTZ <- cbind(data[,c(1:2)], allFreqTZ)

write.table(allFreqFR, file = "SAAFR.txt", sep = "\t", quote = FALSE)
write.table(allFreqCH, file = "SAACH.txt", sep = "\t", quote = FALSE)
write.table(allFreqTZ, file = "SAATZ.txt", sep = "\t", quote = FALSE)

  
# objective: This program will be ran by a script to combine two sets of data
# # Argument inputs
Args <- commandArgs(trailingOnly=TRUE)
# #input variables will need to be a number for the seed and all input files. 
CopyNumber <- read.csv(Args[1], sep = ",", header = TRUE)
Alkaloids <- read.csv(Args[2], sep = ",", header = TRUE)

# # Formatting alklaoid data
samplenum <- nrow(Alkaloids)
for (row in samplenum:1){
  current <- Alkaloids$Agrinostics.ID[row]
  if (is.na(current) ){
    Alkaloids <- Alkaloids[-c(row),]
    print(paste0("destroying row ", row))
  }
}

# # Naming samples the correct name in Alkaloids data
samplenum <- nrow(Alkaloids)
for (row in 1:samplenum){
  current <- Alkaloids$X[row]
  if (is.na(current) ){
    Alkaloids$Sample[row] <- paste0(Alkaloids$UGA..[row])
  } else {
    Alkaloids$Sample[row] <- paste0(Alkaloids$UGA..[row],"-",Alkaloids$X[row],"-",Alkaloids$X.1[row], " ")
  }
}

# # Renaming columns and combining into one data set
names(CopyNumber)[names(CopyNumber) == "Treatment"] <- "Sample"
Alkaloids2 <- subset(Alkaloids, select=c("Sample","Mean.Alk..ng.g"))

CopyNumber2 <- subset(CopyNumber, select=(-c(X)))
for (i in 1:nrow(CopyNumber2)){
  a <- CopyNumber2$Sample[i]
  CopyNumber2$Sample[i] = iconv(a, from = 'UTF-8', to = 'ASCII//TRANSLIT')
}

bothdata <- merge(CopyNumber2,Alkaloids2, by.x = "Sample", by.y = "Sample")
head(bothdata)
write.table(bothdata, file = "comboData.csv",
            sep = ",", row.names = FALSE)

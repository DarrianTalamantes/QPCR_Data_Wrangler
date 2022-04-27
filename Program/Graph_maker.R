#!/usr/bin/env Rscript
library(reshape2)
library(tidyverse)
library(vegan)
library(Rfast)
library(car)
library(ggplot2)
library(plyr)

# # Argument inputs
Args <- commandArgs(trailingOnly=TRUE)
# #input variables will need to be a number for the seed and all input files. 
data_table_tabs <- read.table (Args[1], sep = "\t", header = TRUE)
data_table_comma <- read.table(Args[1], sep = ",", header = TRUE)

# # Delete me when done
data_table_comma <- read.table("/home/drt83172/Documents/QPCR_Data_Wrangler/Program/int_files/edit_me.txt", sep = ",", header = TRUE)

# Checking if columns are tabs or comma seperated
if (columns_comma == 5){
  true_data <- data_table_comma
  print("Your file is comma seperated")
} else if (columns_tabs == 5){
  true_data <- columns_tabs
  print("Your file is tab seperated")
} else {
  print("Something is wrong with the file. Are you sure there is 5 columns seperated by tabs or commas")
}

# # Checking amount of columns in files
# # data_table_comma <- read.table("/home/drt83172/Documents/QPCR_Data_Wrangler/Program/int_files/edit_me.txt", sep = ",", header = TRUE)
columns_comma <- ncol(data_table_comma)
columns_tabs <- ncol(data_table_tabs)

# # calculating average of samples
QPCRmeans = ddply(data_table_comma, .(Treatment), summarize, 
      mean = mean(Cp))
QPCRmeans2 <- subset(QPCRmeans, select = -c(Treatment))
rownames(QPCRmeans2) <- QPCRmeans[,1]

# # Getting the standard curve data
std1x = 150
std2x = 15
std3x = 1.5
std4x = .15
std1y = QPCRmeans2["std1","mean"]
std2y = QPCRmeans2["std2","mean"]
std3y = QPCRmeans2["std3","mean"]
std4y = QPCRmeans2["std4","mean"]

# # Creating a standard curve matrix
stdData = matrix(nrow = 4, ncol = 2)
stdData[1,1]=std1x
stdData[2,1]=std2x
stdData[3,1]=std3x
stdData[4,1]=std4x
stdData[1,2]=std1y
stdData[2,2]=std2y
stdData[3,2]=std3y
stdData[4,2]=std4y

# # finding the sum of X*Y
xy = matrix(nrow = 4, ncol = 1)
i=1
for(i in 1:nrow(stdData)){
  xy[i,1] = stdData[i,1] * stdData[i,2]
}
XY <- sum(xy)
# # finding the sum of X
X <- sum(stdData[,1])
# # finding the sum of Y
Y <- sum(stdData[,2])
# # finding X^2
x2 = matrix(nrow = 4, ncol = 1)
for(i in 1:nrow(stdData)){
  x2[i,1] = stdData[i,1]^2
}
X2 <- sum(x2)

# # finding m and b
n = nrow(stdData)
m = ((n * XY) - (X * Y)) / ((n * X2) - X^2)
b = (Y-(m*X))/n

# # using line of best fit equation to get ng of DNA
PredDNA <- QPCRmeans2 
names(PredDNA)[names(PredDNA) == 'mean'] <- 'ng.of.DNA'
samples <- rownames(PredDNA)
for (i in 1:length(samples)){
  DNAng <- (QPCRmeans2[samples[i],"mean"]-b)/m
  PredDNA[samples[i],"ng.of.DNA"] <- DNAng
}
PredDNA$Cp <- QPCRmeans2$mean



################################## Making Graphs ###############################

#Making plot based on sample names
ggplot(data = true_data, 
  aes(
  x = Pos,
  y = Cp,
  color = Treatment,
  shape = Primer_Set
  )) +
  geom_point(size = 5) +
  theme(legend.position = "right") +
  theme_bw()
ggsave(file="Treatments.png")

#Making plot based on endophyte status
ggplot(data = true_data, 
       aes(
         x = Pos,
         y = Cp,
         color = EndoPos_Neg_Water,
         shape = Primer_Set
       )) +
  geom_point(size = 5) +
  theme(legend.position = "right") +
  theme_bw()
ggsave(file="Endophyte_status.png")


# # Graph of line of best fit
ggplot(data = PredDNA, 
       aes(
         x = ng.of.DNA,
         y = Cp,
       )) +
  geom_point(size = 5) +
  theme(legend.position = "right") +
  geom_text(x=50, y=27.5,label=(paste0("slope = ",m))) +
  theme_bw()





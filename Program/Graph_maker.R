#!/usr/bin/env Rscript
library(reshape2)
library(tidyverse)
library(vegan)
library(Rfast)
library(car)
library(ggplot2)
library(plyr)
library(data.table)

# # Argument inputs
Args <- commandArgs(trailingOnly=TRUE)
# #input variables will need to be a number for the seed and all input files. 
data_table_tabs <- read.table (Args[1], sep = "\t", header = TRUE)
data_table_comma <- read.table(Args[1], sep = ",", header = TRUE)
Length <- strtoi(Args[2])

# # # Testing variables # # # # 
# Length <- 118
# data_table_comma <- read.table("/home/drt83172/Documents/QPCR_Data_Wrangler/Program/int_files/edit_me.txt", sep = ",", header = TRUE)
# # # # #

# # Checking amount of columns in files
columns_comma <- ncol(data_table_comma)
columns_tabs <- ncol(data_table_tabs)

# # Checking if columns are tabs or comma seperated
if (columns_comma == 6){
  true_data <- data_table_comma
  print("Your file is comma seperated")
} else if (columns_tabs == 6){
  true_data <- columns_tabs
  print("Your file is tab seperated")
} else {
  print("Something is wrong with the file. Are you sure there is 5 columns seperated by tabs or commas")
}


# # Splitting data  into samples and standard curve data

std_data <- true_data[true_data$Treatment %like% "std", ]
water_data <- true_data[true_data$Treatment == regex('water', ignore_case = TRUE), ]
std_data <- rbind(std_data,water_data)

sample_data <- true_data[!grepl("std", true_data$Treatment),]
sample_data <- sample_data[sample_data$Treatment != regex('water', ignore_case = TRUE),]

# # calculating average of samples and standards
Sample_means = ddply(sample_data, .(Treatment), summarize, 
      meanCP = mean(Cp))
Std_means = ddply(std_data, .(Treatment), summarize, 
                     meanCP = mean(Cp), ngDNA = mean(Concentration)*5) 
Std_means$CopyNumber <- (Std_means$ngDNA * 6.02214076*10^23)/ (Length *650 * 10^9)

################ Making functions for finding line components ##################
Find_XY <- function(data_col) {
  Std_means$xy <- Std_means$meanCP*data_col
  XY <- sum(Std_means$xy)
  return(XY)
}

Find_X2 <- function(data_col){
  Std_means$x2 <- data_col^2
  X2 <- sum(Std_means$x2)
  return(X2)
}

Find_X <- function(data_col){
  X <- sum(data_col)
  return(X)
}

Find_Y <- function(data_col){
  Y <- sum(data_col)
  return(Y)
}


# # Finding standard curve for ng of DNA
ng_XY <- Find_XY(Std_means$ngDNA)
ng_X2 <- Find_X2(Std_means$ngDNA)
ng_X <- Find_X(Std_means$ngDNA)
ng_Y <- Find_Y(Std_means$ngDNA)
ng_n = nrow(Std_means)
ng_m = ((ng_n * ng_XY) - (ng_X * ng_Y)) / ((ng_n * ng_X2) - ng_X^2)
ng_b = (ng_Y-(ng_m*ng_X))/ng_n

# # using line of best fit equation to get ng of DNA
Sample_means$ngDNA <- round(((Sample_means$meanCP)-ng_b)/ng_m, digits = 5)
ng_r <- round(cor(Std_means$ngDNA,Std_means$meanCP), digits = 4)
print(paste0("r squared value for ng of DNA is  ", ng_r^2))

# # Finding standard curve for Copy Number
CN_XY <- Find_XY(Std_means$CopyNumber)
CN_X2 <- Find_X2(Std_means$CopyNumber)
CN_X <- Find_X(Std_means$CopyNumber)
CN_Y <- Find_Y(Std_means$meanCP)
CN_n = nrow(Std_means)
CN_m = ((CN_n * CN_XY) - (CN_X * CN_Y)) / ((CN_n * CN_X2) - CN_X^2)
CN_b = (CN_Y-(CN_m*CN_X))/CN_n

# # using line of best fit equation to get copy number 
Sample_means$CopyNumber <- round(((Sample_means$meanCP)-CN_b)/CN_m, digits = 0)
CN_r <- round(cor(Std_means$CopyNumber,Std_means$meanCP), digits = 4)
print(paste0("r squared value for Copy Number is ", CN_r^2))

####### creating log of starting DNA to calculate efficiency ###################
std_means_only <- Std_means[Std_means$Treatment %like% "std", ]
std_means_only$ngDNAlog <- log10(std_means_only$ngDNA)
  
# # finding the sum of X*Y
std_means_only$Lxy <- std_means_only$meanCP*std_means_only$ngDNAlog
LXY <- sum(std_means_only$Lxy)

# # finding X^2
std_means_only$Lx2 <- std_means_only$ngDNAlog^2
LX2 <- sum(std_means_only$Lx2)

# # finding the sum of X
LX <- sum(std_means_only$ngDNAlog)

# # finding the sum of Y
LY <- sum(std_means_only$meanCP)

# # finding m and b
Ln = nrow(std_means_only)
Lm = ((Ln * LXY) - (LX * LY)) / ((Ln * LX2) - LX^2)
Lb = (LY-(Lm*LX))/Ln
efficiency <- as.data.frame(std_means_only)
E <- (-1+10^(-1/Lm))*100
print(paste0("efficiency is ", E))
print(paste0("Targets for r squared is >98 and for efficiency is 90-110"))
write.csv(Sample_means,"Sample_Data.csv", row.names = TRUE)
################################## Making Graphs ###############################

#Making plot based on sample names
true_data$Pos <- factor(true_data$Pos, levels = true_data$Pos) # # makes the pos column the order on x axis
ggplot(data = true_data, 
  aes(
  x = Pos,
  y = Cp,
  color = Treatment,
  )) +
  geom_point(size = 5) +
  theme(legend.position = "right") +
  theme_bw()

# # Graph of line of best fit
ggplot() +
  geom_point(data = Std_means, aes(x = Std_means$CopyNumber, y = Std_means$meanCP), color='black', size = 5) +
  geom_point(data = Sample_means, aes(x = Sample_means$CopyNumber, y = Sample_means$meanCP), color='Green') +
  annotate(geom="text",x=150000000, y=30, label= paste0("r squared value of standards for copy number is ", CN_r^2))
ggsave(file="Standard_Curve.png")




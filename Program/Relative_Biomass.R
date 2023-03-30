#!/usr/bin/env Rscript

# This code creates ratios of endophyte/grass
library(reshape2)
library(tidyverse)
library(vegan)
library(Rfast)
library(car)
library(ggplot2)
library(data.table)



# # Argument inputs
Args <- commandArgs(trailingOnly=TRUE)
# #input variables will need to be a number for the seed and all input files. 
Endo_data <- read.csv(Args[1], header = TRUE)
Endo_data = subset(Endo_data, select = -c(X))
Endo_data <- Endo_data[!Endo_data$Treatment=="",]
Fescue_data <- read.csv(Args[2], header = TRUE)
Fescue_data = subset(Fescue_data, select = -c(X))
Fescue_data <- Fescue_data[!Fescue_data$Treatment=="",]


# # Practice Data
# Endo_data <- read.csv("/home/drt06/Documents/QPCR_Data_Wrangler/Program/output/Data_for_project/Seperete_Amplicon_Data/Clone_Parents_2x2_Dec_2022_CPV.csv", header = TRUE)
# Endo_data = subset(Endo_data, select = -c(X))
# Endo_data <- Endo_data[!Endo_data$Treatment=="",]
# Fescue_data <- read.csv("/home/drt06/Documents/QPCR_Data_Wrangler/Program/output/Data_for_project/Seperete_Amplicon_Data/Clone_Parents_G3P4_Dec_2022_CPV.csv", header = TRUE)
# Fescue_data = subset(Fescue_data, select = -c(X))
# Fescue_data <- Fescue_data[!Fescue_data$Treatment=="",]

# Creating Ratios and deleting other data
RelBioMass <-  merge(Endo_data, Fescue_data, by.x = "Treatment",  by.y = "Treatment")
RelBioMass$CP_Ratio <- (RelBioMass$meanCP.x/RelBioMass$meanCP.y)
RelBioMass = subset(RelBioMass, select = -c(meanCP.x,meanCP.y))
RelBioMass$CopyNum <- (RelBioMass$LogCopyNumber.x/RelBioMass$LogCopyNumber.y)
RelBioMass = subset(RelBioMass, select = -c(CopyNumber.x,CopyNumber.y))
RelBioMass = subset(RelBioMass, select = -c(LogCopyNumber.x,LogCopyNumber.y,ngDNA.x,ngDNA.y))

# Save the Data frame
write.csv(RelBioMass,"Biomass_Data.csv", row.names = TRUE)

# Exploring the data
Exploring <- RelBioMass
Exploring$CP_CN <- RelBioMass$CP_Ratio - RelBioMass$CopyNum
mean(Exploring$CP_CN)



#!/usr/bin/env Rscript
library(reshape2)
library(tidyverse)
library(vegan)
library(Rfast)
library(car)
library(ggplot2)


# # Argument inputs
Args <- commandArgs(trailingOnly=TRUE)
# #input variables will need to be a number for the seed and all input files. 
data_table_tabs <- read.table (Args[1], sep = "\t", header = TRUE)
data_table_comma <- read.table(Args[1], sep = ",", header = TRUE)

# # Checking amount of columns in files
# data_table_comma <- read.table("/home/drt83172/Documents/QPCR_Data_Wrangler/Program/int_files/edit_me.txt", sep = ",", header = TRUE)
columns_comma <- ncol(data_table_comma)
columns_tabs <- ncol(data_table_tabs)

if (columns_comma == 5){
  true_data <- data_table_comma
  print("Your file is comma seperated")
} else if (columns_tabs == 5){
  true_data <- columns_tabs
  print("Your file is tab seperated")
} else {
  print("Something is wrong with the file. Are you sure there is 5 columns seperated by tabs or commas")
}

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



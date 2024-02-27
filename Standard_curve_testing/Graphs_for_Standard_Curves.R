#!/usr/bin/env Rscript
library(tidyverse)
library(reshape2)
library(tidyverse)
library(vegan)
library(Rfast)
library(RcppGSL)
library(car)
library(ggplot2)
library(plyr)
library(data.table)

# # # Argument inputs
# Args <- commandArgs(trailingOnly=TRUE)
# # #input variables will need to be a number for the seed and all input files. 
# data_table_tabs <- read.table (Args[1], sep = "\t", header = TRUE)
# data_table_comma <- read.table(Args[1], sep = ",", header = TRUE)
# Length <- strtoi(Args[2])

# # # Testing variables # # # # 
data_table_comma <- read.table("/home/drt06/Documents/QPCR_Data_Wrangler/Standard_curve_testing/Input_Data/run1.csv", sep = ",", header = FALSE)
roche480_data <- read.table("/home/drt06/Documents/QPCR_Data_Wrangler/Standard_curve_testing/int_files/edit_me.txt", sep = ",", header = TRUE)
# # # # # #

#Checking to make sure you have the correct amount of columns 
columns_comma <- ncol(roche480_data)
if (columns_comma != 6){
  print("There are more than 5 columns please fix that.") 
}

# Rename a column 
colnames(roche480_data)[colnames(roche480_data) == "Cp"] <- "CT"

# Separating the water data
water_data <- roche480_data[grepl("water", roche480_data$Treatment), ]
roche480_data <- roche480_data %>% filter(!grepl("water", Treatment))

# Separate the sample data
roche480_data_sample <- roche480_data %>% filter(!grepl("std", Treatment))
roche480_data_standards <- roche480_data[grepl("std", roche480_data$Treatment), ]

# getting the log of the concentreation for standards data
roche480_data_standards$logCon <- log10(roche480_data_standards$Concentration)

# summarize data 
standards_summerized <- roche480_data_standards %>%
  group_by(Primer_Set, logCon) %>%
  summarise(avrage = mean(CT))

sample_summerized <- roche480_data_sample %>%
  group_by(Treatment, Concentration) %>%
  summarise(avrage = mean(CT))

# Getting all variables needed from linear model

standard_lm <- lm(roche480_data_standards$CT ~ roche480_data_standards$logCon)

coefficients <- coef(standard_lm)
y_intercept <- coefficients[1]
slope <- coefficients[2]
summary_stats <- summary(standard_lm)
r_squared <- round(summary_stats$r.squared, digits = 3)
efficiency <- 10^(-1/slope) -1

print(paste0("r^2 is  ", r_squared))
print(paste0("slope is ", slope))
print(paste0("y-intercept is  ", y_intercept))
print(paste0("efficiency is  ", efficiency))


# Find the concentration of the sample data using the linear equation
sample_summerized$logCon <- (sample_summerized$avrage - y_intercept) / slope
sample_summerized


# graph of the standards and samples
ggplot() +
  geom_point(data = standards_summerized, aes(x = logCon, y = avrage), color = "blue", size = 3) +
  geom_text(aes(label = paste0("r^2 is  ", r_squared), x = -2, y = 25)) +
  geom_text(aes(label = paste0("efficiency is  ", efficiency), x = -2, y = 24)) +
  geom_point(data = sample_summerized, aes(x = logCon, y = avrage ), color = "red", size = 3) +
  labs(x = "Log concentration Value", y = "CT Values", title = "Epichloe Standards on GGBC Roche 480") +
  theme_bw()

plot(standards_summerized$logCon, standards_summerized$avrage)













###### Below this is data analysis for table top QPCR ########
# # Checking amount of columns in files
columns_comma <- ncol(data_table_comma)
if (columns_comma != 5){
 print("There are more than 5 columns please fix that.") 
}

#Renaming the data set columns.
new_names <- c("Position", "fluorophore", "Primer", "Concentration", "CT")
colnames(data_table_comma)[1:5] <- new_names

#Removing water
water_data <- data_table_comma[grepl("water", data_table_comma$Concentration), ]
data_table_comma <- data_table_comma %>% filter(!grepl("water", Concentration))

#Logging the Concentration
data_table_comma$Concentration <- as.numeric(data_table_comma$Concentration)
data_table_comma$logCon <- log10(data_table_comma$Concentration)

# Making a different column for every primer set
# data_table_comma2 <- spread(data_table_comma, key = Primer, value = logCon, fill = NA)
data_table_comma <- data_table_comma[!is.na(data_table_comma$CT),]
data_table_comma2 <- data_table_comma %>%
  group_by(Primer, Concentration, logCon) %>%
  summarise(avrage = mean(CT))


# Linear model
modelSet <- lm(data_table_comma$CT ~ data_table_comma$logCon)
y_intercept <- modelSet$coefficients[1]
slope <- modelSet$coe
summary(modelSet)

# ggplot
ggplot(data = data_table_comma2, aes(x=logCon, y=avrage, color=Primer)) +
  geom_point(shape = 18, size = 5) +
  labs(x = "logged concentration", y = "CT Values", title = "Small QPCR machine CT Standard Curves")






#!/usr/bin/env Rscript
library(tidyverse)
library(reshape2)
library(vegan)


# # # Argument inputs
# Args <- commandArgs(trailingOnly=TRUE)
# # #input variables will need to be a number for the seed and all input files. 
# data_table_tabs <- read.table (Args[1], sep = "\t", header = TRUE)
# data_table_comma <- read.table(Args[1], sep = ",", header = TRUE)
# Length <- strtoi(Args[2])

# # # Testing variables # # # # 
Standard_curve_test_epichloe <- read.table("/home/drt06/Documents/QPCR_Data_Wrangler/Standard_curve_testing/int_files/Epichloe_Data_Test.txt", sep = ",", header = TRUE)
roche480_data <- read.table("/home/darrian/Desktop/UGA/QPCR_Data_Wrangler/QPCR_Data_Wrangler/Standard_curve_testing/int_files/edit_me.csv", sep = ",", header = TRUE)
# # # # # #

#Checking to make sure you have the correct amount of columns 
columns_comma <- ncol(roche480_data)
if (columns_comma != 6){
  print("There are more or less than 6 columns please fix that.") 
}

# Rename a column 
colnames(roche480_data)[colnames(roche480_data) == "Cp"] <- "CT"

# Separating the water data
water_data <- roche480_data[grepl("water", roche480_data$Treatment), ]
roche480_data <- roche480_data %>% filter(!grepl("water", Treatment))

# Separate the sample data
roche480_data_sample <- roche480_data %>% filter(!grepl("std", Treatment))
roche480_data_standards <- roche480_data[grepl("std", roche480_data$Treatment), ]

# Removing samples with a concentration of 0
roche480_data_standards <- subset(roche480_data_standards, CT != 0)


# getting the log of the concentreation for standards data
roche480_data_standards$logCon <- log10(roche480_data_standards$Concentration)

# small plot for checking before getting avarages
ggplot() +
  geom_point(data = roche480_data_standards, aes(x = logCon, y = CT, color = Primer_Set), size = 3) +
  labs(x = "Log concentration Value", y = "CT Values", title = "Epichloe Standards on GGBC Roche 480") +
  theme_bw()


# summarize data 
standards_summerized <-
  roche480_data_standards %>%
  group_by(Treatment, Primer_Set, logCon) %>%
  summarise(avrage = mean(CT))

# sample_summerized <- roche480_data_sample %>%
#   group_by(Treatment, Concentration) %>%
#   summarise(avrage = mean(CT))

#Splitting data by different Primers
roche480_data_standards_p1 <- subset(roche480_data_standards, Primer_Set == "1x1")
roche480_data_standards_p2 <- subset(roche480_data_standards, Primer_Set == "1x2")
roche480_data_standards_p3 <- subset(roche480_data_standards, Primer_Set == "2x1")
roche480_data_standards_p4 <- subset(roche480_data_standards, Primer_Set == "2x2")


# Getting all variables needed from linear model p1
standard_lm <- lm(roche480_data_standards_p1$CT ~ roche480_data_standards_p1$logCon)

coefficients_p1 <- coef(standard_lm)
y_intercept_p1 <- coefficients_p1[1]
slope_p1 <- coefficients_p1[2]
summary_stats_p1 <- summary(standard_lm)
r_squared_p1 <- round(summary_stats_p1$r.squared, digits = 3)
efficiency_p1 <- 10^(-1/slope_p1) -1

# Getting all variables needed from linear model p2
standard_lm <- lm(roche480_data_standards$CT ~ roche480_data_standards$logCon)

coefficients_p2 <- coef(standard_lm)
y_intercept_p2 <- coefficients_p2[1]
slope_p2 <- coefficients_p2[2]
summary_stats_p2 <- summary(standard_lm)
r_squared_p2 <- round(summary_stats_p2$r.squared, digits = 3)
efficiency_p2 <- 10^(-1/slope_p2) -1

# Getting all variables needed from linear model p3
standard_lm <- lm(roche480_data_standards$CT ~ roche480_data_standards$logCon)

coefficients_p3 <- coef(standard_lm)
y_intercept_p3 <- coefficients_p3[1]
slope_p3 <- coefficients_p3[2]
summary_stats_p3 <- summary(standard_lm)
r_squared_p3 <- round(summary_stats_p3$r.squared, digits = 3)
efficiency_p3 <- 10^(-1/slope_p3) -1

# Getting all variables needed from linear model p4
standard_lm <- lm(roche480_data_standards$CT ~ roche480_data_standards$logCon)

coefficients_p4 <- coef(standard_lm)
y_intercept_p4 <- coefficients_p4[1]
slope_p4 <- coefficients_p4[2]
summary_stats_p4 <- summary(standard_lm)
r_squared_p4 <- round(summary_stats_p4$r.squared, digits = 3)
efficiency_p4 <- 10^(-1/slope_p4) -1

# print(paste0("r^2 is  ", r_squared))
# print(paste0("slope is ", slope))
# print(paste0("y-intercept is  ", y_intercept))
# print(paste0("efficiency is  ", efficiency))


# Find the concentration of the sample data using the linear equation
# sample_summerized$logCon <- (sample_summerized$avrage - y_intercept) / slope
# sample_summerized


# graph of the standards and samples
ggplot() +
  geom_point(data = standards_summerized, aes(x = logCon, y = avrage, color = Primer_Set), size = 3) +
  
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






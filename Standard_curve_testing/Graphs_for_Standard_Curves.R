#!/usr/bin/env Rscript
library(tidyverse)
library(reshape2)
library(vegan)
library(data.table)
library(ggpubr)


# # # Argument inputs
# Args <- commandArgs(trailingOnly=TRUE)
# # #input variables will need to be a number for the seed and all input files. 
# data_table_tabs <- read.table (Args[1], sep = "\t", header = TRUE)
# data_table_comma <- read.table(Args[1], sep = ",", header = TRUE)
# Length <- strtoi(Args[2])

# # # Testing variables # # # # 
# Standard_curve_test_epichloe <- read.table("/home/drt06/Documents/QPCR_Data_Wrangler/Standard_curve_testing/int_files/Epichloe_Data_Test.txt", sep = ",", header = TRUE)
# roche480_data <- read.table("/home/darrian/Desktop/UGA/QPCR_Data_Wrangler/QPCR_Data_Wrangler/Standard_curve_testing/int_files/edit_me.csv", sep = ",", header = TRUE)
Epichloe_Data_Final <- read.table("/home/darrian/Desktop/UGA/QPCR_Data_Wrangler/QPCR_Data_Wrangler/Standard_curve_testing/int_files/Epichloe_03_11_2024.csv", sep = ",", header = TRUE)
Tall_Fescue_Data_Final <- read.table("/home/darrian/Desktop/UGA/QPCR_Data_Wrangler/QPCR_Data_Wrangler/Standard_curve_testing/int_files/TF_03_12_2024.csv", sep = ",", header = TRUE)
Epichloe_Data_Problem_stds <- read.table("/home/darrian/Desktop/UGA/QPCR_Data_Wrangler/QPCR_Data_Wrangler/Standard_curve_testing/int_files/Epi_Lit_Primers.csv", sep = ",", header = TRUE)
Tall_Fescue_Data_Problem_stds <- read.table("/home/darrian/Desktop/UGA/QPCR_Data_Wrangler/QPCR_Data_Wrangler/Standard_curve_testing/int_files/TF_STD_Curves_Problem_Primers.csv", sep = ",", header = TRUE)
# # # # # #

#Checking to make sure you have the correct amount of columns 
columns_comma <- ncol(Epichloe_Data_Final)
if (columns_comma != 6){
  print("There are more or less than 6 columns please fix that.") 
}

# Rename a column 
colnames(Epichloe_Data_Final)[colnames(Epichloe_Data_Final) == "Cp"] <- "CT"

# Separating the water data
water_data <- Epichloe_Data_Final[grepl("water", Epichloe_Data_Final$Treatment), ]
Epichloe_Data_Final <- Epichloe_Data_Final %>% filter(!grepl("water", Treatment))

# Separate the sample data
Epichloe_Data_Final_sample <- Epichloe_Data_Final %>% filter(!grepl("std", Treatment))
Epichloe_Data_Final_standards <- Epichloe_Data_Final[grepl("std", Epichloe_Data_Final$Treatment), ]

# Removing samples with a concentration of 0
Epichloe_Data_Final_standards <- subset(Epichloe_Data_Final_standards, CT != 0)


# getting the log of the concentreation for standards data
Epichloe_Data_Final_standards$logCon <- log10(Epichloe_Data_Final_standards$Concentration)

# small plot for checking before getting avarages
ggplot() +
  geom_point(data = Epichloe_Data_Final_standards, aes(x = logCon, y = CT, color = Primer_Set), size = 3) +
  labs(x = "Log concentration Value", y = "CT Values", title = "Epichloe Standards on GGBC Roche 480") +
  theme_bw()


# summarize data 
standards_summerized_epi <-
  Epichloe_Data_Final_standards %>%
  group_by(Treatment, Primer_Set, logCon) %>%
  summarise(avrage = mean(CT))

#Splitting data by different Primers
Epichloe_Data_Final_standards_1x1 <- subset(Epichloe_Data_Final_standards, Primer_Set == "1x1")
Epichloe_Data_Final_standards_1x2 <- subset(Epichloe_Data_Final_standards, Primer_Set == "1x2")
Epichloe_Data_Final_standards_2x1 <- subset(Epichloe_Data_Final_standards, Primer_Set == "2x1")
Epichloe_Data_Final_standards_2x2 <- subset(Epichloe_Data_Final_standards, Primer_Set == "2x2")

Epichloe_Data_Final_standards_1x1$Primer_Set <- "DMAW1"
Epichloe_Data_Final_standards_1x2$Primer_Set <- "DMAW2"
Epichloe_Data_Final_standards_2x1$Primer_Set <- "DMAW3"
Epichloe_Data_Final_standards_2x2$Primer_Set <- "DMAW4"



# Getting all variables needed from linear model 1x1
standard_lm <- lm(Epichloe_Data_Final_standards_1x1$CT ~ Epichloe_Data_Final_standards_1x1$logCon)

coefficients_1x1 <- coef(standard_lm)
y_intercept_1x1 <- coefficients_1x1[1]
slope_1x1 <- coefficients_1x1[2]
summary_stats_1x1 <- summary(standard_lm)
r_squared_1x1 <- round(summary_stats_1x1$r.squared, digits = 3)
efficiency_1x1 <- round(10^(-1/slope_1x1) -1, digits = 3)

# Getting all variables needed from linear model 1x2
standard_lm <- lm(Epichloe_Data_Final_standards_1x2$CT ~ Epichloe_Data_Final_standards_1x2$logCon)

coefficients_1x2 <- coef(standard_lm)
y_intercept_1x2 <- coefficients_1x2[1]
slope_1x2 <- coefficients_1x2[2]
summary_stats_1x2 <- summary(standard_lm)
r_squared_1x2 <- round(summary_stats_1x2$r.squared, digits = 3)
efficiency_1x2 <- round(10^(-1/slope_1x2) -1, digits = 3)

# Getting all variables needed from linear model 2x1
standard_lm <- lm(Epichloe_Data_Final_standards_2x1$CT ~ Epichloe_Data_Final_standards_2x1$logCon)

coefficients_2x1 <- coef(standard_lm)
y_intercept_2x1 <- coefficients_2x1[1]
slope_2x1 <- coefficients_2x1[2]
summary_stats_2x1 <- summary(standard_lm)
r_squared_2x1 <- round(summary_stats_2x1$r.squared, digits = 3)
efficiency_2x1 <- round(10^(-1/slope_2x1) -1, digits = 3)

# Getting all variables needed from linear model 2x2
standard_lm <- lm(Epichloe_Data_Final_standards_2x2$CT ~ Epichloe_Data_Final_standards_2x2$logCon)

coefficients_2x2 <- coef(standard_lm)
y_intercept_2x2 <- coefficients_2x2[1]
slope_2x2 <- coefficients_2x2[2]
summary_stats_2x2 <- summary(standard_lm)
r_squared_2x2 <- round(summary_stats_2x2$r.squared, digits = 3)
efficiency_2x2 <- round(10^(-1/slope_2x2) -1, digits = 3)

# Making a final data table 
rs <- c(r_squared_1x1,r_squared_1x2,r_squared_2x1,r_squared_2x2)
efficiency <- c(efficiency_1x1,efficiency_1x2,efficiency_2x1,efficiency_2x2)
epichloe_table <- data.table(r_squared = rs, efficiency = efficiency)
write.csv(epichloe_table, file = "/home/darrian/Desktop/UGA/QPCR_Data_Wrangler/QPCR_Data_Wrangler/Standard_curve_testing/output/epichloe_table.csv", row.names = FALSE)



# graph of the standards and samples
ggplot() +
  geom_point(data = standards_summerized_epi, aes(x = logCon, y = avrage, color = Primer_Set), size = 3) +
  geom_text(aes(label = paste0("r^2   ", r_squared_1x1), x = -2, y = 27),size = 4, color = "red") +
  geom_text(aes(label = paste0("efficiency =  ", efficiency_1x1), x = -2, y = 26),size = 4, color = "red") +
  geom_text(aes(label = paste0("r^2   ", r_squared_1x2), x = -2, y = 25),size = 4, color = "green") +
  geom_text(aes(label = paste0("efficiency =  ", efficiency_1x2), x = -2, y = 24),size = 4, color = "green") +
  geom_text(aes(label = paste0("r^2   ", r_squared_2x1), x = -2, y = 23),size = 4, color = "blue") +
  geom_text(aes(label = paste0("efficiency =  ", efficiency_2x1), x = -2, y = 22),size = 4, color = "blue") +
  geom_text(aes(label = paste0("r^2   ", r_squared_2x2), x = -2, y = 21),size = 4, color = "purple") +
  geom_text(aes(label = paste0("efficiency =  ", efficiency_2x2), x = -2, y = 20),size = 4, color = "purple") +
  labs(x = "Log Concentration Value", y = "CT Values", title = "Epichloe Standard Curves") +
  theme_bw()


####################################################################################################
### Tall Fescue potion of code
####################################################################################################

#Checking to make sure you have the correct amount of columns 
columns_comma <- ncol(Tall_Fescue_Data_Final)
if (columns_comma != 6){
  print("There are more or less than 6 columns please fix that.") 
}

# Rename a column 
colnames(Tall_Fescue_Data_Final)[colnames(Tall_Fescue_Data_Final) == "Cp"] <- "CT"

# Separating the water data
water_data <- Tall_Fescue_Data_Final[grepl("water", Tall_Fescue_Data_Final$Treatment), ]
Tall_Fescue_Data_Final <- Tall_Fescue_Data_Final %>% filter(!grepl("water", Treatment))

# Separate the sample data
Tall_Fescue_Data_Final_sample <- Tall_Fescue_Data_Final %>% filter(!grepl("std", Treatment))
Tall_Fescue_Data_Final_standards <- Tall_Fescue_Data_Final[grepl("std", Tall_Fescue_Data_Final$Treatment), ]

# Removing samples with a concentration of 0
Tall_Fescue_Data_Final_standards <- subset(Tall_Fescue_Data_Final_standards, CT != 0)


# getting the log of the concentreation for standards data
Tall_Fescue_Data_Final_standards$logCon <- log10(Tall_Fescue_Data_Final_standards$Concentration)

# small plot for checking before getting avarages
ggplot() +
  geom_point(data = Tall_Fescue_Data_Final_standards, aes(x = logCon, y = CT, color = Primer_Set), size = 3) +
  labs(x = "Log concentration Value", y = "CT Values", title = "Epichloe Standards on GGBC Roche 480") +
  theme_bw()


# summarize data 
standards_summerized_tf <-
  Tall_Fescue_Data_Final_standards %>%
  group_by(Treatment, Primer_Set, logCon) %>%
  summarise(avrage = mean(CT))

#Splitting data by different Primers
Tall_Fescue_Data_Final_standards_g3p4 <- subset(Tall_Fescue_Data_Final_standards, Primer_Set == "g3p4")
Tall_Fescue_Data_Final_standards_g3p4 <- Tall_Fescue_Data_Final_standards_g3p4[-5, ]
Tall_Fescue_Data_Final_standards_g3p5 <- subset(Tall_Fescue_Data_Final_standards, Primer_Set == "g3p5")
Tall_Fescue_Data_Final_standards_tfef_alpha <- subset(Tall_Fescue_Data_Final_standards, Primer_Set == "tfef_alpha")

# Getting all variables needed from linear model g3p4
standard_lm <- lm(Tall_Fescue_Data_Final_standards_g3p4$CT ~ Tall_Fescue_Data_Final_standards_g3p4$logCon)

coefficients_g3p4 <- coef(standard_lm)
y_intercept_g3p4 <- coefficients_g3p4[1]
slope_g3p4 <- coefficients_g3p4[2]
summary_stats_g3p4 <- summary(standard_lm)
r_squared_g3p4 <- round(summary_stats_g3p4$r.squared, digits = 3)
efficiency_g3p4 <- round(10^(-1/slope_g3p4) -1, digits = 3)

# Getting all variables needed from linear model g3p5
standard_lm <- lm(Tall_Fescue_Data_Final_standards_g3p5$CT ~ Tall_Fescue_Data_Final_standards_g3p5$logCon)

coefficients_g3p5 <- coef(standard_lm)
y_intercept_g3p5 <- coefficients_g3p5[1]
slope_g3p5 <- coefficients_g3p5[2]
summary_stats_g3p5 <- summary(standard_lm)
r_squared_g3p5 <- round(summary_stats_g3p5$r.squared, digits = 3)
efficiency_g3p5 <- round(10^(-1/slope_g3p5) -1, digits = 3)

# Getting all variables needed from linear model tfef_alpha
standard_lm <- lm(Tall_Fescue_Data_Final_standards_tfef_alpha$CT ~ Tall_Fescue_Data_Final_standards_tfef_alpha$logCon)

coefficients_tfef_alpha <- coef(standard_lm)
y_intercept_tfef_alpha <- coefficients_tfef_alpha[1]
slope_tfef_alpha <- coefficients_tfef_alpha[2]
summary_stats_tfef_alpha <- summary(standard_lm)
r_squared_tfef_alpha <- round(summary_stats_tfef_alpha$r.squared, digits = 3)
efficiency_tfef_alpha <- round(10^(-1/slope_tfef_alpha) -1, digits = 3)

# Making a final data table 
rs <- c(r_squared_g3p4,r_squared_g3p5,r_squared_tfef_alpha)
efficiency <- c(efficiency_g3p4,efficiency_g3p5,efficiency_tfef_alpha)
tall_fescue_table <- data.table(r_squared = rs, efficiency = efficiency)
write.csv(tall_fescue_table, file = "/home/darrian/Desktop/UGA/QPCR_Data_Wrangler/QPCR_Data_Wrangler/Standard_curve_testing/output/tall_fescue_table.csv", row.names = FALSE)



# graph of the standards and samples
ggplot() +
  geom_point(data = standards_summerized_tf, aes(x = logCon, y = avrage, color = Primer_Set), size = 3) +
  geom_text(aes(label = paste0("r^2   ", r_squared_g3p4), x = -2, y = 27),size = 4, color = "red") +
  geom_text(aes(label = paste0("efficiency =  ", efficiency_g3p4), x = -2, y = 26),size = 4, color = "red") +
  geom_text(aes(label = paste0("r^2   ", r_squared_g3p5), x = -2, y = 25),size = 4, color = "green") +
  geom_text(aes(label = paste0("efficiency =  ", efficiency_g3p5), x = -2, y = 24),size = 4, color = "green") +
  geom_text(aes(label = paste0("r^2   ", r_squared_tfef_alpha), x = -2, y = 23),size = 4, color = "blue") +
  geom_text(aes(label = paste0("efficiency =  ", efficiency_tfef_alpha), x = -2, y = 22),size = 4, color = "blue") +
  labs(x = "Log Concentration Value", y = "CT Values", title = "Tall Fescue Standard Curve Promising Primers") +
  theme_bw()


####################################################################################################
# Epichloe Lit Primers
####################################################################################################


Epichloe_Data_Problem_stds <- read.table("/home/darrian/Desktop/UGA/QPCR_Data_Wrangler/QPCR_Data_Wrangler/Standard_curve_testing/int_files/Epi_Lit_Primers.csv", sep = ",", header = TRUE)

#Checking to make sure you have the correct amount of columns 
columns_comma <- ncol(Epichloe_Data_Problem_stds)
if (columns_comma != 6){
  print("There are more or less than 6 columns please fix that.") 
}

# Rename a column 
colnames(Epichloe_Data_Problem_stds)[colnames(Epichloe_Data_Problem_stds) == "Cp"] <- "CT"

# Separating the water data
water_data <- Epichloe_Data_Problem_stds[grepl("water", Epichloe_Data_Problem_stds$Treatment), ]
Epichloe_Data_Problem_stds <- Epichloe_Data_Problem_stds %>% filter(!grepl("water", Treatment))

# Separate the sample data
Epichloe_Data_Problem_stds_sample <- Epichloe_Data_Problem_stds %>% filter(!grepl("std", Treatment))
Epichloe_Data_Problem_stds_standards <- Epichloe_Data_Problem_stds[grepl("std", Epichloe_Data_Problem_stds$Treatment), ]

# Removing samples with a concentration of 0
Epichloe_Data_Problem_stds_standards <- subset(Epichloe_Data_Problem_stds_standards, CT != 0)


# getting the log of the concentreation for standards data
Epichloe_Data_Problem_stds_standards$logCon <- log10(Epichloe_Data_Problem_stds_standards$Concentration)

# small plot for checking before getting avarages
ggplot() +
  geom_point(data = Epichloe_Data_Problem_stds_standards, aes(x = logCon, y = CT, color = Primer_Set), size = 3) +
  labs(x = "Log concentration Value", y = "CT Values", title = "Epichloe Standards Literature") +
  theme_bw()


# summarize data 
standards_summerized_epi_lit <-
  Epichloe_Data_Problem_stds_standards %>%
  group_by(Treatment, Primer_Set, logCon) %>%
  summarise(avrage = mean(CT))

#Splitting data by different Primers
Epichloe_Data_Problem_stds_standards_EndoEF <- subset(Epichloe_Data_Problem_stds_standards, Primer_Set == "EndoEF")
Epichloe_Data_Problem_stds_standards_TC35X <- subset(Epichloe_Data_Problem_stds_standards, Primer_Set == "TC35X")
Epichloe_Data_Problem_stds_standards_ProA <- subset(Epichloe_Data_Problem_stds_standards, Primer_Set == "ProA")
Epichloe_Data_Problem_stds_standards_2x2_2 <- subset(Epichloe_Data_Problem_stds_standards, Primer_Set == "2x2")

Epichloe_Data_Problem_stds_standards_TC35X$Primer_Set <- "TC351/TC352"
# Getting all variables needed from linear model EndoEF
standard_lm <- lm(Epichloe_Data_Problem_stds_standards_EndoEF$CT ~ Epichloe_Data_Problem_stds_standards_EndoEF$logCon)

coefficients_EndoEF <- coef(standard_lm)
y_intercept_EndoEF <- coefficients_EndoEF[1]
slope_EndoEF <- coefficients_EndoEF[2]
summary_stats_EndoEF <- summary(standard_lm)
r_squared_EndoEF <- round(summary_stats_EndoEF$r.squared, digits = 3)
efficiency_EndoEF <- round(10^(-1/slope_EndoEF) -1, digits = 3)

# Getting all variables needed from linear model TC35X
standard_lm <- lm(Epichloe_Data_Problem_stds_standards_TC35X$CT ~ Epichloe_Data_Problem_stds_standards_TC35X$logCon)

coefficients_TC35X <- coef(standard_lm)
y_intercept_TC35X <- coefficients_TC35X[1]
slope_TC35X <- coefficients_TC35X[2]
summary_stats_TC35X <- summary(standard_lm)
r_squared_TC35X <- round(summary_stats_TC35X$r.squared, digits = 3)
efficiency_TC35X <- round(10^(-1/slope_TC35X) -1, digits = 3)

# Getting all variables needed from linear model ProA
standard_lm <- lm(Epichloe_Data_Problem_stds_standards_ProA$CT ~ Epichloe_Data_Problem_stds_standards_ProA$logCon)

coefficients_ProA <- coef(standard_lm)
y_intercept_ProA <- coefficients_ProA[1]
slope_ProA <- coefficients_ProA[2]
summary_stats_ProA <- summary(standard_lm)
r_squared_ProA <- round(summary_stats_ProA$r.squared, digits = 3)
efficiency_ProA <- round(10^(-1/slope_ProA) -1, digits = 3)

# Getting all variables needed from linear model 2x2_2
standard_lm <- lm(Epichloe_Data_Problem_stds_standards_2x2_2$CT ~ Epichloe_Data_Problem_stds_standards_2x2_2$logCon)

coefficients_2x2_2 <- coef(standard_lm)
y_intercept_2x2_2 <- coefficients_2x2_2[1]
slope_2x2_2 <- coefficients_2x2_2[2]
summary_stats_2x2_2 <- summary(standard_lm)
r_squared_2x2_2 <- round(summary_stats_2x2_2$r.squared, digits = 3)
efficiency_2x2_2 <- round(10^(-1/slope_2x2_2) -1, digits = 3)

# Making a final data table 
rs <- c(r_squared_EndoEF,r_squared_TC35X,r_squared_ProA,r_squared_2x2_2)
efficiency <- c(efficiency_EndoEF,efficiency_TC35X,efficiency_ProA,efficiency_2x2_2)
epichloe_table <- data.table(r_squared = rs, efficiency = efficiency)
write.csv(epichloe_table, file = "/home/darrian/Desktop/UGA/QPCR_Data_Wrangler/QPCR_Data_Wrangler/Standard_curve_testing/output/epichloe_lit_table.csv", row.names = FALSE)



# graph of the standards and samples
ggplot() +
  geom_point(data = standards_summerized_epi_lit, aes(x = logCon, y = avrage, color = Primer_Set), size = 3) +
  geom_text(aes(label = paste0("r^2   ", r_squared_EndoEF), x = -2, y = 27),size = 4, color = "red") +
  geom_text(aes(label = paste0("efficiency =  ", efficiency_EndoEF), x = -2, y = 26),size = 4, color = "red") +
  geom_text(aes(label = paste0("r^2   ", r_squared_TC35X), x = -2, y = 25),size = 4, color = "green") +
  geom_text(aes(label = paste0("efficiency =  ", efficiency_TC35X), x = -2, y = 24),size = 4, color = "green") +
  geom_text(aes(label = paste0("r^2   ", r_squared_ProA), x = -2, y = 23),size = 4, color = "blue") +
  geom_text(aes(label = paste0("efficiency =  ", efficiency_ProA), x = -2, y = 22),size = 4, color = "blue") +
  geom_text(aes(label = paste0("r^2   ", r_squared_2x2_2), x = -2, y = 21),size = 4, color = "purple") +
  geom_text(aes(label = paste0("efficiency =  ", efficiency_2x2_2), x = -2, y = 20),size = 4, color = "purple") +
  labs(x = "Log Concentration Value", y = "CT Values", title = "Epichloe Standard Curves Lit Primers") +
  theme_bw()


################################################################################
# Tall Fescue Problem Primers
################################################################################

Tall_Fescue_Data_Problem_stds <- read.table("/home/darrian/Desktop/UGA/QPCR_Data_Wrangler/QPCR_Data_Wrangler/Standard_curve_testing/int_files/TF_STD_Curves_Problem_Primers.csv", sep = ",", header = TRUE)

#Checking to make sure you have the correct amount of columns 
columns_comma <- ncol(Tall_Fescue_Data_Problem_stds)
if (columns_comma != 6){
  print("There are more or less than 6 columns please fix that.") 
}

# Rename a column 
colnames(Tall_Fescue_Data_Problem_stds)[colnames(Tall_Fescue_Data_Problem_stds) == "Cp"] <- "CT"

# Separating the water data
water_data <- Tall_Fescue_Data_Problem_stds[grepl("water", Tall_Fescue_Data_Problem_stds$Treatment), ]
Tall_Fescue_Data_Problem_stds <- Tall_Fescue_Data_Problem_stds %>% filter(!grepl("water", Treatment))

# Separate the sample data
Tall_Fescue_Data_Problem_stds_sample <- Tall_Fescue_Data_Problem_stds %>% filter(!grepl("std", Treatment))
Tall_Fescue_Data_Problem_stds_standards <- Tall_Fescue_Data_Problem_stds[grepl("std", Tall_Fescue_Data_Problem_stds$Treatment), ]

# Removing samples with a concentration of 0
Tall_Fescue_Data_Problem_stds_standards <- subset(Tall_Fescue_Data_Problem_stds_standards, CT != 0)


# getting the log of the concentreation for standards data
Tall_Fescue_Data_Problem_stds_standards$logCon <- log10(Tall_Fescue_Data_Problem_stds_standards$Concentration)

# small plot for checking before getting avarages
ggplot() +
  geom_point(data = Tall_Fescue_Data_Problem_stds_standards, aes(x = logCon, y = CT, color = Primer_Set), size = 3) +
  labs(x = "Log concentration Value", y = "CT Values", title = "Tall Fescue Problem Standards") +
  theme_bw()


# summarize data 
standards_summerized_tf_problem <-
  Tall_Fescue_Data_Problem_stds_standards %>%
  group_by(Treatment, Primer_Set, logCon) %>%
  summarise(avrage = mean(CT))


#Splitting data by different Primers
Tall_Fescue_Data_Problem_stds_standards_g3p4_2 <- subset(Tall_Fescue_Data_Problem_stds_standards, Primer_Set == "g3p4")
Tall_Fescue_Data_Problem_stds_standards_Tf_Gap <- subset(Tall_Fescue_Data_Problem_stds_standards, Primer_Set == "Tf-Gap")
Tall_Fescue_Data_Problem_stds_standards_Tf_ACS <- subset(Tall_Fescue_Data_Problem_stds_standards, Primer_Set == "TF-ACS")
Tall_Fescue_Data_Problem_stds_standards_g3p6 <- subset(Tall_Fescue_Data_Problem_stds_standards, Primer_Set == "g3p6")


# Getting all variables needed from linear model g3p4_2
standard_lm <- lm(Tall_Fescue_Data_Problem_stds_standards_g3p4_2$CT ~ Tall_Fescue_Data_Problem_stds_standards_g3p4_2$logCon)

coefficients_g3p4_2 <- coef(standard_lm)
y_intercept_g3p4_2 <- coefficients_g3p4_2[1]
slope_g3p4_2 <- coefficients_g3p4_2[2]
summary_stats_g3p4_2 <- summary(standard_lm)
r_squared_g3p4_2 <- round(summary_stats_g3p4_2$r.squared, digits = 3)
efficiency_g3p4_2 <- round(10^(-1/slope_g3p4_2) -1, digits = 3)

# Getting all variables needed from linear model Tf_Gap
standard_lm <- lm(Tall_Fescue_Data_Problem_stds_standards_Tf_Gap$CT ~ Tall_Fescue_Data_Problem_stds_standards_Tf_Gap$logCon)

coefficients_Tf_Gap <- coef(standard_lm)
y_intercept_Tf_Gap <- coefficients_Tf_Gap[1]
slope_Tf_Gap <- coefficients_Tf_Gap[2]
summary_stats_Tf_Gap <- summary(standard_lm)
r_squared_Tf_Gap <- round(summary_stats_Tf_Gap$r.squared, digits = 3)
efficiency_Tf_Gap <- round(10^(-1/slope_Tf_Gap) -1, digits = 3)

# Getting all variables needed from linear model Tf_ACS
standard_lm <- lm(Tall_Fescue_Data_Problem_stds_standards_Tf_ACS$CT ~ Tall_Fescue_Data_Problem_stds_standards_Tf_ACS$logCon)

coefficients_Tf_ACS <- coef(standard_lm)
y_intercept_Tf_ACS <- coefficients_Tf_ACS[1]
slope_Tf_ACS <- coefficients_Tf_ACS[2]
summary_stats_Tf_ACS <- summary(standard_lm)
r_squared_Tf_ACS <- round(summary_stats_Tf_ACS$r.squared, digits = 3)
efficiency_Tf_ACS <- round(10^(-1/slope_Tf_ACS) -1, digits = 3)

# Getting all variables needed from linear model g3p6
standard_lm <- lm(Tall_Fescue_Data_Problem_stds_standards_g3p6$CT ~ Tall_Fescue_Data_Problem_stds_standards_g3p6$logCon)

coefficients_g3p6 <- coef(standard_lm)
y_intercept_g3p6 <- coefficients_g3p6[1]
slope_g3p6 <- coefficients_g3p6[2]
summary_stats_g3p6 <- summary(standard_lm)
r_squared_g3p6 <- round(summary_stats_g3p6$r.squared, digits = 3)
efficiency_g3p6 <- round(10^(-1/slope_g3p6) -1, digits = 3)

# Making a final data table 
rs <- c(r_squared_g3p4_2,r_squared_Tf_Gap,r_squared_Tf_ACS,r_squared_g3p6)
efficiency <- c(efficiency_g3p4_2,efficiency_Tf_Gap,efficiency_Tf_ACS,efficiency_g3p6)
epichloe_table <- data.table(r_squared = rs, efficiency = efficiency)
write.csv(epichloe_table, file = "/home/darrian/Desktop/UGA/QPCR_Data_Wrangler/QPCR_Data_Wrangler/Standard_curve_testing/output/Tall_Fescue_Problem_Std_table.csv", row.names = FALSE)



# graph of the standards and samples
ggplot() +
  geom_point(data = standards_summerized_tf_problem, aes(x = logCon, y = avrage, color = Primer_Set), size = 3) +
  geom_text(aes(label = paste0("r^2   ", r_squared_g3p4_2), x = -2, y = 27),size = 4, color = "red") +
  geom_text(aes(label = paste0("efficiency =  ", efficiency_g3p4_2), x = -2, y = 26),size = 4, color = "red") +
  geom_text(aes(label = paste0("r^2   ", r_squared_Tf_Gap), x = -2, y = 25),size = 4, color = "green") +
  geom_text(aes(label = paste0("efficiency =  ", efficiency_Tf_Gap), x = -2, y = 24),size = 4, color = "green") +
  geom_text(aes(label = paste0("r^2   ", r_squared_Tf_ACS), x = -2, y = 23),size = 4, color = "blue") +
  geom_text(aes(label = paste0("efficiency =  ", efficiency_Tf_ACS), x = -2, y = 22),size = 4, color = "blue") +
  geom_text(aes(label = paste0("r^2   ", r_squared_g3p6), x = -2, y = 21),size = 4, color = "purple") +
  geom_text(aes(label = paste0("efficiency =  ", efficiency_g3p6), x = -2, y = 20),size = 4, color = "purple") +
  labs(x = "Log Concentration Value", y = "CT Values", title = "Tall Fescue Standard Curves Bad Primers") +
  theme_bw()

###############################################################################################################
# Making the final graphic
###############################################################################################################

# I need a way to label each biological replicate
label_replicates <- function(df) {
  for (i in 1:nrow(df)) {
    # Extract the second column value and remove the first character
    column_value <- substr(df[i, "Pos"], 2, nchar(df[i, "Pos"] ))
    # Convert the modified value to numeric and divide by 3
    remainder <- as.numeric(column_value) %% 3
    # Assign the remainder to a new column
    df$Remainder[i] <- remainder
  }
  return(df)
  
}

Epichloe_Data_Final_standards_1x1 <- label_replicates(Epichloe_Data_Final_standards_1x1)
Epichloe_Data_Final_standards_1x2 <- label_replicates(Epichloe_Data_Final_standards_1x2)
Epichloe_Data_Final_standards_2x1 <- label_replicates(Epichloe_Data_Final_standards_2x1)
Epichloe_Data_Final_standards_2x2 <- label_replicates(Epichloe_Data_Final_standards_2x2)

Tall_Fescue_Data_Final_standards_g3p4 <- label_replicates(Tall_Fescue_Data_Final_standards_g3p4)
Tall_Fescue_Data_Final_standards_g3p5 <- label_replicates(Tall_Fescue_Data_Final_standards_g3p5)
Tall_Fescue_Data_Final_standards_tfef_alpha <- label_replicates(Tall_Fescue_Data_Final_standards_tfef_alpha)

Epichloe_Data_Problem_stds_standards_EndoEF <- label_replicates(Epichloe_Data_Problem_stds_standards_EndoEF)
Epichloe_Data_Problem_stds_standards_TC35X <- label_replicates(Epichloe_Data_Problem_stds_standards_TC35X)
Epichloe_Data_Problem_stds_standards_ProA <- label_replicates(Epichloe_Data_Problem_stds_standards_ProA)
Epichloe_Data_Problem_stds_standards_2x2_2 <- label_replicates(Epichloe_Data_Problem_stds_standards_2x2_2)

Tall_Fescue_Data_Problem_stds_standards_g3p4_2 <- label_replicates(Tall_Fescue_Data_Problem_stds_standards_g3p4_2)
Tall_Fescue_Data_Problem_stds_standards_Tf_Gap <- label_replicates(Tall_Fescue_Data_Problem_stds_standards_Tf_Gap)
Tall_Fescue_Data_Problem_stds_standards_Tf_ACS <- label_replicates(Tall_Fescue_Data_Problem_stds_standards_Tf_ACS)
Tall_Fescue_Data_Problem_stds_standards_g3p6 <- label_replicates(Tall_Fescue_Data_Problem_stds_standards_g3p6)


# Data sets needed
Epichloe_Data_Final_standards_1x1
Epichloe_Data_Final_standards_1x2
Epichloe_Data_Final_standards_2x1
Epichloe_Data_Final_standards_2x2
Epichloe_Data_Problem_stds_standards_EndoEF 
Epichloe_Data_Problem_stds_standards_TC35X
Epichloe_Data_Problem_stds_standards_ProA 

Tall_Fescue_Data_Final_standards_g3p4 
Tall_Fescue_Data_Final_standards_g3p5 
Tall_Fescue_Data_Problem_stds_standards_g3p6 
Tall_Fescue_Data_Final_standards_tfef_alpha  
Tall_Fescue_Data_Problem_stds_standards_Tf_Gap 
Tall_Fescue_Data_Problem_stds_standards_Tf_ACS 

endo_data <- rbind(Epichloe_Data_Final_standards_1x1, Epichloe_Data_Final_standards_1x2, Epichloe_Data_Final_standards_2x1, Epichloe_Data_Final_standards_2x2, Epichloe_Data_Problem_stds_standards_EndoEF, Epichloe_Data_Problem_stds_standards_TC35X,Epichloe_Data_Problem_stds_standards_ProA)
fescue_data <- rbind(Tall_Fescue_Data_Final_standards_g3p4, Tall_Fescue_Data_Final_standards_g3p5, Tall_Fescue_Data_Problem_stds_standards_g3p6, Tall_Fescue_Data_Final_standards_tfef_alpha, Tall_Fescue_Data_Problem_stds_standards_Tf_Gap, Tall_Fescue_Data_Problem_stds_standards_Tf_ACS)

color_palette_endo<- c("red", "orange", "yellow")
color_palette_fescue <- c("blue", "green", "purple")


# Getting the r squared labels for graphs
rsquared_endo <- do.call(rbind, by(endo_data, endo_data$Primer_Set, function(x) {
  fit <- lm(CT ~ logCon, data = x)
  data.frame(Primer_Set = unique(x$Primer_Set), r_squared = summary(fit)$r.squared)
}))
rsquared_endo$label <- sapply(rsquared_endo$r_squared, function(x) {
  bquote(italic(R)^2 == .(round(x, 3)))
})
rsquared_tf <- do.call(rbind, by(fescue_data, fescue_data$Primer_Set, function(x) {
  fit <- lm(CT ~ logCon, data = x)
  data.frame(Primer_Set = unique(x$Primer_Set), r_squared = summary(fit)$r.squared)
}))

rsquared_tf$label <- sapply(rsquared_tf$r_squared, function(x) {
  bquote(italic(R)^2 == .(round(x, 3)))
})

# Endophyte graph
endo_plot <- ggplot(data = endo_data, aes(x = logCon, y = CT)) +
  geom_point(data = endo_data, aes(x = logCon, y = CT, color = factor(Remainder)), shape = 15, size = 3) +
  scale_color_manual(values = color_palette_endo) + 
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_bw() + 
  labs(color = "Replicate") +
  geom_text(data = rsquared_endo, aes(label = label), x = Inf, y = -Inf, hjust = 3, vjust = -1, parse = TRUE) +
  facet_wrap(vars(Primer_Set), scales = "free") 

# Tall Fescue Graph
fescue_plot <- ggplot(data = fescue_data, aes(x = logCon, y = CT)) +
  geom_point(data = fescue_data, aes(x = logCon, y = CT, color = factor(Remainder)), shape = 15, size = 3) +
  scale_color_manual(values = color_palette_fescue) + 
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_bw() + 
  labs(color = "Replicate") +
  geom_text(data = rsquared_tf, aes(label = label), x = Inf, y = -Inf, hjust = 3, vjust = -1, parse = TRUE) +
  facet_wrap(vars(Primer_Set), scales = "free") 


ggarrange(fescue_plot, endo_plot, ncol=1, nrow=2, heights = c(.4, .6))



# 
# plot_DMA1 <- ggplot() +
#   geom_point(data = Epichloe_Data_Final_standards_1x1, aes(x = logCon, y = CT, color = factor(Remainder)), shape = 15, size = 3) +
#   scale_color_manual(values = color_palette_endo) + 
#   geom_text(aes(label = paste0("r^2 = ", r_squared_1x1), x = -2, y = 27),size = 4, color = "red") +
#   geom_text(aes(label = paste0("efficiency = ", efficiency_1x1), x = -2, y = 26),size = 4, color = "red") +
#   labs(x = "Log Concentration Value", y = "CT Values", title = "", color = "DMA1") +
#   theme_bw()
# 
# plot_DMA2 <- ggplot() +
#   geom_point(data = Epichloe_Data_Final_standards_1x2, aes(x = logCon, y = CT, color = factor(Remainder)), shape = 15, size = 3) +
#   scale_color_manual(values = color_palette_endo) + 
#   geom_text(aes(label = paste0("r^2 = ", r_squared_1x2), x = -2, y = 27),size = 4, color = "red") +
#   geom_text(aes(label = paste0("efficiency = ", efficiency_1x2), x = -2, y = 26),size = 4, color = "red") +
#   labs(x = "Log Concentration Value", y = "CT Values", title = "", color = "DMA2") +
#   theme_bw()
# 
# plot_DMA3 <- ggplot() +
#   geom_point(data = Epichloe_Data_Final_standards_2x1, aes(x = logCon, y = CT, color = factor(Remainder)), shape = 15, size = 3) +
#   scale_color_manual(values = color_palette_endo) + 
#   geom_text(aes(label = paste0("r^2 = ", r_squared_2x1), x = -2, y = 27),size = 4, color = "red") +
#   geom_text(aes(label = paste0("efficiency = ", efficiency_2x1), x = -2, y = 26),size = 4, color = "red") +
#   labs(x = "Log Concentration Value", y = "CT Values", title = "", color = "DMA3") +
#   theme_bw()
# 
# plot_DMA4 <- ggplot() +
#   geom_point(data = Epichloe_Data_Final_standards_2x2, aes(x = logCon, y = CT, color = factor(Remainder)), shape = 15, size = 3) +
#   scale_color_manual(values = color_palette_endo) + 
#   geom_text(aes(label = paste0("r^2 = ", r_squared_2x2), x = -2, y = 27),size = 4, color = "red") +
#   geom_text(aes(label = paste0("efficiency = ", efficiency_2x2), x = -2, y = 26),size = 4, color = "red") +
#   labs(x = "Log Concentration Value", y = "CT Values", title = "", color = "DMA4") +
#   theme_bw()
# 
# plot_EndoEF <- ggplot() +
#   geom_point(data = Epichloe_Data_Problem_stds_standards_EndoEF, aes(x = logCon, y = CT, color = factor(Remainder)), shape = 15, size = 3) +
#   scale_color_manual(values = color_palette_endo) + 
#   geom_text(aes(label = paste0("r^2 = ", r_squared_EndoEF), x = -2, y = 27),size = 4, color = "red") +
#   geom_text(aes(label = paste0("efficiency = ", efficiency_EndoEF), x = -2, y = 26),size = 4, color = "red") +
#   labs(x = "Log Concentration Value", y = "CT Values", title = "", color = "Endo-EF1") +
#   theme_bw()
# 
# plot_TC35X <- ggplot() +
#   geom_point(data = Epichloe_Data_Problem_stds_standards_TC35X, aes(x = logCon, y = CT, color = factor(Remainder)), shape = 15, size = 3) +
#   scale_color_manual(values = color_palette_endo) + 
#   geom_text(aes(label = paste0("r^2 = ", r_squared_TC35X), x = -2, y = 27),size = 4, color = "red") +
#   geom_text(aes(label = paste0("efficiency = ", efficiency_TC35X), x = -2, y = 26),size = 4, color = "red") +
#   labs(x = "Log Concentration Value", y = "CT Values", title = "", color = "TC351/TC352") +
#   theme_bw()
# 
# plot_ProA <- ggplot() +
#   geom_point(data = Epichloe_Data_Problem_stds_standards_ProA, aes(x = logCon, y = CT, color = factor(Remainder)), shape = 15, size = 3) +
#   scale_color_manual(values = color_palette_endo) + 
#   geom_text(aes(label = paste0("r^2 = ", r_squared_ProA), x = -2, y = 27),size = 4, color = "red") +
#   geom_text(aes(label = paste0("efficiency = ", efficiency_ProA), x = -2, y = 26),size = 4, color = "red") +
#   labs(x = "Log Concentration Value", y = "CT Values", title = "", color = "ProA.5/ProA.1") +
#   theme_bw()
# 
# 
# 
# ggplot(mpg, aes(displ, hwy)) +
#   geom_point() +
#   facet_wrap(vars(class), scales = "free")
# 
# 
# 
# # Tall Fescue Primers
# 
# plot_g3p4 <- ggplot() +
#   geom_point(data = Tall_Fescue_Data_Final_standards_g3p4, aes(x = logCon, y = CT, color = factor(Remainder)), shape = 15, size = 3) +
#   scale_color_manual(values = color_palette_fescue) + 
#   geom_text(aes(label = paste0("r^2 = ", r_squared_g3p4), x = -2, y = 27),size = 4, color = "blue") +
#   geom_text(aes(label = paste0("efficiency = ", efficiency_g3p4), x = -2, y = 26),size = 4, color = "blue") +
#   labs(x = "Log Concentration Value", y = "CT Values", title = "", color = "G3P4") +
#   theme_bw()
# 
# plot_g3p5  <- ggplot() +
#   geom_point(data = Tall_Fescue_Data_Final_standards_g3p5 , aes(x = logCon, y = CT, color = factor(Remainder)), shape = 15, size = 3) +
#   scale_color_manual(values = color_palette_fescue) + 
#   geom_text(aes(label = paste0("r^2 = ", r_squared_g3p5 ), x = -2, y = 27),size = 4, color = "blue") +
#   geom_text(aes(label = paste0("efficiency = ", efficiency_g3p5 ), x = -2, y = 26),size = 4, color = "blue") +
#   labs(x = "Log Concentration Value", y = "CT Values", title = "", color = "G3P5") +
#   theme_bw()
# 
# plot_g3p6  <- ggplot() +
#   geom_point(data = Tall_Fescue_Data_Problem_stds_standards_g3p6 , aes(x = logCon, y = CT, color = factor(Remainder)), shape = 15, size = 3) +
#   scale_color_manual(values = color_palette_fescue) + 
#   geom_text(aes(label = paste0("r^2 = ", r_squared_g3p6 ), x = -2, y = 27),size = 4, color = "blue") +
#   geom_text(aes(label = paste0("efficiency = ", efficiency_g3p6 ), x = -2, y = 26),size = 4, color = "blue") +
#   labs(x = "Log Concentration Value", y = "CT Values", title = "", color = "G3P6") +
#   theme_bw()
# 
# plot_tfef_alpha  <- ggplot() +
#   geom_point(data = Tall_Fescue_Data_Final_standards_tfef_alpha , aes(x = logCon, y = CT, color = factor(Remainder)), shape = 15, size = 3) +
#   scale_color_manual(values = color_palette_fescue) + 
#   geom_text(aes(label = paste0("r^2 = ", r_squared_tfef_alpha ), x = -2, y = 27),size = 4, color = "blue") +
#   geom_text(aes(label = paste0("efficiency = ", efficiency_tfef_alpha ), x = -2, y = 26),size = 4, color = "blue") +
#   labs(x = "Log Concentration Value", y = "CT Values", title = "", color = "Tf-EF1-1α") +
#   theme_bw()
# 
# plot_Tf_ACS  <- ggplot() +
#   geom_point(data = Tall_Fescue_Data_Problem_stds_standards_Tf_ACS , aes(x = logCon, y = CT, color = factor(Remainder)), shape = 15, size = 3) +
#   scale_color_manual(values = color_palette_fescue) + 
#   geom_text(aes(label = paste0("r^2 = ", r_squared_Tf_ACS ), x = -2, y = 27),size = 4, color = "blue") +
#   geom_text(aes(label = paste0("efficiency = ", efficiency_Tf_ACS ), x = -2, y = 26),size = 4, color = "blue") +
#   labs(x = "Log Concentration Value", y = "CT Values", title = "", color = "Tf_ACS") +
#   theme_bw()
# 
# plot_Tf_Gap  <- ggplot() +
#   geom_point(data = Tall_Fescue_Data_Problem_stds_standards_Tf_Gap , aes(x = logCon, y = CT, color = factor(Remainder)), shape = 15, size = 3) +
#   scale_color_manual(values = color_palette_fescue) + 
#   geom_text(aes(label = paste0("r^2 = ", r_squared_Tf_Gap ), x = -2, y = 27),size = 4, color = "blue") +
#   geom_text(aes(label = paste0("efficiency = ", efficiency_Tf_Gap ), x = -2, y = 26),size = 4, color = "blue") +
#   labs(x = "Log Concentration Value", y = "CT Values", title = "", color = "Tf_Gap") +
#   theme_bw() 
# 
# # Making the final large graph
# 
# # (plot_DMA1 + plot_DMA2 + plot_DMA3 + plot_DMA4  + plot_TC35X + plot_ProA + plot_annotation(title = 'My Combined Plot')) / (plot_g3p4 + plot_g3p5 + plot_g3p6 + plot_tfef_alpha + plot_Tf_ACS + plot_Tf_Gap) 

#!/usr/bin/env Rscript
library(tidyverse)
library(reshape2)
library(vegan)
library(data.table)


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
Epichloe_Data_Final_standards_p1 <- subset(Epichloe_Data_Final_standards, Primer_Set == "1x1")
Epichloe_Data_Final_standards_p2 <- subset(Epichloe_Data_Final_standards, Primer_Set == "1x2")
Epichloe_Data_Final_standards_p3 <- subset(Epichloe_Data_Final_standards, Primer_Set == "2x1")
Epichloe_Data_Final_standards_p4 <- subset(Epichloe_Data_Final_standards, Primer_Set == "2x2")


# Getting all variables needed from linear model p1
standard_lm <- lm(Epichloe_Data_Final_standards_p1$CT ~ Epichloe_Data_Final_standards_p1$logCon)

coefficients_p1 <- coef(standard_lm)
y_intercept_p1 <- coefficients_p1[1]
slope_p1 <- coefficients_p1[2]
summary_stats_p1 <- summary(standard_lm)
r_squared_p1 <- round(summary_stats_p1$r.squared, digits = 3)
efficiency_p1 <- round(10^(-1/slope_p1) -1, digits = 3)

# Getting all variables needed from linear model p2
standard_lm <- lm(Epichloe_Data_Final_standards_p2$CT ~ Epichloe_Data_Final_standards_p2$logCon)

coefficients_p2 <- coef(standard_lm)
y_intercept_p2 <- coefficients_p2[1]
slope_p2 <- coefficients_p2[2]
summary_stats_p2 <- summary(standard_lm)
r_squared_p2 <- round(summary_stats_p2$r.squared, digits = 3)
efficiency_p2 <- round(10^(-1/slope_p2) -1, digits = 3)

# Getting all variables needed from linear model p3
standard_lm <- lm(Epichloe_Data_Final_standards_p3$CT ~ Epichloe_Data_Final_standards_p3$logCon)

coefficients_p3 <- coef(standard_lm)
y_intercept_p3 <- coefficients_p3[1]
slope_p3 <- coefficients_p3[2]
summary_stats_p3 <- summary(standard_lm)
r_squared_p3 <- round(summary_stats_p3$r.squared, digits = 3)
efficiency_p3 <- round(10^(-1/slope_p3) -1, digits = 3)

# Getting all variables needed from linear model p4
standard_lm <- lm(Epichloe_Data_Final_standards_p4$CT ~ Epichloe_Data_Final_standards_p4$logCon)

coefficients_p4 <- coef(standard_lm)
y_intercept_p4 <- coefficients_p4[1]
slope_p4 <- coefficients_p4[2]
summary_stats_p4 <- summary(standard_lm)
r_squared_p4 <- round(summary_stats_p4$r.squared, digits = 3)
efficiency_p4 <- round(10^(-1/slope_p4) -1, digits = 3)

# Making a final data table 
rs <- c(r_squared_p1,r_squared_p2,r_squared_p3,r_squared_p4)
efficiency <- c(efficiency_p1,efficiency_p2,efficiency_p3,efficiency_p4)
epichloe_table <- data.table(r_squared = rs, efficiency = efficiency)
write.csv(epichloe_table, file = "/home/darrian/Desktop/UGA/QPCR_Data_Wrangler/QPCR_Data_Wrangler/Standard_curve_testing/output/epichloe_table.csv", row.names = FALSE)



# graph of the standards and samples
ggplot() +
  geom_point(data = standards_summerized_epi, aes(x = logCon, y = avrage, color = Primer_Set), size = 3) +
  geom_text(aes(label = paste0("r^2   ", r_squared_p1), x = -2, y = 27),size = 4, color = "red") +
  geom_text(aes(label = paste0("efficiency =  ", efficiency_p1), x = -2, y = 26),size = 4, color = "red") +
  geom_text(aes(label = paste0("r^2   ", r_squared_p2), x = -2, y = 25),size = 4, color = "green") +
  geom_text(aes(label = paste0("efficiency =  ", efficiency_p2), x = -2, y = 24),size = 4, color = "green") +
  geom_text(aes(label = paste0("r^2   ", r_squared_p3), x = -2, y = 23),size = 4, color = "blue") +
  geom_text(aes(label = paste0("efficiency =  ", efficiency_p3), x = -2, y = 22),size = 4, color = "blue") +
  geom_text(aes(label = paste0("r^2   ", r_squared_p4), x = -2, y = 21),size = 4, color = "purple") +
  geom_text(aes(label = paste0("efficiency =  ", efficiency_p4), x = -2, y = 20),size = 4, color = "purple") +
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
Tall_Fescue_Data_Final_standards_p1 <- subset(Tall_Fescue_Data_Final_standards, Primer_Set == "g3p4")
Tall_Fescue_Data_Final_standards_p2 <- subset(Tall_Fescue_Data_Final_standards, Primer_Set == "g3p5")
Tall_Fescue_Data_Final_standards_p3 <- subset(Tall_Fescue_Data_Final_standards, Primer_Set == "tfef_alpha")

# Getting all variables needed from linear model p1
standard_lm <- lm(Tall_Fescue_Data_Final_standards_p1$CT ~ Tall_Fescue_Data_Final_standards_p1$logCon)

coefficients_p1 <- coef(standard_lm)
y_intercept_p1 <- coefficients_p1[1]
slope_p1 <- coefficients_p1[2]
summary_stats_p1 <- summary(standard_lm)
r_squared_p1 <- round(summary_stats_p1$r.squared, digits = 3)
efficiency_p1 <- round(10^(-1/slope_p1) -1, digits = 3)

# Getting all variables needed from linear model p2
standard_lm <- lm(Tall_Fescue_Data_Final_standards_p2$CT ~ Tall_Fescue_Data_Final_standards_p2$logCon)

coefficients_p2 <- coef(standard_lm)
y_intercept_p2 <- coefficients_p2[1]
slope_p2 <- coefficients_p2[2]
summary_stats_p2 <- summary(standard_lm)
r_squared_p2 <- round(summary_stats_p2$r.squared, digits = 3)
efficiency_p2 <- round(10^(-1/slope_p2) -1, digits = 3)

# Getting all variables needed from linear model p3
standard_lm <- lm(Tall_Fescue_Data_Final_standards_p3$CT ~ Tall_Fescue_Data_Final_standards_p3$logCon)

coefficients_p3 <- coef(standard_lm)
y_intercept_p3 <- coefficients_p3[1]
slope_p3 <- coefficients_p3[2]
summary_stats_p3 <- summary(standard_lm)
r_squared_p3 <- round(summary_stats_p3$r.squared, digits = 3)
efficiency_p3 <- round(10^(-1/slope_p3) -1, digits = 3)

# Making a final data table 
rs <- c(r_squared_p1,r_squared_p2,r_squared_p3)
efficiency <- c(efficiency_p1,efficiency_p2,efficiency_p3)
tall_fescue_table <- data.table(r_squared = rs, efficiency = efficiency)
write.csv(tall_fescue_table, file = "/home/darrian/Desktop/UGA/QPCR_Data_Wrangler/QPCR_Data_Wrangler/Standard_curve_testing/output/tall_fescue_table.csv", row.names = FALSE)



# graph of the standards and samples
ggplot() +
  geom_point(data = standards_summerized_tf, aes(x = logCon, y = avrage, color = Primer_Set), size = 3) +
  geom_text(aes(label = paste0("r^2   ", r_squared_p1), x = -2, y = 27),size = 4, color = "red") +
  geom_text(aes(label = paste0("efficiency =  ", efficiency_p1), x = -2, y = 26),size = 4, color = "red") +
  geom_text(aes(label = paste0("r^2   ", r_squared_p2), x = -2, y = 25),size = 4, color = "green") +
  geom_text(aes(label = paste0("efficiency =  ", efficiency_p2), x = -2, y = 24),size = 4, color = "green") +
  geom_text(aes(label = paste0("r^2   ", r_squared_p3), x = -2, y = 23),size = 4, color = "blue") +
  geom_text(aes(label = paste0("efficiency =  ", efficiency_p3), x = -2, y = 22),size = 4, color = "blue") +
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
standards_summerized_epi <-
  Epichloe_Data_Problem_stds_standards %>%
  group_by(Treatment, Primer_Set, logCon) %>%
  summarise(avrage = mean(CT))

#Splitting data by different Primers
Epichloe_Data_Problem_stds_standards_p1 <- subset(Epichloe_Data_Problem_stds_standards, Primer_Set == "EndoEF")
Epichloe_Data_Problem_stds_standards_p2 <- subset(Epichloe_Data_Problem_stds_standards, Primer_Set == "TC35X")
Epichloe_Data_Problem_stds_standards_p3 <- subset(Epichloe_Data_Problem_stds_standards, Primer_Set == "ProA")
Epichloe_Data_Problem_stds_standards_p4 <- subset(Epichloe_Data_Problem_stds_standards, Primer_Set == "2x2")


# Getting all variables needed from linear model p1
standard_lm <- lm(Epichloe_Data_Problem_stds_standards_p1$CT ~ Epichloe_Data_Problem_stds_standards_p1$logCon)

coefficients_p1 <- coef(standard_lm)
y_intercept_p1 <- coefficients_p1[1]
slope_p1 <- coefficients_p1[2]
summary_stats_p1 <- summary(standard_lm)
r_squared_p1 <- round(summary_stats_p1$r.squared, digits = 3)
efficiency_p1 <- round(10^(-1/slope_p1) -1, digits = 3)

# Getting all variables needed from linear model p2
standard_lm <- lm(Epichloe_Data_Problem_stds_standards_p2$CT ~ Epichloe_Data_Problem_stds_standards_p2$logCon)

coefficients_p2 <- coef(standard_lm)
y_intercept_p2 <- coefficients_p2[1]
slope_p2 <- coefficients_p2[2]
summary_stats_p2 <- summary(standard_lm)
r_squared_p2 <- round(summary_stats_p2$r.squared, digits = 3)
efficiency_p2 <- round(10^(-1/slope_p2) -1, digits = 3)

# Getting all variables needed from linear model p3
standard_lm <- lm(Epichloe_Data_Problem_stds_standards_p3$CT ~ Epichloe_Data_Problem_stds_standards_p3$logCon)

coefficients_p3 <- coef(standard_lm)
y_intercept_p3 <- coefficients_p3[1]
slope_p3 <- coefficients_p3[2]
summary_stats_p3 <- summary(standard_lm)
r_squared_p3 <- round(summary_stats_p3$r.squared, digits = 3)
efficiency_p3 <- round(10^(-1/slope_p3) -1, digits = 3)

# Getting all variables needed from linear model p4
standard_lm <- lm(Epichloe_Data_Problem_stds_standards_p4$CT ~ Epichloe_Data_Problem_stds_standards_p4$logCon)

coefficients_p4 <- coef(standard_lm)
y_intercept_p4 <- coefficients_p4[1]
slope_p4 <- coefficients_p4[2]
summary_stats_p4 <- summary(standard_lm)
r_squared_p4 <- round(summary_stats_p4$r.squared, digits = 3)
efficiency_p4 <- round(10^(-1/slope_p4) -1, digits = 3)

# Making a final data table 
rs <- c(r_squared_p1,r_squared_p2,r_squared_p3,r_squared_p4)
efficiency <- c(efficiency_p1,efficiency_p2,efficiency_p3,efficiency_p4)
epichloe_table <- data.table(r_squared = rs, efficiency = efficiency)
write.csv(epichloe_table, file = "/home/darrian/Desktop/UGA/QPCR_Data_Wrangler/QPCR_Data_Wrangler/Standard_curve_testing/output/epichloe_lit_table.csv", row.names = FALSE)



# graph of the standards and samples
ggplot() +
  geom_point(data = standards_summerized_epi, aes(x = logCon, y = avrage, color = Primer_Set), size = 3) +
  geom_text(aes(label = paste0("r^2   ", r_squared_p1), x = -2, y = 27),size = 4, color = "red") +
  geom_text(aes(label = paste0("efficiency =  ", efficiency_p1), x = -2, y = 26),size = 4, color = "red") +
  geom_text(aes(label = paste0("r^2   ", r_squared_p2), x = -2, y = 25),size = 4, color = "green") +
  geom_text(aes(label = paste0("efficiency =  ", efficiency_p2), x = -2, y = 24),size = 4, color = "green") +
  geom_text(aes(label = paste0("r^2   ", r_squared_p3), x = -2, y = 23),size = 4, color = "blue") +
  geom_text(aes(label = paste0("efficiency =  ", efficiency_p3), x = -2, y = 22),size = 4, color = "blue") +
  geom_text(aes(label = paste0("r^2   ", r_squared_p4), x = -2, y = 21),size = 4, color = "purple") +
  geom_text(aes(label = paste0("efficiency =  ", efficiency_p4), x = -2, y = 20),size = 4, color = "purple") +
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
standards_summerized_epi <-
  Tall_Fescue_Data_Problem_stds_standards %>%
  group_by(Treatment, Primer_Set, logCon) %>%
  summarise(avrage = mean(CT))


#Splitting data by different Primers
Tall_Fescue_Data_Problem_stds_standards_p1 <- subset(Tall_Fescue_Data_Problem_stds_standards, Primer_Set == "g3p4")
Tall_Fescue_Data_Problem_stds_standards_p2 <- subset(Tall_Fescue_Data_Problem_stds_standards, Primer_Set == "Tf-Gap")
Tall_Fescue_Data_Problem_stds_standards_p3 <- subset(Tall_Fescue_Data_Problem_stds_standards, Primer_Set == "TF-ACS")
Tall_Fescue_Data_Problem_stds_standards_p4 <- subset(Tall_Fescue_Data_Problem_stds_standards, Primer_Set == "g3p6")


# Getting all variables needed from linear model p1
standard_lm <- lm(Tall_Fescue_Data_Problem_stds_standards_p1$CT ~ Tall_Fescue_Data_Problem_stds_standards_p1$logCon)

coefficients_p1 <- coef(standard_lm)
y_intercept_p1 <- coefficients_p1[1]
slope_p1 <- coefficients_p1[2]
summary_stats_p1 <- summary(standard_lm)
r_squared_p1 <- round(summary_stats_p1$r.squared, digits = 3)
efficiency_p1 <- round(10^(-1/slope_p1) -1, digits = 3)

# Getting all variables needed from linear model p2
standard_lm <- lm(Tall_Fescue_Data_Problem_stds_standards_p2$CT ~ Tall_Fescue_Data_Problem_stds_standards_p2$logCon)

coefficients_p2 <- coef(standard_lm)
y_intercept_p2 <- coefficients_p2[1]
slope_p2 <- coefficients_p2[2]
summary_stats_p2 <- summary(standard_lm)
r_squared_p2 <- round(summary_stats_p2$r.squared, digits = 3)
efficiency_p2 <- round(10^(-1/slope_p2) -1, digits = 3)

# Getting all variables needed from linear model p3
standard_lm <- lm(Tall_Fescue_Data_Problem_stds_standards_p3$CT ~ Tall_Fescue_Data_Problem_stds_standards_p3$logCon)

coefficients_p3 <- coef(standard_lm)
y_intercept_p3 <- coefficients_p3[1]
slope_p3 <- coefficients_p3[2]
summary_stats_p3 <- summary(standard_lm)
r_squared_p3 <- round(summary_stats_p3$r.squared, digits = 3)
efficiency_p3 <- round(10^(-1/slope_p3) -1, digits = 3)

# Getting all variables needed from linear model p4
standard_lm <- lm(Tall_Fescue_Data_Problem_stds_standards_p4$CT ~ Tall_Fescue_Data_Problem_stds_standards_p4$logCon)

coefficients_p4 <- coef(standard_lm)
y_intercept_p4 <- coefficients_p4[1]
slope_p4 <- coefficients_p4[2]
summary_stats_p4 <- summary(standard_lm)
r_squared_p4 <- round(summary_stats_p4$r.squared, digits = 3)
efficiency_p4 <- round(10^(-1/slope_p4) -1, digits = 3)

# Making a final data table 
rs <- c(r_squared_p1,r_squared_p2,r_squared_p3,r_squared_p4)
efficiency <- c(efficiency_p1,efficiency_p2,efficiency_p3,efficiency_p4)
epichloe_table <- data.table(r_squared = rs, efficiency = efficiency)
write.csv(epichloe_table, file = "/home/darrian/Desktop/UGA/QPCR_Data_Wrangler/QPCR_Data_Wrangler/Standard_curve_testing/output/Tall_Fescue_Problem_Std_table.csv", row.names = FALSE)



# graph of the standards and samples
ggplot() +
  geom_point(data = standards_summerized_epi, aes(x = logCon, y = avrage, color = Primer_Set), size = 3) +
  geom_text(aes(label = paste0("r^2   ", r_squared_p1), x = -2, y = 27),size = 4, color = "red") +
  geom_text(aes(label = paste0("efficiency =  ", efficiency_p1), x = -2, y = 26),size = 4, color = "red") +
  geom_text(aes(label = paste0("r^2   ", r_squared_p2), x = -2, y = 25),size = 4, color = "green") +
  geom_text(aes(label = paste0("efficiency =  ", efficiency_p2), x = -2, y = 24),size = 4, color = "green") +
  geom_text(aes(label = paste0("r^2   ", r_squared_p3), x = -2, y = 23),size = 4, color = "blue") +
  geom_text(aes(label = paste0("efficiency =  ", efficiency_p3), x = -2, y = 22),size = 4, color = "blue") +
  geom_text(aes(label = paste0("r^2   ", r_squared_p4), x = -2, y = 21),size = 4, color = "purple") +
  geom_text(aes(label = paste0("efficiency =  ", efficiency_p4), x = -2, y = 20),size = 4, color = "purple") +
  labs(x = "Log Concentration Value", y = "CT Values", title = "Tall Fescue Standard Curves Bad Primers") +
  theme_bw()




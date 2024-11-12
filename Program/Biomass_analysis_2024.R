# Objective : This code uses the alkaloid data and the biomass data from 2024 to find correlations between them.

# Import libraries
library(tidyverse)
library(ggplot2)
library(ggpubr)

#Importing data sets
# This first set of data is all the data corresponding to the first set of standards
all_2x2_loc <- "/home/darrian/Desktop/UGA/QPCR_Data_Wrangler/QPCR_Data_Wrangler/Program/output/Final_Data_2024/all_2x2_2024.csv"
all_g3p4_loc <- "/home/darrian/Desktop/UGA/QPCR_Data_Wrangler/QPCR_Data_Wrangler/Program/output/Final_Data_2024/all_g3p4_2024.csv"
alkaloid <- "/home/darrian/Desktop/UGA/Tall_Fescue/Alkaloid_Data/Alkaloid_concentrations_2024.csv"
parents_loc <- "/home/darrian/Desktop/UGA/Wallace_Lab/Plant_Info/Match_Parent_Data/usable_predicted_parents_double.csv"

all_2x2 <- read.csv(all_2x2_loc, header = TRUE)
all_g3p4 <- read.csv(all_g3p4_loc, header = TRUE)
alkaloid <- read.csv(alkaloid, header = TRUE)
parents <- read.csv(parents_loc, header = TRUE)


# Fixing data to be agglomerated together
all_2x2 <- subset(all_2x2, select = -X)
all_g3p4 <- subset(all_g3p4, select = -X)


all_2x2 <-
  all_2x2 %>% rename(meanCP_TF = meanCP, ngDNA_TF = ngDNA, CopyNumber_TF = CopyNumber, LogCopyNumber_TF = LogCopyNumber)

all_g3p4 <-
  all_g3p4 %>% rename(meanCP_E = meanCP, ngDNA_E = ngDNA, CopyNumber_E = CopyNumber, LogCopyNumber_E = LogCopyNumber)

raw_all <- merge(all_2x2, all_g3p4, by = c("Treatment", "Plate"))
raw_all <- merge(raw_all, alkaloid, by = c("Treatment"))
raw_all$Plate <- as.character(raw_all$Plate)


raw_all$EoverTF <- (raw_all$ngDNA_E / raw_all$ngDNA_TF)


#############################################################
# linear models
#############################################################

AlkaloidbyDNA <- lm(ppbAlk_con ~ EoverTF ,data = raw_all) 
summary_model <- summary(AlkaloidbyDNA)
rsquared <- summary_model$r.squared

#####################################
# Residual data extraction
#####################################

# We want to use the entire data set and take out everything that effects CT values.
normalized_EoTF <- lm(EoverTF ~ Plate, data = raw_all)
summary(normalized_EoTF) # The filtered data

# removing the residuals from the data by getting the mean then subtracting the residuals
adjustedEoTF <- resid(normalized_EoTF)
adjusted_EoverTF <- mean(raw_all$EoverTF) + adjustedEoTF # Optionally adding back the mean to make numbers meaningful. Could just use residuals

# Removing residuals from data
data1 <- data.frame(raw_all$Treatment)
data1$adjusted_EoverTF <- adjusted_EoverTF
data2 <- data.frame(raw_all$Treatment, raw_all$ppbAlk_con)

names(data1)[names(data1) == "raw_all.Treatment"] <- "Treatment"
names(data2)[names(data2) == "raw_all.Treatment"] <- "Treatment"

data_no_residuals <- merge(data1, data2, by = "Treatment", all = TRUE)

# Renaming Columns
new_names <- c("Treatment", "EoTF_ResidualRemoved", "ppbAlkiloids")
names(data_no_residuals) <- new_names

# Adding in metadata and parents
metadata <- subset(raw_all, select = (c(Plate, Treatment)))
data_no_residuals <- merge(data_no_residuals, metadata, by = c("Treatment"))
head(parents)
parents$Origin <- paste(parents$X0, parents$X1, sep = "X")
parents2 <- parents[, c(1, 2, 3, ncol(parents))]
# Split the string by "X"
parts <- strsplit(parents2$Origin, "X")
# Sort the numeric parts and paste them back together with "X" in between
parents2$Origin <- sapply(parts, function(x) paste(sort(as.numeric(x)), collapse = "X"))
parents2 <- parents2[, c(1, 4)]
colnames(parents2)[colnames(parents2) == "X"] <- "Treatment"
# Adding in parents to data
df <- data.frame(
  Treatment = c("310", "314", "312"),
  Origin = c("Parent", "Parent", "Parent")
)
parents2 <- rbind(parents2,df)
data_no_residuals <- merge(data_no_residuals, parents2, by = c("Treatment"))




# getting r squared 
model_residual <- lm(ppbAlkiloids ~ EoTF_ResidualRemoved ,data = data_no_residuals) 
summary_model_residual <- summary(model_residual)
rsquared_residual <- summary_model_residual$r.squared

#####################################
# Making correlation graphs
#####################################



ggplot() +
  geom_point(data = raw_all, aes(x = ngDNA_E, y = ngDNA_TF, color = Plate), size = 2) +
  labs(x = "Endophyte ng DNA", y = "Tall Fescue ng DNA", title = "Biomass Correlation") +
  theme_bw()

ggplot(data = raw_all, aes(x = EoverTF, y = ppbAlk_con)) +
  geom_point(data = raw_all, aes(x = EoverTF, y = ppbAlk_con, color = Plate), size = 2) +
  labs(x = "Endophyte over Fescue ng DNA", y = "Alkaloid ppb", title = "Alkaloid Concentration vs Endophyte Biomass Estimation") +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  annotate("text", x = Inf, y = -Inf, label = paste("R² =", round(rsquared, 3)), 
           hjust = 1.1, vjust = -0.5, size = 5, color = "black") +
  theme_bw()

# residual plot
ggplot(data = data_no_residuals, aes(x = EoTF_ResidualRemoved, y = ppbAlkiloids)) +
  geom_point(data = data_no_residuals, aes(x = EoTF_ResidualRemoved, y = ppbAlkiloids, color = Plate), size = 2) +
  labs(x = "Endophyte over Fescue ng DNA", y = "Alkaloid ppb", title = "Alkaloid Concentration vs Endophyte Biomass Estimation") +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  annotate("text", x = Inf, y = -Inf, label = paste("R² =", round(rsquared_residual, 3)), 
           hjust = 1.1, vjust = -0.5, size = 5, color = "black") +
  theme_bw()















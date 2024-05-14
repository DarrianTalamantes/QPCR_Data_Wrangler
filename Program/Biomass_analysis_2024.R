# Objective : This code will combine all the data using my original standards find a correlation between them

# Import libraries
library(tidyverse)
library(ggplot2)
library(ggpubr)

#Importing data sets
# This first set of data is all the data corresponding to the first set of standards
all_2x2_loc <- "/home/darrian/Desktop/UGA/QPCR_Data_Wrangler/QPCR_Data_Wrangler/Program/output/Final_Data_2024/all_2x2_2024.csv"
all_g3p4_loc <- "/home/darrian/Desktop/UGA/QPCR_Data_Wrangler/QPCR_Data_Wrangler/Program/output/Final_Data_2024/all_g3p4_2024.csv"

all_2x2 <- read.csv(all_2x2_loc, header = TRUE)
all_g3p4 <- read.csv(all_g3p4_loc, header = TRUE)


# Fixing data to be agglomerated together
all_2x2 <- subset(all_2x2, select = -X)
all_g3p4 <- subset(all_g3p4, select = -X)


all_2x2 <-
  all_2x2 %>% rename(meanCP_TF = meanCP, ngDNA_TF = ngDNA, CopyNumber_TF = CopyNumber, LogCopyNumber_TF = LogCopyNumber)

all_g3p4 <-
  all_g3p4 %>% rename(meanCP_E = meanCP, ngDNA_E = ngDNA, CopyNumber_E = CopyNumber, LogCopyNumber_E = LogCopyNumber)

raw_all <- merge(all_2x2, all_g3p4, by = c("Treatment", "Plate"))
raw_all$Plate<-as.character(raw_all$Plate)

####################################
# Making correlation graphs
#####################################



ggplot() +
  geom_point(data = raw_all, aes(x = ngDNA_E, y = ngDNA_TF, color = Plate), size = 2) +
  labs(x = "Endophyte ng DNA", y = "Tall Fescue ng DNA", title = "Biomass Correlation") +
  theme_bw()


















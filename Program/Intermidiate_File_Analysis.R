# Objective: This script will allow for the analysis of standards from various runs and see if samples from different runs have a large difference in CP value
# Loading packages

library(tidyverse)
library(rapport)


## Loading Data 2x2
data_315x320 <- read.table("/home/drt/Desktop/UGA/QPCR_Data_Wrangler/QPCR_Data_Wrangler/Program/int_files/315x320_2x2.csv", sep = ",", header = TRUE)
data_481_560 <- read.table("/home/drt/Desktop/UGA/QPCR_Data_Wrangler/QPCR_Data_Wrangler/Program/int_files/cp_values_481-560_2x2.csv", sep = ",", header = TRUE)
data_961_1040 <- read.table("/home/drt/Desktop/UGA/QPCR_Data_Wrangler/QPCR_Data_Wrangler/Program/int_files/cp_values_961-1040_2x2_2.csv", sep = ",", header = TRUE)
data_1041_1120_bad_stds <- read.table("/home/drt/Desktop/UGA/QPCR_Data_Wrangler/QPCR_Data_Wrangler/Program/int_files/cp_values_1041-1120_2x2_FAILED.csv", sep = ",", header = TRUE)

# Standard Data 2x2
data_standards3 <- read.table("/home/drt/Desktop/UGA/QPCR_Data_Wrangler/QPCR_Data_Wrangler/Program/int_files/Standard_Creation_2x2_Try3.csv", sep = ",", header = TRUE)
data_standards4 <- read.table("/home/drt/Desktop/UGA/QPCR_Data_Wrangler/QPCR_Data_Wrangler/Program/int_files/Standard_Creation_2x2_Try4.csv", sep = ",", header = TRUE)

# loading data G3p4
data_161_240_G3P4 <- read.table("/home/drt/Desktop/UGA/QPCR_Data_Wrangler/QPCR_Data_Wrangler/Program/int_files/cp_values_161-240_g3p4.csv", sep = ",", header = TRUE)
data_481_560_G3P4 <- read.table("/home/drt/Desktop/UGA/QPCR_Data_Wrangler/QPCR_Data_Wrangler/Program/int_files/cp_values_481-560_g3p4.csv", sep = ",", header = TRUE)
data_961_1040_G3P4 <- read.table("/home/drt/Desktop/UGA/QPCR_Data_Wrangler/QPCR_Data_Wrangler/Program/int_files/cp_values_961-1040_g3p4.csv", sep = ",", header = TRUE)

# standard data g3p4
data_g3p4_standards_test <- read.table("/home/drt/Desktop/UGA/QPCR_Data_Wrangler/QPCR_Data_Wrangler/Program/int_files/G3P4_Standards_Red_Blue_Purple_Yellow_allStriped_Black_072623.csv", sep = ",", header = TRUE)

# Getting rid of empty rows

delete_NA <- function(data){
  
  for (x in nrow(data):1){
    if (is.na(data[x,3]) | (data[x,3] == '' | data[x,3]=="water")){
      data <- data[-x,]
      }
  }
  return(data)
}


data_315x320 <- delete_NA(data_315x320)
data_481_560 <- delete_NA(data_481_560)
data_961_1040 <-  delete_NA(data_961_1040)
data_1041_1120_bad_stds <- delete_NA(data_1041_1120_bad_stds)

data_standards3 <- delete_NA(data_standards3)
data_standards4 <- delete_NA(data_standards4)


data_161_240_G3P4 <- delete_NA(data_161_240_G3P4)
data_481_560_G3P4 <- delete_NA(data_481_560_G3P4)
data_961_1040_G3P4 <- delete_NA(data_961_1040_G3P4)

data_g3p4_standards_test <- delete_NA(data_g3p4_standards_test)

### Extracting standards from other data
data_standards1 <- slice(data_315x320, 1:30,  223:225)
data_standards2 <- data_481_560[1:30,]
data_standards11 <- data_961_1040[1:30,]
data_standards22 <- data_1041_1120_bad_stds[1:30,]

### Extracting samples from other data
data_samples1 <- slice(data_315x320, 31:225, 229:240)
data_samples2 <- data_481_560[31:270,]
data_samples11 <- data_961_1040[31:270,]
data_samples22 <- data_1041_1120_bad_stds[31:270,]


data_samples1$Treatment <- "Sample1"
data_samples2$Treatment <- "Sample2"
data_samples11$Treatment <- "Sample3"
data_samples22$Treatment <- "Sample4"

#Extracting standards from g3p4 data
G3P4_standards1 <- data_161_240_G3P4[1:15,]
G3P4_standards2 <- data_481_560_G3P4[1:30,]
G3P4_standards3 <- data_961_1040_G3P4[1:30,]

# Extracting samples from g3p4 data
G3P4_samples1 <- data_161_240_G3P4[16:255,]
G3P4_samples2 <- data_481_560_G3P4[31:270,]
G3P4_samples3 <- data_961_1040_G3P4[31:270,]

G3P4_samples1$Treatment <- "Sample1"
G3P4_samples2$Treatment <- "Sample2"
G3P4_samples3$Treatment <- "Sample3"

# Sepereating Standard Testing Data, to be used seperated as is
data_g3p4_standards_testRB <- data_g3p4_standards_test[1:30,]
data_g3p4_standards_testPY <- data_g3p4_standards_test[31:60,]
data_g3p4_standards_testBlack <- data_g3p4_standards_test[61:81,]

### Combining data again

data2_315x320 <- rbind(data_standards1,data_samples1)
data2_481_560 <- rbind(data_standards2,data_samples2)
data2_961_1040 <- rbind(data_standards11,data_samples11)
data2_1041_1120_bad_stds <- rbind(data_standards22,data_samples22)

# Combining G3P4 data 
data2_161_240_G3P4 <- rbind(G3P4_standards1, G3P4_samples1)
data2_481_560_G3P4 <- rbind(G3P4_standards2, G3P4_samples2)
data2_961_1040_G3P4 <- rbind(G3P4_standards2, G3P4_samples3)

### Combine all data into one big data set
data2_315x320$Data_Set <- "data2_315x320"
data2_481_560$Data_Set <- "data2_481_560"
data2_961_1040$Data_Set <- "data2_961_1040"
data2_1041_1120_bad_stds$Data_Set <- "data2_1041_1120_bad_stds"
data_standards3$Data_Set <- "data_standards3"
data_standards4$Data_Set <- "data_standards4"

alldata<-rbind(data2_315x320,data2_481_560,data2_961_1040,data2_1041_1120_bad_stds,data_standards3,data_standards4)

# Combining the G3P4 data

data2_161_240_G3P4$Data_Set <-"data2_161_240"
data2_481_560_G3P4$Data_Set <- "data2_481_560"
data2_961_1040_G3P4$Data_Set <- "data2_961_1040"
data_g3p4_standards_testPY$Data_Set <- "HG"
data_g3p4_standards_testRB$Data_Set <- "DRT"
data_g3p4_standards_testBlack$Data_Set <- "Black"

alldata_G3P4 <- rbind(data2_161_240_G3P4,data2_481_560_G3P4,data2_961_1040_G3P4,data_g3p4_standards_testPY,data_g3p4_standards_testRB,data_g3p4_standards_testBlack)


### Plotting the data

ggplot(alldata, aes(x=Treatment, y=Cp, color=Data_Set)) +
  geom_point(alpha=.2, size = 2)

ggplot(alldata_G3P4, aes(x=Treatment, y=Cp, color=Data_Set)) +
  geom_point(alpha=.3, size = 2)


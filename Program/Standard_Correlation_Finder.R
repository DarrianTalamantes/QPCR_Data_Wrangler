# Objective : This code will combine all the data using my original standards find a correlation between them

# Import libraries
library(tidyverse)
library(gvlma)
library(mbQTL)
library(ggpubr)

#Importing data sets
# This first set of data is all the data corresponding to the first set of standards
all_2x2_loc <- "/home/drt06/Documents/QPCR_Data_Wrangler/Program/int_files/Data_for_Project/all_2x2_Red_Blue.csv"
all_g3p4_loc <- "/home/drt06/Documents/QPCR_Data_Wrangler/Program/int_files/Data_for_Project/all_g3p4_correspond_to_Red_Blue_Standards.csv"

all_2x2 <- read.csv(all_2x2_loc, header = TRUE)
all_g3p4 <- read.csv(all_g3p4_loc, header = TRUE)

# This set of data corresponds to the second set of standards and a few redoes. Will be added towrds end
all_2x2_2_loc <- "/home/drt06/Documents/QPCR_Data_Wrangler/Program/int_files/Data_for_Project/1041-end-redoes_2x2.csv"
all_g3p4_2_loc <- "/home/drt06/Documents/QPCR_Data_Wrangler/Program/int_files/Data_for_Project/1041-end-redoes_g3p4.csv"  

# This set of data is for parents and cross 305x320, will be added twords end
all_315x320_2x2_loc <- "/home/drt06/Documents/QPCR_Data_Wrangler/Program/int_files/Preliminary_Project_Data/315x320_2x2.txt"
all_315x320_g3p4_loc <- "/home/drt06/Documents/QPCR_Data_Wrangler/Program/int_files/Preliminary_Project_Data/315x320_G3P4.txt"



# Commands I used to make the data set and then I added labels to the data for each set it came from 
# cat cp_values_1-80_g3p4.csv cp_values_81-160_g3p4.csv cp_values_161-240_g3p4.csv cp_values_241-320_g3p4.csv cp_values_321-400_g3p4.csv cp_values_401-480_g3p4.csv cp_values_481-560_g3p4.csv cp_values_561-640_g3p4.csv cp_values_641-720_g3p4.csv cp_values_721-800_g3p4.csv cp_values_801-880_g3p4.csv cp_values_881-960_g3p4.csv cp_values_961-1040_g3p4.csv > all_g3p4_correspond_to_Red_Blue_Standards.csv 
# cat cp_values_1-80_2x2.csv cp_values_81-160_2x2.csv cp_values_161-240_2x2.csv cp_values_241-320_2x2.csv cp_values_321-400_2x2.csv cp_values_401-480_2x2.csv cp_values_481-560_2x2.csv cp_values_561-640_2x2.csv cp_values_641-720_2x2.csv cp_values_721-800_2x2.csv cp_values_801-880_2x2.csv cp_values_881-960_2x2.csv cp_values_961-1040_2x2_2.csv  > all_2x2_Red_Blue.csv

# Cleaning up Standard data 
# Remove any 0's or anything with more than 4 char in Treatment column from the stds data

for (x in nrow(all_2x2):1){
  if (all_2x2[x,2]=="" | all_2x2[x,2]=="0" | nchar(all_2x2[x,3])!=4){
    all_2x2 <- all_2x2[-x,]
  }
}
for (x in nrow(all_g3p4):1){
  if (all_g3p4[x,2]=="" | all_g3p4[x,2]=="0" | nchar(all_g3p4[x,3])!=4){
    all_g3p4 <- all_g3p4[-x,]
  }
}

all_2x2$Cp <- as.numeric(all_2x2$Cp)
all_g3p4$Cp <- as.numeric(all_g3p4$Cp)

# Count 
all_2x2 %>% count(Treatment, sort = TRUE)


# ddply to make the std# and the data set the parameters
stds_g3p4 <-
  all_g3p4 %>%
  group_by(Data_Set,Treatment) %>%
  summarise(mean = mean(Cp))

stds_2x2 <-
  all_2x2 %>%
  group_by(Data_Set,Treatment) %>%
summarise(mean = mean(Cp))

# Graphing out all currently used standard sets to learn which ones to re-evaluate
stds_2x2 <- stds_2x2 %>% 
  rename("Epichloe_Standard_CP_Values" = "mean")
stds_g3p4 <- stds_g3p4 %>% 
  rename("Tall_Fescue_Standard_CP_Values" = "mean")

Graph_Data <- merge(stds_2x2, stds_g3p4, by.x = c("Data_Set", "Treatment"), by.y = c("Data_Set", "Treatment"))

ggplot(Graph_Data, aes(x=Tall_Fescue_Standard_CP_Values, y=Epichloe_Standard_CP_Values, color =Data_Set, labels = Data_Set)) +
  geom_point( size = 2) +
  geom_text(label = Graph_Data$Data_Set, 
            nudge_x = 0.25, nudge_y = 0.25, 
            check_overlap = T) +
  ggtitle(" 2x2 Standards")

ggplot(Graph_Data, aes(x=Tall_Fescue_Standard_CP_Values, y=Epichloe_Standard_CP_Values, color = Treatment, labels = Data_Set)) +
  geom_point( size = 2) +
  geom_text(label = Graph_Data$Data_Set, 
            nudge_x = 0.25, nudge_y = 0.25, 
            check_overlap = T) +
  ggtitle(" 2x2 Standards")

# Remove samples 1-80 for now because something is wrong with the 2x2 part of it.
# Removing samples 1-80 made the model not work 
stds_2x2 <- stds_2x2[!stds_2x2$Data_Set == "1-80",]
stds_g3p4 <- stds_g3p4[!stds_g3p4$Data_Set == "1-80",]

# Retrying the models but with a quick log transformation 
stds_2x2$LogEpichloe_Standard_CP_Values <- log(stds_2x2$Epichloe_Standard_CP_Values)
stds_g3p4$LogTall_Fescue_Standard_CP_Values <- log(stds_g3p4$Tall_Fescue_Standard_CP_Values)

# removing standards from 2x2 to make equal to g3p4 standards
stds_2x2$combo <- paste(stds_2x2$Data_Set, stds_2x2$Treatment, sep="_")
stds_g3p4$combo <- paste(stds_g3p4$Data_Set, stds_g3p4$Treatment, sep="_")
stdsall <- merge(stds_2x2, stds_g3p4, by.x = "combo", by.y = "combo")
stds_2x2 <- subset(stdsall, select = c("Data_Set.x", "Treatment.x", "Epichloe_Standard_CP_Values", "LogEpichloe_Standard_CP_Values"))
# Fixing names 
stds_2x2 <- 
  stds_2x2 %>% 
  rename(Data_Set = Data_Set.x, Treatment = Treatment.x)

stds_g3p4 <- subset(stdsall, select = c("Data_Set.y", "Treatment.y", "Tall_Fescue_Standard_CP_Values", "LogTall_Fescue_Standard_CP_Values"))
stds_g3p4 <- 
  stds_g3p4 %>% 
  rename(Data_Set = Data_Set.y, Treatment = Treatment.y)

# Linear model, We must remove 1-80 to get this to work
model_logged <- lm(stds_2x2$LogEpichloe_Standard_CP_Values ~ stds_g3p4$LogTall_Fescue_Standard_CP_Values)
summary(model_logged)
gvlma(model_logged)

y_intercept <- model_logged$coefficients[1] # the y intercept
TF_Coefficient <- model_logged$coefficients[2] # The tall Fescue coefficient

Graph_Data_logged <- merge(stds_2x2, stds_g3p4, by.x = c("Data_Set", "Treatment"), by.y = c("Data_Set", "Treatment"))

ggplot(Graph_Data_logged, aes(x=LogTall_Fescue_Standard_CP_Values, y=LogEpichloe_Standard_CP_Values, color = Treatment, labels = Data_Set)) +
  geom_point( size = 2) +
  geom_text(label = Graph_Data_logged$Data_Set,
            check_overlap = T)


############################ This function uses that formula to predict 2x2 based on g3p4 numbers  ########################
# The formula is y-intercept + coefficient*value
# This function is based off of logged values so any 2x2 standard i put in here must be logged.
EpiPredictor <- 
  function(x){
    y <- y_intercept + TF_Coefficient*x
    CP <- 2.71828^y  
    return(CP)
      }

EpiPredictor(2.74)


################# Calculating efficiency for each data set and creating adjusted CP values #####################

# Recreating all standard data for use in the rest of the script 
Standards_2x2 <- all_2x2
Standards_g3p4 <- all_g3p4
for (x in nrow(Standards_2x2):1){
  if (Standards_2x2[x,2]=="" | Standards_2x2[x,2]=="0" | nchar(Standards_2x2[x,3])!=4){
    Standards_2x2 <- Standards_2x2[-x,]
  }
}
for (x in nrow(Standards_g3p4):1){
  if (Standards_g3p4[x,2]=="" | Standards_g3p4[x,2]=="0" | nchar(Standards_g3p4[x,3])!=4){
    Standards_g3p4 <- Standards_g3p4[-x,]
  }
}
Standards_2x2$Concentration <- as.numeric(Standards_2x2$Concentration)
Standards_g3p4$Concentration <- as.numeric(Standards_g3p4$Concentration)


# seperating the standards and getting the means of them, calculating ng of DNA
# function takes full data set and seperates standards to find their means
# The function also calculates the copy number and log copy number
FindStandardMeans <- function(stdData, length){
  stdData$Cp <- as.numeric(stdData$Cp)
  stdData$Concentration <- as.numeric(stdData$Concentration)
  for (x in nrow(stdData):1){
    if (stdData[x,2]=="" | stdData[x,2]=="0" | nchar(stdData[x,3])!=4){
      stdData <- stdData[-x,]}
  }
  Standard_Means <-
    stdData %>%
    group_by(Data_Set,Treatment) %>%
    summarise(MeanCP = mean(Cp), NgDNA = mean(Concentration)*5) # we add 5 ul of sample
  Standard_Means$CopyNumber <- (Standard_Means$NgDNA * 6.02214076*10^23/ (length * 650 * 10^9))
  # NgDNA*1 mole / length * 1 mole of base pairs in ng
  Standard_Means$LogCopyNumber <- log10(Standard_Means$CopyNumber)
  return(Standard_Means)
  }

Standard_means_2x2 <- FindStandardMeans(Standards_2x2,118)
Standard_Means_g3p4 <- FindStandardMeans(Standards_g3p4,230)

# # Test with small subset to see if I get different corelation coefficients with cor() and lm()
# Epi1_80 <- subset(Standard_means_2x2, Data_Set == "1-80")
# modelEpi1_80 <- lm(Epi1_80$MeanCP2x2 ~ Epi1_80$LogCopyNumber) # yeild different correlation results
# summary(modelEpi1_80)
# round(cor(Epi1_80$MeanCP2x2,Epi1_80$LogCopyNumber), digits = 4) # yeild different correlation results

# Making model to use to predict copy num from CP for Epichloe
Standard_means_2x2 <- subset(Standard_means_2x2, Data_Set !="1-80")

ggplot(Standard_means_2x2, aes(x=LogCopyNumber, y=MeanCP, color = Treatment, labels = Data_Set)) +
  geom_point( size = 2) +
  geom_text(label = Standard_means_2x2$Data_Set,
            check_overlap = T)

modelEpiStandardsVCopynum <- lm(Standard_means_2x2$LogCopyNumber ~ Standard_means_2x2$MeanCP) # yeild different correlation results
summary(modelEpiStandardsVCopynum)
y_intercept <- modelEpiStandardsVCopynum$coefficients[1] # the y intercept
CP_Coefficient <- modelEpiStandardsVCopynum$coefficients[2] # The tall Fescue coefficient
gvlma(modelEpiStandardsVCopynum)


Standard_means_2x2 <- FindStandardMeans(Standards_2x2,118)

# This function allows you to put in a CP value and it will output an adjusted CP value that takes efficiency into consideration
# Step 1 Calculate line for a single set of standards
# Step 2 Calculate efficiency for standards
# Step 2.5 Calculate the outliars and remove them
# Step 3 Use efficiency to find adjusted CP value
CpAdjuster <- 
  function(model, Data){
    Data$adjCP <- 0
    slope <- model$coefficients[2] # This is not the slope I should use, need to use log copy number slope
    E = -1+10^(-1/slope)
    
    for (i in 1:nrow(Data)){
      ogCP <- Data[i,3]
      ratio =  (2^ogCP) / (2*E)^ogCP # Finding ratio of expected vs actual
      adjuster = log(ratio, base = 2) # Finding amount needed to adjust CP videa
      adjCP = round(ogCP - adjuster, 3) #finalizes CP adjustment
      Data$adjCP[i] <- adjCP
    }
    Data$Efficiency <- E
    return(Data)
  }

# This function takes the data set I combined using cat, all_2x2 is an example
# It will filter out NA's, 0's, standards
# It spits out a table of means for all samples that passed filtering. 
Data_Filtering <- function(Data){
  Data$Cp <- as.numeric(Data$Cp)
  Data[Data == ''] <- NA
# filtering data 2x2
for (x in nrow(Data):1){
    if (is.na(Data[x,3])){
      Data <- Data[-x,]
    }
    else if (is.na(Data[x,2])){
      Data <- Data[-x,]
    }
    else if (Data[x,2]=="0"){
      Data <- Data[-x,]
    }
    else if (Data[x,1]=="Pos"){
      Data <- Data[-x,]
    }
    else if (str_detect(Data[x,3], "std") | Data[x,3] == 0 | Data[x,3] == "water"){
      Data <- Data[-x,]
    }
  }
  # Getting rid of extrenous 0's 2x2
  # TreatmentSplit <- as.data.frame(str_split_fixed(Data$Treatment, "-", 3))
  # TreatmentSplit$V1 <- sub("^0+", "", TreatmentSplit$V1)
  # TreatmentSplit$V2 <- sub("^0+", "", TreatmentSplit$V2)
  # TreatmentSplit$V3 <- sub("^0+", "", TreatmentSplit$V3)
  # TreatmentSplit$x <- paste0(TreatmentSplit$V1, "-", TreatmentSplit$V2, "-", TreatmentSplit$V3)
  # Data$Treatment <- TreatmentSplit$x 
  
  Data_Mean <-
    Data %>%
    group_by(Data_Set,Treatment) %>%
    summarise(MeanCP = mean(Cp))
  return(Data_Mean)
}

# Using the filtering function on both data sets
Data_2x2_Means <- Data_Filtering(all_2x2)
Data_g3p4_Means <- Data_Filtering(all_g3p4)

# Looking for differences in data sets
for (i in 1:nrow(Data_g3p4_Means)){
  if(Data_g3p4_Means[i,2] %in% Data_2x2_Means$Treatment){
  } else 
    print(paste0(Data_g3p4_Means[i,2], " is missing in 2x2 data"))
}
for (i in 1:nrow(Data_2x2_Means)){
  if(Data_2x2_Means[i,2] %in% Data_g3p4_Means$Treatment){
  } else 
    print(paste0(Data_2x2_Means[i,2], " is missing in g3p4 data"))
}

# Using Standards and Data to adjust CP Values 
Data_Sets <- unique(Data_2x2_Means$Data_Set)

# This next function works with 2 data sets. Mean standards and mean data. 
# It will use the standards to calculate efficiency and use that to spit out the data with adjusted CP values
# This function uses another function called "CpAdjuster"
# THis function also removes outliars who are out of standard range or 3x std away from mean
CpAdjusterP2 <- function(Mean_Standards, Mean_Data, Data_Sets){
  Data_Sets <- as.matrix(Data_Sets)
  All_Data = list()
  for (i in 1:nrow(Data_Sets)){
    
    Set <- Data_Sets[i]
    Current_Standard <- subset(Mean_Standards, Data_Set == Set)
    Current_Data <- subset(Mean_Data, Data_Set == Set)
    
    #Removing outliars
    meanCP <- mean(Current_Data$MeanCP)
    SdCP <- sd(Current_Data$MeanCP)
    cutoff1 <- meanCP + 3*SdCP
    cutoff2 <- meanCP - 3*SdCP
    cutoffmin <- min(Current_Standard$MeanCP)
    cutoffmax <- max(Current_Standard$MeanCP)
    Before <- nrow(Current_Data)
    Current_Data <- subset(Current_Data, MeanCP < cutoff1 & MeanCP > cutoff2 & MeanCP > cutoffmin & MeanCP < cutoffmax)
    after <- nrow(Current_Data)
    
    print(paste0(Set, " Dropped ", Before-after))
    
    modelSet <- lm(Current_Standard$MeanCP ~ Current_Standard$LogCopyNumber)
    y_intercept <- modelSet$coefficients[1]
    slope <- modelSet$coefficients[2]
    summary(modelSet)
    Current_Data$LogCopyNumber <- round(((Current_Data$MeanCP)-y_intercept)/slope, digits = 3)
    Current_Data$CopyNumber <- 10^Current_Data$LogCopyNumber # Reverse the log
    
    All_DataC <- CpAdjuster(modelSet, Current_Data) # Need this to remove outliars too
    All_Data[[i]] <-  All_DataC
  }
  All_Data <- bind_rows(All_Data)
  return(All_Data)
  }
Data_2x2_AdjCP <- CpAdjusterP2(Standard_means_2x2, Data_2x2_Means, Data_Sets)
Data_g3p4_AdjCP <- CpAdjusterP2(Standard_Means_g3p4, Data_g3p4_Means, Data_Sets)

############# End function here
# Graphically showing how adjusted CP makes this much more variable t 

ggplot(Data_2x2_AdjCP, aes(x = MeanCP, y = as.numeric(adjCP))) +
  geom_point(size = 2)

EpiAdjCP <- ggplot(Data_2x2_AdjCP, aes(x = Data_Set, y = as.numeric(adjCP))) +
  geom_point(size = 2) +
  ggtitle("Adjusted CP Values Epichloe")

EpiOgCP <- ggplot(Data_2x2_AdjCP, aes(x = Data_Set, y = as.numeric(MeanCP))) +
  geom_point(size = 2) +
  ggtitle("Origional CP Values Epichloe")

FescueAdjCP <- ggplot(Data_g3p4_AdjCP, aes(x = Data_Set, y = as.numeric(adjCP))) +
  geom_point(size = 2) +
  ggtitle("Adjusted CP Values Fescue")

FescueOgCP <- ggplot(Data_g3p4_AdjCP, aes(x = Data_Set, y = as.numeric(MeanCP))) +
  geom_point(size = 2) +
  ggtitle("Origional CP Values Fescue")

ggarrange(EpiAdjCP, EpiOgCP, ncol = 1, nrow = 2)
ggarrange(FescueAdjCP, FescueOgCP, ncol = 1, nrow = 2)

########## Inputting the rest of the data (315, data not using red blue 2x2 standards)  ###############
###### Introducing the rest of the data 
###################################################################################################
# commands used to make data sets, then added data set column
# cat cp_values_1041-1120_2x2.csv cp_values_1121-1200_2x2.csv cp_values_1201-1280_2x2.csv cp_values_1281-END-Redoes_2x2.csv > 1041-end-redoes_2x2.csv
# cat cp_values_1041-1121_g3p4.csv cp_values_1121-1200_g3p4.csv cp_values_1201-1280_g3p4.csv cp_values_1281-END-Redoes_g3p4.csv > 1041-end-redoes_g3p4.csv
all_2x2_1 <- read.csv(all_2x2_loc, header = TRUE)
all_g3p4_1 <- read.csv(all_g3p4_loc, header = TRUE)
all_2x2_2 <- read.csv(all_2x2_2_loc, header = TRUE)
all_g3p4_2 <- read.csv(all_g3p4_2_loc, header = TRUE)
all_315x320_2x2 <- read.csv(all_315x320_2x2_loc, header = TRUE)
all_315x320_g3p4 <- read.csv(all_315x320_g3p4_loc, header = TRUE)

all_2x2_1$Standard <- "First_combo"
all_2x2_2$Standard <- "Second_combo"
all_g3p4_1$Standard <- "First_combo"
all_g3p4_2$Standard <- "Second_combo"
all_315x320_2x2$Standard <- "Pre_combo"
all_315x320_g3p4$Standard <- "Pre_combo"

all_2x2_samples <- rbind(all_2x2_1, all_2x2_2, all_315x320_2x2)
all_g3p4_samples <- rbind(all_g3p4_1, all_g3p4_2, all_315x320_g3p4)


#Removing Redone samples
redoes_loc <- "/home/drt06/Documents/QPCR_Data_Wrangler/Program/int_files/Data_for_Project/Sample_Redo_List.csv"
redoes <- read.csv(redoes_loc, header = FALSE)
redoner<- function(redoes,sampledata){
  for (x in nrow(reodes):1){
    for (y in nrow(sampledata):1){
      if(reodes[x,] == sampledata[y,3]){
        if(sampledata[y,7] != "1281-1325-R"){
          sampledata <- sampledata[-y,]
        }
      }  
    }
  }
  return(sampledata)
}

all_2x2_samples <- redoner(redoes,all_2x2_samples)
all_g3p4_samples <- redoner(redoes,all_g3p4_samples)

# Removing specific redone samples
all_2x2_samples <- all_2x2_samples[!(all_2x2_samples$Treatment == "303-1-41" & all_2x2_samples$Data_Set == "81-160" ),]

all_2x2_samples$Concentration <- as.numeric(all_2x2_samples$Concentration)
all_g3p4_samples$Concentration <- as.numeric(all_g3p4_samples$Concentration)

standard_group <- subset(all_2x2_samples, select = c("Data_Set","Standard"))
# Filtering and adjusting cp values fo the rest of the data
# takes out obviously bad samples, 0's, missing, ect.
Data_g3p4_Means_all <- Data_Filtering(all_g3p4_samples) # method filters data and leaves only samples
Data_2x2_Means_all <- Data_Filtering(all_2x2_samples)

# Method seperates standards from data,  to get their mean and adds copy number and log copy number
std_rest_g3p4_means <- FindStandardMeans(all_g3p4_samples,230)
# Needed to double filter 2x2 for some reason
for (x in nrow(all_2x2_samples):1){
  if (all_2x2_samples[x,2]=="" | all_2x2_samples[x,2]=="0" | nchar(all_2x2_samples[x,3])!=4){
    all_2x2_samples <- all_2x2_samples[-x,]
  }
}
std_rest_2x2_means <- FindStandardMeans(all_2x2_samples,118)

# method adjusts the CP value based on efficiency
# Also filters out samples that are off the standard curve, too low, too high, ect. 
Data_Sets <- unique(std_rest_g3p4_means$Data_Set) #Creating list of data sets for next method

Data_g3p4_AdjCP_rest <- CpAdjusterP2(std_rest_g3p4_means, Data_g3p4_Means_rest,Data_Sets) 
Data_2x2_AdjCP_rest <- CpAdjusterP2(std_rest_2x2_means, Data_2x2_Means_rest,Data_Sets) 


################ Combining data and getting the delta CT values  ###############

combined_2x2 <- rbind(Data_2x2_AdjCP,Data_2x2_AdjCP_rest)
combined_g3p4 <- rbind(Data_g3p4_AdjCP, Data_g3p4_AdjCP_rest)

#Combining data set and naming columns
all_Data <- merge(combined_2x2,combined_g3p4, by.x = c("Data_Set", "Treatment"), by.y = c("Data_Set", "Treatment"))
Columns <- colnames(all_Data)
Columns <- gsub("\\.y$", ".Fescue", Columns)
Columns <- gsub("\\.x$", ".Epichloe", Columns)
colnames(all_Data) <- Columns 
all_Data$adjCP.Fescue <- as.numeric(all_Data$adjCP.Fescue)
all_Data$adjCP.Epichloe <- as.numeric(all_Data$adjCP.Epichloe)

# Dela CT is Fescue - EPichloe
all_Data$Delta_CT <- all_Data$adjCP.Fescue - all_Data$adjCP.Epichloe

# Delta Ratio is Fescue / Epichloe
all_Data$Fes_to_Epi_Ratio <- all_Data$adjCP.Fescue / all_Data$adjCP.Epichloe

# Introducing Standard Combination Data
all_Data <- merge(all_Data,standard_group, by.x = c("Data_Set"), by.y =c("Data_Set"))
##### Making graphs to visualize the Fescue to Epi differences ######

ggplot(all_Data, aes(x=Data_Set, y=Fes_to_Epi_Ratio, color = Data_Set, shape = Standard)) +
  geom_point( size = 4) +
  ggtitle("Ratios of Fescue / Epichloe")

ggplot(all_Data, aes(x=Data_Set, y=Delta_CT, color = Data_Set)) +
  geom_point( size = 4) +
  ggtitle("Fescue CT - Epichloe CT")

ggplot(all_Data, aes(x=Data_Set, y=Delta_CT, color = Standard)) +
  geom_point( size = 4) +
  ggtitle("Fescue CT - Epichloe CT") +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1),
  panel.background = element_rect(fill = 'white', color = 'grey'))

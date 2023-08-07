# Objective : This code will combine all the data using my original standards find a correlation between them

# Import libraries
library(tidyverse)
library(gvlma)
library(mbQTL)

#Importing data sets
all_2x2_loc <- "/home/drt06/Documents/QPCR_Data_Wrangler/Program/int_files/Data_for_Project/all_2x2_Red_Blue.csv"
all_g3p4_loc <- "/home/drt06/Documents/QPCR_Data_Wrangler/Program/int_files/Data_for_Project/all_g3p4_correspond_to_Red_Blue_Standards.csv"

all_2x2 <- read.csv(all_2x2_loc, header = TRUE)
all_g3p4 <- read.csv(all_g3p4_loc, header = TRUE)


# Cleaning up Standard data 
# Remove any 0's from the stds data

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
            check_overlap = T)

ggplot(Graph_Data, aes(x=Tall_Fescue_Standard_CP_Values, y=Epichloe_Standard_CP_Values, color = Treatment, labels = Data_Set)) +
  geom_point( size = 2) +
  geom_text(label = Graph_Data$Data_Set, 
            nudge_x = 0.25, nudge_y = 0.25, 
            check_overlap = T)

# Remove samples 1-80 for now because something is wrong with the 2x2 part of it.
# Removing samples 1-80 made the model not work 
stds_2x2 <- stds_2x2[!stds_2x2$Data_Set == "1-80",]
stds_g3p4 <- stds_g3p4[!stds_g3p4$Data_Set == "1-80",]





# looking at different correlation methods and creating a model
# These models do not work I have to log transform them

# cor(stds_g3p4$Tall_Fescue_Standard_CP_Values,stds_2x2$Epichloe_Standard_CP_Values, method=c("spearman"))
# cor(stds_g3p4$Tall_Fescue_Standard_CP_Values,stds_2x2$Epichloe_Standard_CP_Values, method=c("pearson"))
# 
# model <- lm(stds_2x2$Epichloe_Standard_CP_Values ~ stds_g3p4$Tall_Fescue_Standard_CP_Values)
# summary(model)
# model2 <- lm(stds_g3p4$Tall_Fescue_Standard_CP_Values ~ stds_2x2$Epichloe_Standard_CP_Values)
# summary(model2)
# 
# gvlma(model)
# gvlma(model2) # THis model fails

# Retrying the models but with a quick log transformation 
stds_2x2$LogEpichloe_Standard_CP_Values <- log(stds_2x2$Epichloe_Standard_CP_Values)
stds_g3p4$LogTall_Fescue_Standard_CP_Values <- log(stds_g3p4$Tall_Fescue_Standard_CP_Values)

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
# Then I have to 
EpiPredictor <- 
  function(x){
    y <- y_intercept + TF_Coefficient*x
    CP <- 2.71828^y  
    return(CP)
      }

EpiPredictor(2.74)


################# Calculating efficiency for each data set and creating adjusted CP values #####################
Data_Sets <- unique(all_2x2$Data_Set)

Standards_2x2 <- all_2x2
Standards_g3p4 <- all_g3p4
all_2x2 <- read.csv(all_2x2_loc, header = TRUE)
all_g3p4 <- read.csv(all_g3p4_loc, header = TRUE)

Standards_2x2$Concentration <- as.numeric(Standards_2x2$Concentration)
Standards_g3p4$Concentration <- as.numeric(Standards_g3p4$Concentration)

# seperating the standards and getting the means of them, calculating ng of DNA
FindStandardMeans <- function(stdData, length){
  Standard_Means <-
    stdData %>%
    group_by(Data_Set,Treatment) %>%
    summarise(MeanCP = mean(Cp), NgDNA = mean(Concentration)*5) # we add 5 ul of sample
  Standard_Means$CopyNumber <- (Standard_Means$NgDNA * 6.02214076*10^23/ (length * 650 * 10^9))
  # NgDNA*1 mole / length * 1 mole of base pairs in ng
  Standard_Means$LogCopyNumber <- log10(Standard_Means$CopyNumber)
  return(Standard_Means)}

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
gvlma(modelEpiStandards)

# This function allows you to put in a CP value and it will output an adjusted CP value that takes efficiency into consideration
# Step 1 Calculate line for a single set of standards
# Step 2 Calculate efficiency for standards
# Step 3 Use efficiency to find adjusted CP value
Standard_means_2x2 <- FindStandardMeans(Standards_2x2,118)

CpAdjuster <- # This is here more just for convince, move it up later
  function(model, Data){
    Data$adjCP <- 0
    slope <- model$coefficients[2] # This is not the slope I should use, need to use log copy number slope
    E = -1+10^(-1/slope)
    print(slope) 
    print(E)
    
    for (i in 1:nrow(Data)){
      ogCP <- Data[i,3]
      ratio =  (2^ogCP) / (2*E)^ogCP # Finding ratio of expected vs actual
      adjuster = log(ratio, base = 2) # Finding amount needed to adjust CP videa
      adjCP = round(ogCP + adjuster, 3) #finalizes CP adjustment
      Data$adjCP[i] <- adjCP
    }
    
    return(Data)
  }

# Seperating Data From Standards and getting means
Data_2x2 <- all_2x2
Data_2x2$Cp <- as.numeric(Data_2x2$Cp)
for (x in nrow(Data_2x2):1){
  if (Data_2x2[x,2]=="" | Data_2x2[x,2]=="0" | nchar(all_g3p4[x,3])<7 | Data_2x2[x,4]=="EndoPos_Neg_Water" ){
    Data_2x2 <- Data_2x2[-x,]
  }
}
Data_2x2_Mean <-
  Data_2x2 %>%
  group_by(Data_Set,Treatment) %>%
  summarise(MeanCP = mean(Cp))
            
Data_g3p4 <- all_g3p4
Data_g3p4$Cp <- as.numeric(Data_g3p4$Cp)
for (x in nrow(Data_g3p4):1){
  if (Data_g3p4[x,2]=="" | Data_g3p4[x,2]=="0" | nchar(all_g3p4[x,3])<8 | Data_g3p4[x,4]=="EndoPos_Neg_Water" ){
    Data_g3p4 <- Data_g3p4[-x,]
  }
}
Data_g3p4_Mean <-
  Data_g3p4 %>%
  group_by(Data_Set,Treatment) %>%
  summarise(MeanCP = mean(Cp))

# Using Standards and Data to adjust CP Values 
Data_Sets <- as.matrix(Data_Sets)
All_Data = NULL
for (i in 1:nrow(Data_Sets)){
  
  Set <- Data_Sets[i]
  # Need to test if this loop works now
  Current_Standard <- subset(Standard_means_2x2, Data_Set == Set)
  Current_Data <- subset(Data_2x2_Mean, Data_Set == Set)
  
  # Something is wrong here that is causing my efficiency to not work. It could also be in the efficiency 
  # Equation itself.
  modelSet <- lm(Current_Standard$MeanCP ~ Current_Standard$LogCopyNumber)
  y_intercept <- modelSet$coefficients[1]
  slope <- modelSet$coefficients[2]
  summary(modelSet)
  Current_Data$LogCopyNumber <- round(((Current_Data$MeanCP)-y_intercept)/slope, digits = 3)
  Current_Data$CopyNumber <- 10^Current_Data$LogCopyNumber # Reverse the log
  
  All_DataC <- CpAdjuster(modelSet, Current_Data) 
  All_Data <- rbind(All_Data, All_DataC)
}

# Looking at the data to see if there are any descrepencies  
ggplot(Standard_means_2x2, aes(x=LogCopyNumber, y=MeanCP, color = Treatment, labels = Data_Set)) +
  geom_point( size = 2) +
  geom_text(label = Standard_means_2x2$Data_Set,
            check_overlap = T)




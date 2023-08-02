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
cor(stds_g3p4$Tall_Fescue_Standard_CP_Values,stds_2x2$Epichloe_Standard_CP_Values, method=c("spearman"))
cor(stds_g3p4$Tall_Fescue_Standard_CP_Values,stds_2x2$Epichloe_Standard_CP_Values, method=c("pearson"))
stds_2x2$LogEpichloe_Standard_CP_Values <- log(stds_2x2$Epichloe_Standard_CP_Values)
stds_g3p4$LogTall_Fescue_Standard_CP_Values <- log(stds_g3p4$Tall_Fescue_Standard_CP_Values)

model <- lm(stds_2x2$Epichloe_Standard_CP_Values ~ stds_g3p4$Tall_Fescue_Standard_CP_Values)
summary(model)
model2 <- lm(stds_g3p4$Tall_Fescue_Standard_CP_Values ~ stds_2x2$Epichloe_Standard_CP_Values)
summary(model2)

gvlma(model)
gvlma(model2) # THis model fails

# Retrying the models but with a quick log transformation 
stds_2x2$LogEpichloe_Standard_CP_Values <- log(stds_2x2$Epichloe_Standard_CP_Values)
stds_g3p4$LogTall_Fescue_Standard_CP_Values <- log(stds_g3p4$Tall_Fescue_Standard_CP_Values)

model_logged <- lm(stds_2x2$LogEpichloe_Standard_CP_Values ~ stds_g3p4$LogTall_Fescue_Standard_CP_Values)
summary(model_logged)
gvlma(model_logged)

Graph_Data_logged <- merge(stds_2x2, stds_g3p4, by.x = c("Data_Set", "Treatment"), by.y = c("Data_Set", "Treatment"))

ggplot(Graph_Data_logged, aes(x=LogTall_Fescue_Standard_CP_Values, y=LogEpichloe_Standard_CP_Values, color = Treatment, labels = Data_Set)) +
  geom_point( size = 2) +
  geom_text(label = Graph_Data_logged$Data_Set,
            check_overlap = T)

############################ This function uses that formula to predict 2x2 based on g3p4 numbers  ########################
# The formula is y-intercept + coefficient*value
EpiPredictor <- 
  function(x){5.31162 + 0.90407*x
  }
# This function does not work because it fails tests in the gvlma. 
# FescuePredictor <-
#   function(x){-4.45678 + 1.05716*x
#   }

EpiPredictor(21)

# Now we will use the concentration and Cp value of 2x2 to use new Cp values and assiagn them new concentrations.
# This assumes concentrations are wrong but keeps things consistant with my old standards
stds_2x2$LogConcentration <- 0
for (x in nrow(stds_2x2):1){
  if(stds_2x2[x,2] == "std1"){
    stds_2x2[x,4] = log(.01) 
  }else if (stds_2x2[x,2] == "std2"){
    stds_2x2[x,4] = log(.001) 
  }else if (stds_2x2[x,2] == "std3"){
    stds_2x2[x,4] = log(.0001)
  }else if (stds_2x2[x,2] == "std4"){
    stds_2x2[x,4] = log(.00001)
  }else if (stds_2x2[x,2] == "std5"){
    stds_2x2[x,4] = log(.000001)
  }
}

stds_2x2 <- stds_2x2[!stds_2x2$Data_Set == "1-80",]
stds_2x2 <- stds_2x2[!stds_2x2$Data_Set == "81-160",]
model_2x2_Conc <- lm(stds_2x2$LogConcentration ~ stds_2x2$Epichloe_Standard_CP_Values)
summary(model_2x2_Conc)
gvlma(model_2x2_Conc)

ggplot(stds_2x2, aes(x=LogConcentration, y=Epichloe_Standard_CP_Values, color = Data_Set)) +
  geom_point( size = 2) 


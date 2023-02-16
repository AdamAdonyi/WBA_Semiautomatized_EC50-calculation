#Automatization of Whole Blood Assay - 2022 Q4 - 

#Library:
#0_Load packages
#1_Load data
#2_rename
#3_Qc_transform class-as.double 
#4_add names
#5_ pivot longer & Values fitts to 0-100 range or EC50 calculation
#6_assign to treatment and remove std1 and Buffer - STING_A and STING_B
#7_separate treatment into Molecule, Unit, Concentration and add Molecule ID
#8_split() to separate diff citokines
#9_assign separated dataframs for each Analyte based on Treatment(Molecule)
#10_EC50_A
#Combination of A_EC50 values
#11_EC50_B
#Combination of B_EC50 values

rm(list = ls())
setwd("/Desktop/WBA automatized analysis/GitHub_upload/Real_data")


#0_Load packages----
library(tidyverse)
library(ggplot2)
library(readxl)
library(drc)


#1_Load data-----
train_data <- read_excel("2BM2HAMM84.xlsx")


#2_interaction-----
#View(train_data)

#3_Qc_transform class-as.double------
for (col in paste0("c", 1:12)) {
  train_data[[col]] <- as.double(train_data[[col]])
}



#4_add names------
train_data <- train_data %>% 
  rename(std = c(c1, c2)) %>% 
  rename(Buffer_然_0 = c3) %>% 
  rename(STING_然_1 = c4) %>% 
  rename(STING_然_3 = c5) %>% 
  rename(STING_然_5 = c6) %>% 
  rename(STING_然_8 = c7) %>% 
  rename(STING_然_10 = c8) %>% 
  rename(STING_然_15 = c9) %>% 
  rename(STING_然_30 = c10) %>% 
  rename(STING_然_60 = c11) %>% 
  rename(STING_然_100 = c12)



#5_ pivot longer & Values fitts to 0-100 range or EC50 calculation------
train_data <- train_data %>% 
  pivot_longer(c("std1", "std2", "Buffer_然_0", "STING_然_1", "STING_然_3", "STING_然_5", "STING_然_8", "STING_然_10", "STING_然_15", "STING_然_30", "STING_然_60", "STING_然_100"), names_to = "Treatment", values_to = "Values") %>% 
  mutate(Values = Values/1000)


#6_assign to treatment and remove std1 and Buffer - STING_A and STING_B------
all_STING_A <- train_data %>% 
  filter(Position %in% c("B", "C", "D") & grepl('STING', Treatment))  
  
all_STING_B <- train_data %>% 
  filter(Position %in% c("E", "F", "G") & grepl('STING', Treatment))


#7_separate treatment into Molecule, Unit, Concentration and add Molecule ID-----
all_STING_A <- all_STING_A %>% 
  separate(Treatment, into = c("Molecule", "Unit", "Concentration"), sep = "_") %>% 
  mutate(Molecule = 'BI A')

all_STING_B <- all_STING_B %>% 
  separate(Treatment, into = c("Molecule", "Unit", "Concentration"), sep = "_") %>% 
  mutate(Molecule = 'BI B')



#8_split to separate diff citokines----
all_STING_A_list <- split(all_STING_A, f=all_STING_A$Analyte)
all_STING_B_list <- split(all_STING_B, f=all_STING_B$Analyte)





#9_assign separated dataframs for each Analyte based on Treatment(Molecule) EVERYTHING MUST BE NUMERIC SPEC CONC!!!!!!! CHANGE------
##9.1_A_as dataframe-----
cols <- c("IL-10", "IL-6", "IFNg", "IFNa2a", "IP-10", "MCP-1", "IL-2", "GRZMb", "IFNb", "TNFa")

for (col in cols) {
  assign(paste0(gsub("-", "_", col), "_A"), as.data.frame(all_STING_A_list[[col]]))
}




###9.1.2_A_as.double-----
cols_2 <- c("IL_10_A", "IL_6_A", "IFNg_A", "IFNa2a_A", "IP_10_A", "MCP_1_A", "IL_2_A", "GRZMb_A", "IFNb_A", "TNFa_A")

for (col in cols_2) {
  assign(col, mutate(.data = get(col), Concentration = as.double(Concentration)))
}


##9.2_B_as dataframe (longer solution)-----
IL_10_B <- as.data.frame(all_STING_B_list$`IL-10`)
IL_6_B <- as.data.frame(all_STING_B_list$`IL-6`)
IFNg_B <- as.data.frame(all_STING_B_list$IFNg)
IFNa2a_B <- as.data.frame(all_STING_B_list$IFNa2a)
IP_10_B <- as.data.frame(all_STING_B_list$`IP-10`)
MCP_1_B <- as.data.frame(all_STING_B_list$`MCP-1`)
IL_2_B <- as.data.frame(all_STING_B_list$`IL-2`)
GRZMb_B <- as.data.frame(all_STING_B_list$GRZMb)
IFNb_B <- as.data.frame(all_STING_B_list$IFNb)
TNFa_B <- as.data.frame(all_STING_B_list$TNFa)

###9.2.2_B_as.double (longer solution)-----
IL_10_B$Concentration <- as.double(IL_10_B$Concentration)
IL_6_B$Concentration <- as.double(IL_6_B$Concentration)
IFNg_B$Concentration <- as.double(IFNg_B$Concentration)
IFNa2a_B$Concentration <- as.double(IFNa2a_B$Concentration)
IP_10_B$Concentration <- as.double(IP_10_B$Concentration)
MCP_1_B$Concentration <- as.double(MCP_1_B$Concentration)
IL_2_B$Concentration <- as.double(IL_2_B$Concentration)
GRZMb_B$Concentration <- as.double(GRZMb_B$Concentration)
IFNb_B$Concentration <- as.double(IFNb_B$Concentration)
TNFa_B$Concentration <- as.double(TNFa_B$Concentration)


#10_EC50_A-----
##10.1_IL_10_A----
curve_fit_IL_10_A <- drm(
  formula = Values ~ Concentration,
  data = IL_10_A,
  fct = LL.4(fixed = c(NA, NA, NA, NA), names = c('hill', 'min_value', 'max_value', 'ec_50'))
)

ec50_IL_10_A <- curve_fit_IL_10_A$coefficients['ec_50:(Intercept)']



##10.2_IL_6_A(adopted to the dataset)-----

IL_6_A <- na.omit(IL_6_A)

curve_fit_IL_6_A <- drm(
  formula = Values ~ Concentration,
  data = IL_6_A,
  fct = LL.4(fixed = c(NA, NA, NA, NA), names = c('hill', 'min_value', 'max_value', 'ec_50')),
  start = c(1, 0, max(IL_6_A$Values), median(IL_6_A$Concentration)),
  na.action = na.omit
)

ec50_IL_6_A <- curve_fit_IL_6_A$coefficients['ec_50:(Intercept)']



##10.3_IFNg_A------
curve_fit_IFNg_A  <- drm(
  formula = Values ~ Concentration,
  data = IFNg_A,
  fct = LL.4(fixed = c(NA, NA, NA, NA), names = c('hill', 'min_value', 'max_value', 'ec_50'))
)

ec50_IFNg_A  <- curve_fit_IFNg_A $coefficients['ec_50:(Intercept)']


##10.4_IFNa2a_A-----
curve_fit_IFNa2a_A <- drm(
  formula = Values ~ Concentration,
  data = IFNa2a_A,
  fct = LL.4(fixed = c(NA, NA, NA, NA), names = c('hill', 'min_value', 'max_value', 'ec_50'))
)

ec50_IFNa2a_A <- curve_fit_IFNa2a_A$coefficients['ec_50:(Intercept)']


##10.5_IP_10_A------
curve_fit_IP_10_A <- drm(
  formula = Values ~ Concentration,
  data = IP_10_A,
  fct = LL.4(fixed = c(NA, NA, NA, NA), names = c('hill', 'min_value', 'max_value', 'ec_50'))
)

ec50_IP_10_A <- curve_fit_IP_10_A$coefficients['ec_50:(Intercept)']


##10.6_MCP_1_A-------
curve_fit_MCP_1_A  <- drm(
  formula = Values ~ Concentration,
  data = MCP_1_A,
  fct = LL.4(fixed = c(NA, NA, NA, NA), names = c('hill', 'min_value', 'max_value', 'ec_50'))
)

ec50_MCP_1_A  <- curve_fit_MCP_1_A $coefficients['ec_50:(Intercept)']


##10.7_IL_2_A-----
curve_fit_IL_2_A <- drm(
  formula = Values ~ Concentration,
  data = IL_2_A,
  fct = LL.4(fixed = c(NA, NA, NA, NA), names = c('hill', 'min_value', 'max_value', 'ec_50'))
)

ec50_IL_2_A <- curve_fit_IL_2_A$coefficients['ec_50:(Intercept)']



##10.8_GRZMb_A-----
curve_fit_GRZMb_A <- drm(
  formula = Values ~ Concentration,
  data = GRZMb_A,
  fct = LL.4(fixed = c(NA, NA, NA, NA), names = c('hill', 'min_value', 'max_value', 'ec_50'))
)

ec50_GRZMb_A <- curve_fit_GRZMb_A$coefficients['ec_50:(Intercept)']


##10.9_IFNb_A-----
curve_fit_IFNb_A <- drm(
  formula = Values ~ Concentration,
  data = IFNb_A,
  fct = LL.4(fixed = c(NA, NA, NA, NA), names = c('hill', 'min_value', 'max_value', 'ec_50'))
)

ec50_IFNb_A <- curve_fit_IFNb_A$coefficients['ec_50:(Intercept)']



##10.10_TNFa_A------
curve_fit_TNFa_A <- drm(
  formula = Values ~ Concentration,
  data = TNFa_A,
  fct = LL.4(fixed = c(NA, NA, NA, NA), names = c('hill', 'min_value', 'max_value', 'ec_50'))
)

ec50_TNFa_A <- curve_fit_TNFa_A$coefficients['ec_50:(Intercept)']



#Combination of A_EC50 values----

EC50_A <- data.frame(
  cytokine = c("IL_6", "GRZMb", "IFNa", "IFNb", "IFNg", "IL_10", "IL_2", "IP_10", "MCP_1", "TNFa"),
  ec50_of_A_然= c(ec50_IL_6_A, ec50_GRZMb_A, ec50_IFNa2a_A, ec50_IFNb_A, ec50_IFNg_A, ec50_IL_10_A, ec50_IL_2_A, ec50_IP_10_A, ec50_MCP_1_A, ec50_TNFa_A)
)

View(EC50_A)

#11_EC50_B----
##11.1_IL_10_B----
curve_fit_IL_10_B <- drm(
  formula = Values ~ Concentration,
  data = IL_10_B,
  fct = LL.4(fixed = c(NA, NA, NA, NA), names = c('hill', 'min_value', 'max_value', 'ec_50'))
)

ec50_IL_10_B <- curve_fit_IL_10_B$coefficients['ec_50:(Intercept)']



##11.2_IL_6_B-----
curve_fit_IL_6_B <- drm(
  formula = Values ~ Concentration,
  data = IL_6_B,
  fct = LL.4(fixed = c(NA, NA, NA, NA), names = c('hill', 'min_value', 'max_value', 'ec_50'))
)

ec50_IL_6_B <- curve_fit_IL_6_B$coefficients['ec_50:(Intercept)']



##11.3_IFNg_B------
curve_fit_IFNg_B  <- drm(
  formula = Values ~ Concentration,
  data = IFNg_B,
  fct = LL.4(fixed = c(NA, NA, NA, NA), names = c('hill', 'min_value', 'max_value', 'ec_50'))
)

ec50_IFNg_B  <- curve_fit_IFNg_B $coefficients['ec_50:(Intercept)']


##11.4_IFNa2a_B-----
curve_fit_IFNa2a_B <- drm(
  formula = Values ~ Concentration,
  data = IFNa2a_B,
  fct = LL.4(fixed = c(NA, NA, NA, NA), names = c('hill', 'min_value', 'max_value', 'ec_50'))
)

ec50_IFNa2a_B <- curve_fit_IFNa2a_B$coefficients['ec_50:(Intercept)']


##11.5_IP_10_B------
curve_fit_IP_10_B <- drm(
  formula = Values ~ Concentration,
  data = IP_10_B,
  fct = LL.4(fixed = c(NA, NA, NA, NA), names = c('hill', 'min_value', 'max_value', 'ec_50'))
)

ec50_IP_10_B <- curve_fit_IP_10_B$coefficients['ec_50:(Intercept)']


##11.6_MCP_1_B-------
curve_fit_MCP_1_B  <- drm(
  formula = Values ~ Concentration,
  data = MCP_1_B,
  fct = LL.4(fixed = c(NA, NA, NA, NA), names = c('hill', 'min_value', 'max_value', 'ec_50'))
)

ec50_MCP_1_B  <- curve_fit_MCP_1_B $coefficients['ec_50:(Intercept)']


##11.7_IL_2_B-----
curve_fit_IL_2_B <- drm(
  formula = Values ~ Concentration,
  data = IL_2_B,
  fct = LL.4(fixed = c(NA, NA, NA, NA), names = c('hill', 'min_value', 'max_value', 'ec_50'))
)

ec50_IL_2_B <- curve_fit_IL_2_B$coefficients['ec_50:(Intercept)']



##11.8_GRZMb_B-----
curve_fit_GRZMb_B <- drm(
  formula = Values ~ Concentration,
  data = GRZMb_B,
  fct = LL.4(fixed = c(NA, NA, NA, NA), names = c('hill', 'min_value', 'max_value', 'ec_50'))
)

ec50_GRZMb_B <- curve_fit_GRZMb_B$coefficients['ec_50:(Intercept)']


##11.9_IFNb_B-----
curve_fit_IFNb_B <- drm(
  formula = Values ~ Concentration,
  data = IFNb_B,
  fct = LL.4(fixed = c(NA, NA, NA, NA), names = c('hill', 'min_value', 'max_value', 'ec_50'))
)

ec50_IFNb_B <- curve_fit_IFNb_B$coefficients['ec_50:(Intercept)']



##11.10_TNFa_B------
curve_fit_TNFa_B <- drm(
  formula = Values ~ Concentration,
  data = TNFa_B,
  fct = LL.4(fixed = c(NA, NA, NA, NA), names = c('hill', 'min_value', 'max_value', 'ec_50'))
)

ec50_TNFa_B <- curve_fit_TNFa_B$coefficients['ec_50:(Intercept)']





#Combination of B_EC50 values----

EC50_B <- data.frame(
  cytokine = c("IL_6", "GRZMb", "IFNa", "IFNb", "IFNg", "IL_10", "IL_2", "IP_10", "MCP_1", "TNFa"),
  ec50_of_B_然 = c(ec50_IL_6_B, ec50_GRZMb_B, ec50_IFNa2a_B, ec50_IFNb_B, ec50_IFNg_B, ec50_IL_10_B, ec50_IL_2_B, ec50_IP_10_B, ec50_MCP_1_B, ec50_TNFa_B)
)

View(EC50_B)


# EC50 calculation - final - 
#UPDATE: EC50: 2.091µM where GraphPadPrism calculated 2.02µM 

#GOOD AND WORKS PERFECTLY!
#Let use it on the real dataset



rm(list = ls())
setwd("C:/Users/adonyi/OneDrive - Boehringer Ingelheim/Desktop/WBA automatized analysis/")



#0_Load packages----
library(tidyverse)
library(ggplot2)
library(readxl)
library(drc)




#1_load real data-----
library(readxl)
data_measured <- read_excel("076_EC50_training_data.xlsx")

#Modify dataset formate to be prepared for the calculation
data_measured <- data_measured %>% 
  mutate(viable_cells = viable_cells/1000000) %>% 
  filter(replicate != "A")

#Interact with the dataset for double check---
View(data_measured)


#2_curve fit-----
curve_fit <- drm(
  formula = viable_cells ~ concentration,
  data = data_measured,
  fct = LL.4(fixed = c(NA, NA, NA, NA), names = c('hill', 'min_value', 'max_value', 'ec_50'))
)

ec50 <- curve_fit$coefficients['ec_50:(Intercept)']

curve_fit #hill/min/max/ec_50
ec50   #2.09099




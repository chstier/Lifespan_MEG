### This script creates a dataset containing cortical thickness values of all subjects of interest (CamCAN data) 
### and concatenates with frequency-specific global MEG connectivity
### Written and modified by Christina Stier, 2021/2022

### R Studio 3.6.3 (2020-02-29)
### RStudio 2021.09.0+351 "Ghost Orchid" Release (077589bcad3467ae79f318afe8641a1899a51606, 2021-09-20) for macOS
### Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) QtWebEngine/5.12.10 Chrome/69.0.3497.128 Safari/537.36

rm(list = ls())

## load packages 
install.packages("corrplot")
install.packages("igraph")
install.packages("qgraph")
install.packages("car")
install.packages("compute.es")
install.packages("effects")
install.packages("compute.es")
install.packages("ggplot2")
install.packages("multcomp")
install.packages("pastecs")
install.packages("psych")
install.packages("Hmisc")
install.packages("Rcmdr")
install.packages("splines")
install.packages("gridExtra")
install.packages("grid")	
install.packages("ggpubr")
install.packages("cowplot")
install.packages("optimx")
install.packages("plyr")
install.packages("doBy")
install.packages("boot")
install.packages("lmPerm")
install.packages('R.matlab')
install.packages('abind')
install.packages("readxl")

library(qgraph)
library(corrplot)
library(igraph)
require(igraph)
library(compute.es)
library(effects)
library(ggplot2)
library(multcomp)
library(pastecs)
library(WRS)
library(psych)
library(Hmisc)
library(car)
library(grid)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(optimx)
library(plyr)
library(doBy)
library(boot)
library(lmPerm)
library(R.matlab)
library(plyr)
library(abind)
library(reshape2)
library(tidyr)
library(readxl)

# set working directory
setwd("~/Documents/Projects/Age/R_camcan/files_pipeline/")

# load existing dataset
load("~/Documents/Projects/Age/R_camcan/datasets_pipeline/coh_img_data_absmean_global_allcohorts.Rda")

# view dataset
head(all_imcoh)
str(all_imcoh)

# load global thickness values
data = read.csv("~/Documents/Projects/Age/R_camcan/files_pipeline/metric_summary_05-Oct-2022.csv", header = TRUE, sep="\t")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
data = as.data.frame(data)
str(data)

# concatenate with global MEG values and rename variable
full_imcoh_th = cbind(all_imcoh, data$Global)
head(full_imcoh_th) 
names(full_imcoh_th)[15] = 'thickness' 

# get intracranial volume
icv = read.csv("~/Documents/Projects/Age/R_camcan/files_pipeline/icv_all_subjects.csv", header = FALSE)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
icv = as.data.frame(icv)
str(icv)

# concatenate with dataset and rename variable
full_imcoh_th = cbind(full_imcoh_th, icv$V1)
head(full_imcoh_th) 
names(full_imcoh_th)[16] = 'tiv' 

str(full_imcoh_th)
head(full_imcoh_th)

# save datasets 
setwd("~/Documents/Projects/Age/R_camcan/datasets_pipeline/")
save(full_imcoh_th, file = "coh_img_thickness_global_allcohorts.Rda") 

  
  
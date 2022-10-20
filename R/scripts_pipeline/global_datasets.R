### This script creates a dataset containing all subjects of interest (CamCAN data) including
### frequency-specific global values for MEG connectivity and power and demographic information
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

# load metric summary (global connectivity and power for each subject in the CamCAN-dataset)
global_files = Sys.glob("cohort*.csv") 

# load demographics for the respective cohort
demo_files = Sys.glob("Demo*.xlsx") 

datasetfolder = ("~/Documents/Projects/Age/R_camcan/datasets_pipeline")
cohorts = c("cohortX") # modify if needed

# create directory for global datasets
dir.create("~/Documents/Projects/Age/R_camcan/datasets_pipeline")

# loop over cohorts included
all_cohorts_imcoh = list()
all_cohorts_power = list()

for ( co in 1:length(demo_files)){	
  
  # get demographics of subjects per cohort
  f_name = demo_files[co]
  subj_data = as.data.frame(read_excel(f_name, col_names=FALSE))
  names(subj_data)[1] = 'subject_id'
  names(subj_data)[2] = 'age'
  names(subj_data)[3] = 'sex'
  names(subj_data)[4] = 'handedness'
  
  cohort = rep(co,50)
  
  # get global values per cohort and split for power and imcoh
  g_name = global_files[co]
  glob_data = as.data.frame(read.csv(g_name, header = TRUE, sep="\t"))
  data_power_absmean = subset(glob_data, glob_data$Metric == 'power' & glob_data$Stat == 'abs_mean')
  data_imcoh_absmean = subset(glob_data, glob_data$Metric == 'coh_img' & glob_data$Stat == 'abs_mean')
  
  #### now merge demographics and data for IMCOH
  data_imcoh_wide_mean= spread(data_imcoh_absmean, Freq, Global)
  
  sorted_demo = subj_data[match(data_imcoh_wide_mean$Subject_ID, subj_data$subject_id),]
  data_imcoh_full = cbind(sorted_demo, cohort, data_imcoh_wide_mean[,-(1)]) # get rid of ID-variable and cbind with metric_summary data frame
  
  # save full dataset for each cohort
  file_cohort = paste(datasetfolder, "/", cohorts[co], "_coh_img_data_absmean.Rda", sep="")
  save(data_imcoh_full, file = file_cohort) 
  
  all_cohorts_imcoh[[co]] = data_imcoh_full

  #### now merge demographics and data for POWER
  data_power_wide_mean= spread(data_power_absmean, Freq, Global)
  
  sorted_demo = subj_data[match(data_power_wide_mean$Subject_ID, subj_data$subject_id),]
  data_power_full = cbind(sorted_demo,cohort, data_power_wide_mean[,-(1)]) # get rid of ID-variable and cbind with metric_summary data frame
  
  # save full dataset for each cohort
  file_cohort = paste(datasetfolder, "/", cohorts[co], "_power_data_absmean.Rda", sep="")
  save(data_power_full, file = file_cohort) 
  
  all_cohorts_power[[co]] = data_power_full
}
  
all_imcoh = as.data.frame(do.call(rbind, all_cohorts_imcoh))  
all_power = as.data.frame(do.call(rbind, all_cohorts_power)) 

# save datasets 
setwd("~/Documents/Projects/Age/R_camcan/datasets_pipeline")
save(all_imcoh, file = "coh_img_data_absmean_global_allcohorts.Rda") 
save(all_power, file = "power_data_absmean_global_allcohorts.Rda")   
  
# for visualization of global value distribution use "global_plots_imcoh/power.R"  
### This script plots the correlation between age effects on MEG markers and cortical thickness
### Needs the output files from the statistical analysis using PALM (t-values fpr MEG connectivity/power and cortical thickness)
### Written and modified by Christina Stier, 2022

### R version 3.6.3 (2020-02-29)
### RStudio 2021.09.0+351 "Ghost Orchid" Release (077589bcad3467ae79f318afe8641a1899a51606, 2021-09-20) for macOS
### Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) QtWebEngine/5.12.10 Chrome/69.0.3497.128 Safari/537.36

rm(list = ls())

# install and load packages
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

library(qgraph)
library(corrplot)
library(igraph)
require(igraph)
library(compute.es)
library(effects)
library(ggplot2)
library(multcomp)
library(pastecs)
library(psych)
library(Hmisc)
library(car)
library(grid)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(lme4)
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

# create directory for storing the t-value maps
dir.create("~/Documents/Projects/Age/R_camcan/files_pipeline/vertexwise")
dir.create("~/Documents/Projects/Age/R_camcan/files_pipeline/vertexwise/dpv_tstat_age_age2_sex_tiv")

# Load all files in the folder
setwd("~/Documents/Projects/Age/R_camcan/files_pipeline/vertexwise/dpv_tstat_age_age2_sex_tiv")

filenames_imcoh = Sys.glob("*coh_img*.mat") 
filenames_power = Sys.glob("*power*.mat") 

# load t-statistics for thickness (here decrease with age: c2)
thickness_file = 'th_dpv_tstat_c2.mat'
thickness_filename = "t cortical thickness (decrease)" # change contrast if needed

# read cortical thickness file and name variable
t_th = data = readMat(thickness_file)
th = as.data.frame(array(unlist(t_th)))
names(th) = 't_th'

## set combination for correlations between age effects on MEG markers and cortical thickness (C2) that were significant after applying the spin test
# for connectivity
imcoh_set = c('Delta_coh_img_dpv_tstat_c1.mat')
label_imcoh = c('t delta connectivity (increase)')

# for power
power_set = c('Delta_power_dpv_tstat_c2.mat','Beta1_power_dpv_tstat_c4.mat', 'Beta2_power_dpv_tstat_c4.mat')
label_power = c('t delta power (decrease)', 't low beta power (∩-shape)', 't high beta power (∩-shape)')

# prepare lists
results_imcoh = list()
results_power = list()
plot_list_imcoh = list()
plot_list_power = list()

### plot for connectivity
for ( s in 1:length(imcoh_set)){
  
  func_name = imcoh_set[s]
  t_func = readMat(func_name)
  func = as.data.frame(array(unlist(t_func)))
  names(func) = 't_func'
  label_imcoh2 = label_imcoh[s]
  
  t_data = cbind(th, func)
  t_data = t_data[apply(t_data, 1, function(row) all (row !=0)), ]
  
  # plot correlation
  g = ggplot(t_data, aes(x = t_func, y = t_th)) + geom_point(size = 1) + stat_smooth(method = "lm", formula = y ~ x, se = FALSE) + labs(y=thickness_filename, x=label_imcoh2)
  g = g + theme(text = element_text(size = 23, family = "Calibri"), legend.title = element_blank(), panel.background = element_rect(fill = "white", colour = "grey50")) 
  g
  
  plot_list_imcoh[[s]] = g
  
  }

# plot and save
n = length(plot_list_imcoh)			 # save for all frequencies
nCol <- floor(sqrt(n))
finalplot = do.call("grid.arrange", c(plot_list_imcoh, ncol=3))
ggsave(file="/Users/christinastier/Documents/Projects/Age/R_camcan/plots_pipeline/Corr_ageeffects_imcoh_thickness.png", plot=finalplot, dpi = 150, limitsize = TRUE, width = 15, height = 6)

###  same for power sets
for ( s in 1:length(power_set)){
  
  func_name = power_set[s]
  t_func = readMat(func_name)
  label_power2 = label_power[s]
  func = as.data.frame(array(unlist(t_func)))
  names(func) = 't_func'
  
  t_data = cbind(th, func)
  t_data = t_data[apply(t_data, 1, function(row) all (row !=0)), ]
  
  # plot correlation
  g = ggplot(t_data, aes(x = t_func, y = t_th)) + geom_point(size = 1) + 
    # stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, aes(color = "black")) + 
    stat_smooth(method = "lm", formula = y ~ x, se = FALSE) + labs(y=thickness_filename, x=label_power2)
  g = g + theme(text = element_text(size = 23, family = "Calibri"), legend.title = element_blank(), panel.background = element_rect(fill = "white", colour = "grey50"), plot.margin=unit(c(0.5,0.5,0.5,0.5), "cm")) 
  g
  
  plot_list_power[[s]] = g
  
}

# plot and save
n = length(plot_list_power)			 # save for all frequencies
nCol <- floor(sqrt(n))
finalplot = do.call("grid.arrange", c(plot_list_power, ncol=3))
ggsave(file="/Users/christinastier/Documents/Projects/Age/R_camcan/plots_pipeline/Corr_ageeffects_power_thickness.png", plot=finalplot, dpi = 150, limitsize = TRUE, width = 15, height = 6)

####################  do the same with c4 contrast for cortical thickness (quadratic effects: inverted U)

# load t-statistics for thickness
thickness_file = 'th_dpv_tstat_c4.mat'
thickness_filename = "t cortical thickness (∩-shape)"

# read cortical thickness file and name variable
t_th = data = readMat(thickness_file)
th = as.data.frame(array(unlist(t_th)))
names(th) = 't_th'

# set combination for correlations between age effects on MEG markers and cortical thickness (C4) that were significant after applying the spin test
imcoh_set = c('Delta_coh_img_dpv_tstat_c3.mat', 'Alpha_coh_img_dpv_tstat_c2.mat', 'Beta1_coh_img_dpv_tstat_c2.mat')
label_imcoh = c('t delta connectivity (U-shape)', 't alpha connectivity (decrease)', 't low beta connectivity (decrease)')

results_imcoh = list()
plot_list_imcoh = list()

### compute for connectivity
for ( s in 1:length(imcoh_set)){
  
  func_name = imcoh_set[s]
  t_func = readMat(func_name)
  func = as.data.frame(array(unlist(t_func)))
  names(func) = 't_func'
  label_imcoh2 = label_imcoh[s]
  
  t_data = cbind(th, func)
  t_data = t_data[apply(t_data, 1, function(row) all (row !=0)), ]
  
  # plot correlation
  g = ggplot(t_data, aes(x = t_func, y = t_th)) + geom_point(size = 1) + 
    # stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, aes(color = "black")) + 
    stat_smooth(method = "lm", formula = y ~ x, se = FALSE) + labs(y=thickness_filename, x=label_imcoh2)
  g = g + theme(text = element_text(size = 22, family = "Calibri"), legend.title = element_blank(), panel.background = element_rect(fill = "white", colour = "grey50"), plot.margin=unit(c(0.5,0.5,0.5,0.5), "cm")) 
  g
  
  plot_list_imcoh[[s]] = g
  
}

# plot and save
n = length(plot_list_imcoh)			 # save for all frequencies
nCol <- floor(sqrt(n))
finalplot = do.call("grid.arrange", c(plot_list_imcoh, ncol=3))
ggsave(file="/Users/christinastier/Documents/Projects/Age/R_camcan/plots_pipeline/Corr_ageeffects_imcoh_thickness^2.png", plot=finalplot, dpi = 150, limitsize = TRUE, width = 15, height = 6)


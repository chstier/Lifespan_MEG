### This script plots the results for the predictor cortical thickness on MEG connectivity levels 
### (after regressing out the effects of sex, age, age2, intracranial vol. using the PALM toolbox)
### Written and modified by Christina Stier, 2022

### R version 3.6.3 (2020-02-29)
### RStudio 2021.09.0+351 "Ghost Orchid" Release (077589bcad3467ae79f318afe8641a1899a51606, 2021-09-20) for macOS
### Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) QtWebEngine/5.12.10 Chrome/69.0.3497.128 Safari/537.36

rm(list = ls())

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
#install.packages("remotes")
#remotes::install_github("ggseg/ggseg")
#install.packages('ggseg')

install.packages("remotes")
remotes::install_github("LCBC-UiO/ggsegYeo2011")

library(tidyverse)
library(qgraph)
library(corrplot)
library(igraph)
require(igraph)
library(compute.es)
library(effects)
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
library(ggseg)
library(ggplot2)
library(multcomp)
library(ggsegYeo2011)

# store results of PALM analysis for each functional network plus regnmaes_yeo_manual.mat-file in a folder and set directory accordingly:
setwd("~/Documents/Projects/Age/R_camcan/files_pipeline/yeo/palm_th_age_age2_sex_icv")

# load names of the networks
rois = readMat('regnames_yeo_manual.mat')
rois = as.data.frame(rois)
rois = as.data.frame(matrix(unlist(rois))) 

# extract information needed, name variables, and concatenate
hemi = as.data.frame(rois[17:32,1])
reg = as.data.frame(rois[1:16,])
region = as.data.frame(rois[33:48,])

roi_label = cbind(hemi, region, reg)
names(roi_label) = c('hemi', 'region', 'reg')
roi_label = roi_label[-c(1),] # exclude medial wall
roi_label = roi_label[-c(9),] # exclude medial wall

# choose which results you want to plot
# here an example of results for contrast 1 (positive prediction of MEG connectivity by cortical thickness after regressing out the co-variables)
# load file
file = '~/Documents/Projects/Age/R_camcan/files_pipeline/yeo/palm_th_age_age2_sex_icv/yeo_func_struct_allfreq_coh_img_c1.csv'
res = as.data.frame(read.csv(file, header = TRUE, sep="\t"))

# name frequency bands
freq = c('Delta', 'Theta', 'Alpha', 'Beta1','Beta2', 'Gamma') 
freq2 =  c('delta', 'theta', 'alpha', 'low beta','high beta', 'gamma') 

# prepare lists for storing the results 
# (t- and p-value for the predictor cortical thickness, corrected p (FDR) for the numbers of networks tested, filtered t: shows only t-values for significant networks)
plot_list_imcoh_t = list()
plot_list_imcoh_p = list()
plot_list_imcoh_pfdr = list()
plot_list_imcoh_t_filtered1 = list()
plot_list_imcoh_t_filtered2 = list()

for ( f in 1:length(freq)){
  freqname = freq[f]
  dat_freq = subset(res,res$freq == freqname)
  dat_freq = cbind(dat_freq, roi_label)
  
  # correct for number of networks
  p_fdr = p.adjust(dat_freq$p, method = "fdr", n = length(dat_freq$p))

  dat_freq = cbind(dat_freq, p_fdr)
  freqname2 = freq2[f]
  
  # plot original t-values
  g = ggplot(dat_freq) +
    geom_brain(atlas = yeo7, 
               position = position_brain(hemi ~ side),
               aes(fill = t)) +
    scale_fill_viridis_c(option = "cividis", direction = 1, limits = c(-3,3.25)) +
    theme_void() +
    labs(title = freqname2, size = 12)
  
  g1 = ggplot(dat_freq) +
    geom_brain(atlas = yeo7, 
               position = position_brain(hemi ~ side),
               aes(fill = -log10(p))) +
    scale_fill_viridis_c(option = "inferno", direction = 1, limits = c(1.3,3)) +
    theme_void() +
    labs(title = freqname2, size = 12)
  
  g11 = ggplot(dat_freq) +
    geom_brain(atlas = yeo7, 
               position = position_brain(hemi ~ side),
               aes(fill = -log10(p_fdr))) +
    scale_fill_viridis_c(option = "inferno", direction = 1, limits = c(1.3,1.7)) +
    theme_void() +
    labs(title = freqname2, size = 12)
  
  # plot filtered t (uncorrected)
  filtered1 = dat_freq
  filtered1$t[dat_freq$p > 0.05] = NA
  
  g_f1 = ggplot(filtered1) +
    geom_brain(atlas = yeo7, 
               position = position_brain(hemi ~ side),
               aes(fill = t)) +
    scale_fill_viridis_c(option = "cividis", direction = 1, limits = c(-3,3.25)) +
    theme_void() +
    labs(title = freqname2, size = 12)
  
  # plot filtered t (fdr corrected)
  filtered2 = dat_freq
  filtered2$t[dat_freq$p_fdr > 0.05] = NA
  
  g_f2 = ggplot(filtered2) +
    geom_brain(atlas = yeo7, 
               position = position_brain(hemi ~ side),
               aes(fill = t)) +
    scale_fill_viridis_c(option = "cividis", direction = 1, limits = c(-3,3.25)) +
    theme_void() +
    labs(title = freqname2, size = 12)
  
  
  plot_list_imcoh_t[[f]] = g
  plot_list_imcoh_p[[f]] = g1
  plot_list_imcoh_pfdr[[f]] = g11
  plot_list_imcoh_t_filtered1[[f]] = g_f1
  plot_list_imcoh_t_filtered2[[f]] = g_f2
}

# set working directory
setwd("~/Documents/Projects/Age/R_camcan/plots_pipeline")
 
# make plot for original t-values and save
n = length(plot_list_imcoh_t)			 # save for all frequencies
nCol <- floor(sqrt(n))
finalplot = do.call("grid.arrange", c(plot_list_imcoh_t, ncol=3))
ggsave(file="Imcoh_thickness_all_t_c1_300.png", plot=finalplot, dpi = 300, limitsize = TRUE, width = 12, height = 12)

# make plot for original p-values and save
n = length(plot_list_imcoh_p)			 # save for all frequencies
nCol <- floor(sqrt(n))
finalplot = do.call("grid.arrange", c(plot_list_imcoh_p, ncol=3))
ggsave(file="Imcoh_thickness_all_p_c1.png", plot=finalplot, dpi = 150, limitsize = TRUE, width = 12, height = 12)

# make plot for FDR-corrected p-values and save
n = length(plot_list_imcoh_pfdr)			 # save for all frequencies
nCol <- floor(sqrt(n))
finalplot = do.call("grid.arrange", c(plot_list_imcoh_pfdr, ncol=3))
ggsave(file="Imcoh_thickness_all_pfdr_c1.png", plot=finalplot, dpi = 150, limitsize = TRUE, width = 12, height = 12)  

# make plot for t-values in significant networks (uncorrected) and save
n = length(plot_list_imcoh_t_filtered1)			 # save for all frequencies
nCol <- floor(sqrt(n))
finalplot = do.call("grid.arrange", c(plot_list_imcoh_t_filtered1, ncol=3))
ggsave(file="Imcoh_thickness_all_t_filtered_c1_300.png", plot=finalplot, dpi = 300, limitsize = TRUE, width = 12, height = 12)

# make plot for t-values in significant networks (FDR corrected) and save
n = length(plot_list_imcoh_t_filtered2)			 # save for all frequencies
nCol <- floor(sqrt(n))
finalplot = do.call("grid.arrange", c(plot_list_imcoh_t_filtered2, ncol=3))
ggsave(file="Imcoh_thickness_all_t_fdrfiltered_c1.png", plot=finalplot, dpi = 150, limitsize = TRUE, width = 12, height = 12)

  
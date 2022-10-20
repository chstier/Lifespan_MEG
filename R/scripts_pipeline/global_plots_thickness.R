### This script plot global cortical thickness values across age ranges
### Written and modified by Christina Stier, 2021/2022

### R version 3.6.3 (2020-02-29)
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
#install.packages("Rcmdr")
install.packages("splines")
install.packages("gridExtra")
install.packages("grid")	
install.packages("ggpubr")
install.packages("cowplot")
install.packages("devtools")
install.packages("plotrix")
install.packages("gcookbook")

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
library(devtools)
library(reshape2)
library(tidyr)
library(plyr)
library(gcookbook) 
library(zoo)

# create directory for the plots
dir.create("~/Documents/Projects/Age/R_camcan/plots_pipeline")

# set working directory
setwd("~/Documents/Projects/Age/R_camcan/plots_pipeline")

# load global dataset
load("~/Documents/Projects/Age/R_camcan/datasets_pipeline/coh_img_thickness_global_allcohorts.Rda")

head(full_imcoh_th)
str(full_imcoh_th)

# change format of a few variables
full_imcoh_th$age = as.numeric(full_imcoh_th$age)
full_imcoh_th$sex = as.factor(full_imcoh_th$sex)
full_imcoh_th$handedness = as.numeric(full_imcoh_th$handedness)
full_imcoh_th$cohort = as.factor(full_imcoh_th$cohort)

dimnames(full_imcoh_th)
th_data = full_imcoh_th

# reorder data frame regarding frequencies
th_data = th_data[c(1,2,3,4,5,6,7,8,12,14,9,10,11,13, 15, 16)]  

# plot linear & quadratic effects of age on cortical thickness
name = "Cortical thickness" 
g = ggplot(th_data, aes(x = age, y = thickness)) + geom_point(aes(col=cohort), alpha = 1) + stat_smooth(method = "lm", formula = y ~ x + I(x^2), se= FALSE, size = 1.5, color = "gray") + 
  stat_smooth(method = "lm", formula = y ~ x, se = FALSE, size = 1, color = "#525252") + labs(y=name, x="Age")
g = g + scale_color_viridis_d(option = "inferno", guide = FALSE) 
g = g + scale_x_continuous(breaks=c(20,30,40,50,60,70,80))
g <- g + theme(text = element_text(size = 22, family = "Calibri"), legend.title = element_blank(), panel.background = element_rect(fill = "white", colour = "grey50"), plot.margin=unit(c(0.5,0.5,0.5,0.5), "cm"))	
g

ggsave(file="Age_linear_quadr_thickness_raw.png", plot=g, dpi = 300, limitsize = TRUE, width = 8, height = 6)	



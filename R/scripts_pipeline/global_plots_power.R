### This script plot global power values across age ranges for each frequency band separately
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
load("~/Documents/Projects/Age/R_camcan/datasets_pipeline/power_data_absmean_global_allcohorts.Rda")

# view the dataset
head(all_power)
str(all_power)

# change format of a few variables
all_power$age = as.numeric(all_power$age)
all_power$sex = as.factor(all_power$sex)
all_power$handedness = as.numeric(all_power$handedness)
all_power$cohort = as.factor(all_power$cohort)

# check and rename
dimnames(all_power)
power_data = all_power

# reorder data frame regarding frequencies change spelling
power_data = power_data[c(1,2,3,4,5,6,7,8,12,14,9,10,11,13)]  # reorder data frame regarding frequencies

names(power_data)[9] = 'delta'
names(power_data)[10] = 'theta'
names(power_data)[11] = 'alpha'
names(power_data)[12] = 'low beta'
names(power_data)[13] = 'high beta'
names(power_data)[14] = 'gamma'

# create long format 
data_long <- gather(power_data, freq, power, delta:gamma, factor_key=TRUE)

# take log10 of the power values for better visualization
data_long$power = log10(data_long$power)

# plot linear & quadratic age effects across the adult lifespan
freq = c('delta', 'theta','alpha', 'low beta', 'high beta','gamma') 	
plot_list = list()

for (i in 1:length(freq)){
  freqname = freq[i]
  data = data_long[data_long$freq == freqname,]
  
  #data$value = log(data$value)
  name = paste("Power (", freqname, ")", sep="") 
  g = ggplot(data, aes(x = age, y = power)) + geom_point(aes(col=cohort), alpha = 1) + stat_smooth(method = "lm", formula = y ~ x + I(x^2), se= FALSE, size = 1.5, color = "gray") + 
    stat_smooth(method = "lm", formula = y ~ x, se = FALSE, size = 1, color = "#525252") + labs(y=name, x="Age")
  g = g + scale_color_viridis_d(option = "inferno", guide = FALSE) 
  g = g + scale_x_continuous(breaks=c(20,30,40,50,60,70,80))
  g <- g + theme(text = element_text(size = 22, family = "Calibri"), legend.title = element_blank(), panel.background = element_rect(fill = "white", colour = "grey50"), plot.margin=unit(c(0.5,0.5,0.5,0.5), "cm"))	
  g
  plot_list[[i]] = g
}

n = length(plot_list)			 # save for all frequencies
nCol <- floor(sqrt(n))
finalplot = do.call("grid.arrange", c(plot_list, ncol=3))
ggsave(file="Age_linear_quadr_power.png", plot=finalplot, dpi = 300, limitsize = TRUE, width = 18, height = 10)	


# plot interaction with sex using facet grid
levels(data_long$sex)[1] = c('male')
levels(data_long$sex)[2] = c('female')

freq = c('delta', 'theta','alpha', 'low beta', 'high beta','gamma') 	
plot_list = list()

for (i in 1:length(freq)){
  freqname = freq[i]
  data = data_long[data_long$freq == freqname,]
  
  #data$value = log(data$value)
  name = paste("Power (", freqname, ")", sep="") 
  g = ggplot(data, aes(x = age, y = power)) + geom_point(aes(col=cohort), alpha = 1) + stat_smooth(method = "lm", formula = y ~ x + I(x^2), se= FALSE, size = 1.5, color = "gray") +
    stat_smooth(method = "lm", formula = y ~ x, se = FALSE, size = 1, color = "#525252") + labs(y=name, x="Age")
  g = g + scale_color_viridis_d(option = "inferno", guide = FALSE) 
  g = g + scale_x_continuous(breaks=c(20,30,40,50,60,70,80))
  g <- g + theme(text = element_text(size = 22, family = "Calibri"), legend.title = element_blank(), panel.background = element_rect(fill = "white", colour = "grey50"), plot.margin=unit(c(0.5,0.5,0.5,0.5), "cm")) + facet_grid(.~sex, switch="both")
  g
  plot_list[[i]] = g
}

n = length(plot_list)			 # save for all frequencies
nCol <- floor(sqrt(n))
finalplot = do.call("grid.arrange", c(plot_list, ncol=nCol))
ggsave(file="Age_linear_quadr_power_sex.png", plot=finalplot, dpi = 300, limitsize = TRUE, width = 12, height = 10)	


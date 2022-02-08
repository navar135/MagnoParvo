## This code solely tranforms MP data and computes analyses
## from the csv file "PminusMnorm_data_20200724.csv"
## By: Karen Navarro 
## 07/27/2020
####set environment####
install.packages("dplyr")
library(dplyr)
install.packages("reshape")
library(reshape)
install.packages("tidyr")
library(tidyr)
install.packages("car")
library(car)
install.packages("DescTools")
library(DescTools)
install.packages("tidyverse")
library(tidyverse)
install.packages("RColorBrewer")
library(RColorBrewer)
install.packages("plotly")
library(plotly)
install.packages("lsmeans")
library(lsmeans)
install.packages("ez")
library(ez)
install.packages("Hmisc")
library(Hmisc)
install.packages("sjstats")
library(sjstats)

####clear out current environment and set path to csv####
rm(list = ls(all.names = TRUE))
setwd("/Users/navar135/Documents/UMN_Research/Projects/MagnoParvo/Stimuli")


####Process MP data####
#Set & transform Dataframes w Veins
MPtrue <- read.csv("PminusMnorm_data_20200724.csv", sep = ",", header = TRUE)
#remove useless columns
MPtrue<-MPtrue[ , !(names(MPtrue) %in% "mask")]
#transformdata from wide to long for stats
temp <- melt(MPtrue, id=c("subj","voxels"))
#rename variables
names(temp) <- c("subject","voxels", "roi", "signal")
#filter data by eccentricity
FOV<- filter(temp,roi=="FOV_depth0" | roi=="FOV_depth2" | roi=="FOV_depth4"  | 
               roi=="FOV_depth6"  | roi=="FOV_depth8"  | roi=="FOV_depth10")
MID <- filter(temp,roi=="MID_depth0" | roi=="MID_depth2" | roi=="MID_depth4"  | 
                roi=="MID_depth6"  | roi=="MID_depth8"  | roi=="MID_depth10")
PER <-filter(temp,roi=="PER_depth0" | roi=="PER_depth2" | roi=="PER_depth4"  | 
               roi=="PER_depth6"  | roi=="PER_depth8"  | roi=="PER_depth10")

#create a 4th column for depths & rename roi column in each datafram
sequence<- seq(0, 10, by=2) #create depth labels
FOV$depth<-rep(sequence, each=10)
FOV$roi <- "FOV"
MID$depth<-rep(sequence, each=10)
MID$roi <- "MID"
PER$depth<-rep(sequence, each=10)
PER$roi <- "PER"
#combine all three dataframes
MP <-rbind(FOV,MID,PER)
## set all variables as factors and set levels for stats
## there are two factors (roi and depth) and 3 levels
## for roi and 6 levels for depth 
MP$roi <- as.factor(MP$roi)
MP$depth <- as.factor(MP$depth)
levels(MP$roi)<- list("FOV"=1,"MID"=2, "PER"=3)

#### Selected depths statistics ####
MPsub<- subset(MP, depth %in% c("0", "10"))

#check for sphericity violations
(sphericity<-ezANOVA(data=MPsub,
                     within=.(roi, depth),
                     wid=.(subject),
                     dv=.(signal),
                     return_aov = TRUE))
(justAnova2<-aov(signal~(roi *depth)+Error(subject/(roi * depth)),
    data=MPsub))
summary(justAnova2)
eta_sq(justAnova2, partial = TRUE)
anova_stats(justAnova2)
#perform 2way ANOVA 
justAnova1<-aov(signal ~ roi * depth,data=MPsub) #same results as model1
eta_sq(justAnova1)
summary(justAnova1)
anova_stats(justAnova1)
# different way to compute anova 
model1 <- lm(signal ~ roi*depth,data=MPsub) #this is what we are using in the paper! 
anova_stats(model1)

#perform some posthoc pairwise comparisons against other conditions
(TukeyHSD(justAnova1, "roi"))
(TukeyHSD(justAnova1, "depth"))
TukeyHSD(justAnova1,"roi:depth")
#different way to compute tukeysHSD (gives t-ratio) using model1
(int<-lsmeans(model1,
              pairwise ~ roi:depth,
              adjust="tukey"))
(tukROi<-lsmeans(model1,
                 pairwise ~ roi,
                 adjust="tukey"))
(tukdepth<-lsmeans(model1,
                   pairwise ~ depth,
                   adjust="tukey"))

#perform a bonferroni using model1
(bonRoi<-lsmeans(justAnova2,
                 pairwise ~ roi,
                 adjust="bonferroni"))
(bonDepth<-lsmeans(justAnova2,
                   pairwise ~ depth,
                   adjust="bonferroni"))
#different way to compute bonferroni using the data
(pairwise.t.test(MPsub$signal, MPsub$roi, p.adj = "bonferroni"))
(pairwise.t.test(MPsub$signal, MPsub$depth, p.adj = "bonferroni"))

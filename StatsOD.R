## By: Karen Navarro 
## 10/08/2019

## This code analyses Left Eye-Right Eye statistics from the csv file called OD_data.csv
#set environment
library(ggplot2)
library(dplyr)
library(reshape)
library(tidyr)
library(car)
library(DescTools)
setwd("/Users/navar135/Documents/UMN_Research/Projects/MagnoParvo/Stimuli")

####Set & transform Dataframes w/o Veins####
OD <- read.csv("OD_data.csv", sep = ",", header = TRUE)
#choose included datasets
included <-c("pnr161_20190124","pnr161_20190404","pnr492_20190204",
             "pnr492_20190429","pnr521_20190321","pnr739_20190131",
             "pnr739_20190415","pnr521_20190808","pnr328_20190808")
OD<- OD[OD$subj %in% included, ]
ODveins <- subset(OD,mask=="GMsigveins_75pct.nii")
#remove useless columns
ODveins<-ODveins[ , !(names(ODveins) %in% "mask")]
ODveins<- ODveins[ , !(names(ODveins) %in% "voxels")]
#transformdata from wide to long for stats
ODveins<- melt(ODveins)
#rename variables
names(ODveins) <- c("subject", "roi", "signal")
#filter data by eccentricity
V1<- filter(ODveins,roi=="V1_depth0" | roi=="V1_depth2" | roi=="V1_depth4"  | 
               roi=="V1_depth6"  | roi=="V1_depth8"  | roi=="V1_depth10")
V2 <- filter(ODveins,roi=="V2_depth0" | roi=="V2_depth2" | roi=="V2_depth4"  | 
                roi=="V2_depth6"  | roi=="V2_depth8"  | roi=="V2_depth10")

#crete a 4th column for depths & rename roi column in each datafram
sequence<- seq(0, 10, by=2)
V1$depth<-rep(sequence, each=9)
V1$roi <- "V1"
V2$depth<-rep(sequence, each=9)
V2$roi <- "V2"
#combine all three dataframes
ODveins <-rbind(V1,V2)

###PLOT STUFF###
library(tidyverse)
#
sumV1<- ODveins %>%
  group_by(depth) %>% 
  summarise(average = mean(signal),sd=sd(signal))
seV1 <- sumV1$sd/sqrt(length(included))
#standardize depths
#sumV1$depth<-c(0,2,4,6,8,10)
sumV1$depth<-sumV1$depth/max(sumV1$depth)
#put it into 1 dataframe by depths
drawnRL_V1<-data.frame(sumV1,seV1)

####Make Plots####
library(ggplot2)
library(RColorBrewer)
library(plotly)
## plot data collapsed across eccentricity
gg2<- ggplot(drawnRL_V1[order(drawnRL_V1$depth),], 
             aes(x=depth,
                 ymin = average-seV1,ymax = average+seV1, 
                 lty = '          V1 \n (PF + Mid + Per)'))
#plot shaded SE and flip coordinate system
gg2<- gg2 + geom_ribbon(alpha=.4,fill="#41ae76") + 
  coord_flip() 
#plot line and points and set size
gg2<- gg2+geom_point(aes(y = average),size=2,color="#41ae76") 
gg2<- gg2+ geom_line(aes(y = average), size=1, color="#41ae76")
#change axis limits
gg2<- gg2+ scale_x_continuous(limits=c(0,1),breaks = c(0,.2,.4,.6,.8,1))
gg2<- gg2+ scale_y_continuous(limits=c(0,.3),
                              breaks = c(0,.1,.2,.3))
#change labels name, size, and font 
gg2<- gg2 + labs(y="% Signal change \n |L-R|/(L+R)",
                 x= "Relative (equivolume) distance from WM")
gg2<-gg2 + theme(panel.background = element_rect(fill = "white", colour = "grey50"),
                 legend.title = element_blank(),
                 axis.text = element_text(colour = "black", 
                                          size = 16),
                 legend.text = element_text(size=16))
gg2<-gg2 + theme(axis.title.x = element_text(size=20)) #x title
gg2<-gg2 + theme(axis.title.y = element_text(size =20)) # y title
gg2
ggsave("RL_depthLaminarV1.pdf", device="pdf",width = 6, height = 8)

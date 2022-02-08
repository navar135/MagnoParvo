## This code analyses Magno-Parvo statistics from the csv files that are both 
## normilized and not-normalized
## and create laminar profiles for all data types
## By: Karen Navarro 
## 07/27/2020

####set environment####
install.packages("ggplot2")
library(ggplot2)
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
#clear out current environment
rm(list = ls(all.names = TRUE))
setwd("/Users/navar135/Documents/UMN_Research/Projects/MagnoParvo/Stimuli")

####Process MP data####
#Set & transform Dataframes w/o Veins
MPtrue <- read.csv("PM_data_20200724.csv", sep = ",", header = TRUE)
#remove useless columns
MPtrue<-MPtrue[ , !(names(MPtrue) %in% "mask")]

#transformdata from wide to long for stats
temp <- melt(MPtrue, id=c("subj","voxels"))
#temp<- melt(MPtrue)

#rename variables
names(temp) <- c("subject","voxels", "roi", "signal")
#filter data by eccentricity
FOV<- filter(temp,roi=="FOV_depth0" | roi=="FOV_depth2" | roi=="FOV_depth4"  | 
               roi=="FOV_depth6"  | roi=="FOV_depth8"  | roi=="FOV_depth10")
MID <- filter(temp,roi=="MID_depth0" | roi=="MID_depth2" | roi=="MID_depth4"  | 
                roi=="MID_depth6"  | roi=="MID_depth8"  | roi=="MID_depth10")
PER <-filter(temp,roi=="PER_depth0" | roi=="PER_depth2" | roi=="PER_depth4"  | 
               roi=="PER_depth6"  | roi=="PER_depth8"  | roi=="PER_depth10")
ALL <-filter(temp,roi=="ALL_depth0" | roi=="ALL_depth2" | roi=="ALL_depth4"  | 
               roi=="ALL_depth6"  | roi=="ALL_depth8"  | roi=="ALL_depth10")
#create a 5th column for depths & rename roi column in each datafram
sequence<- seq(0, 10, by=2)
FOV$depth<-rep(sequence, each=10)
FOV$roi <- "FOV"
MID$depth<-rep(sequence, each=10)
MID$roi <- "MID"
PER$depth<-rep(sequence, each=10)
PER$roi <- "PER"
ALL$depth <- rep(sequence, each=10)
ALL$roi <- "ALL"
#combine all three dataframes
MP <-rbind(FOV,MID,PER,ALL)
## set all variables as factors and set levels for stats
## there are two factors (roi and depth) and 3 levels
## for roi and 6 levels for depth 
MP$roi <- as.factor(MP$roi)
MP$depth <- as.factor(MP$depth)
levels(MP$roi)<- list("FOV"=1,"MID"=2, "PER"=3,"ALL"=4)


####Make Plots####
#group data and compute mean and sd for veins mask
sumFOV<- FOV %>% 
  group_by(depth) %>% 
  dplyr::summarise(average = mean(signal),sd=sd(signal))
sumMID<-MID %>% 
  group_by(depth) %>% 
  dplyr::summarise(average = mean(signal),sd=sd(signal))
sumPER <- PER %>% 
  group_by(depth) %>% 
  dplyr::summarise(average = mean(signal),sd=sd(signal))
sumALL <- ALL %>% 
  group_by(depth) %>% 
  dplyr::summarise(average = mean(signal),sd=sd(signal))
sumSignal <- rbind(sumFOV, sumMID, sumPER)
roinew <- c("FOV","FOV","FOV","FOV","FOV","FOV",
            "MID","MID","MID","MID","MID","MID",
            "PER","PER","PER","PER","PER","PER")

#compute standard error
se <- sumSignal$sd/sqrt(length(MPtrue$subj))
#standardize depths
sumSignal$depth<-sumSignal$depth/max(sumSignal$depth)
#put it into 1 dataframe by roi
laminar<-data.frame(sumSignal,roinew,se)
##Plot data across eccentricity depth vs % signal change
#plot data and switch default axis for ploting SE
gg<- ggplot(laminar[order(laminar$depth),], 
            aes(x=depth,ymin = average-se,ymax = average+se, 
                group=roinew))
#plot shaded SE and flip coordinate system
gg<- gg + geom_ribbon(aes(fill=roinew),alpha=.4) + 
  coord_flip() 
#plot lines and pick colors
gg<- gg+geom_point(aes(y = average,color = roinew),size=2) 
  #add a second layer of data to the graph - here we plot single subject data points 
  #geom_point(data = MPplot,aes(y = average,color = roinew),size=1)
gg<- gg+ geom_line(aes(y = average,color=roinew), size=1) 
#change axis limits
gg<- gg+ scale_x_continuous(limits=c(0,1),breaks = c(0,.2,.4,.6,.8,1))
gg<- gg+ scale_y_continuous(limits=c(-.2,.5),
                            breaks = c(-.2,-.1,0,.1,.2,.3,.4,.5,.6))
#change labels name, size, and font 
gg<- gg + labs(y="% Signal change \n PM",
               x= "Relative (equivolume) distance from WM")
gg<-gg + theme(panel.background = element_rect(fill = "white", colour = "grey50"),
               legend.title = element_blank(),
               axis.text = element_text(colour = "black", 
                                        size = 10),
               legend.text = element_text(size=10))
gg<-gg + theme(axis.title.x = element_text(size=12)) #x title
gg<-gg + theme(axis.title.y = element_text(size =12)) # y title
#change line and fill colors & legend labels
gg<-gg+scale_colour_manual(values = c("#045a8d" ,"#3690c0", "#a6bddb","#2171b5"), 
                           labels = c("Parafovea", "Middle", "Periphery","V1")
                           #values = c("0.3","0.6","0.9")
)
gg<- gg+scale_fill_manual(values = c("#045a8d" ,"#3690c0", "#a6bddb","#2171b5"), 
                          labels = c("Parafovea", "Middle", "Periphery","V1")
                          #values = c("0.3","0.6","0.9")
)

gg

#### plot just M vs P data ####
#Set & transform Dataframes w/o Veins
Ptrue <- read.csv("P_data_20200724.csv", sep = ",", header = TRUE)
#remove useless columns
Ptrue<-Ptrue[ , !(names(Ptrue) %in% "mask")]

#transform data from wide to long for stats
temp <- melt(Ptrue, id=c("subj","voxels"))
#temp<- melt(MPtrue)

#rename variables
names(temp) <- c("subject","voxels", "roi", "signal")
#filter data by eccentricity
FOV<- filter(temp,roi=="FOV_depth0" | roi=="FOV_depth2" | roi=="FOV_depth4"  | 
               roi=="FOV_depth6"  | roi=="FOV_depth8"  | roi=="FOV_depth10")
MID <- filter(temp,roi=="MID_depth0" | roi=="MID_depth2" | roi=="MID_depth4"  | 
                roi=="MID_depth6"  | roi=="MID_depth8"  | roi=="MID_depth10")
PER <-filter(temp,roi=="PER_depth0" | roi=="PER_depth2" | roi=="PER_depth4"  | 
               roi=="PER_depth6"  | roi=="PER_depth8"  | roi=="PER_depth10")
ALL <-filter(temp,roi=="ALL_depth0" | roi=="ALL_depth2" | roi=="ALL_depth4"  | 
               roi=="ALL_depth6"  | roi=="ALL_depth8"  | roi=="ALL_depth10")
#create a 5th column for depths & rename roi column in each datafram
sequence<- seq(0, 10, by=2)
FOV$depth<-rep(sequence, each=10)
FOV$roi <- "FOV"
MID$depth<-rep(sequence, each=10)
MID$roi <- "MID"
PER$depth<-rep(sequence, each=10)
PER$roi <- "PER"
ALL$depth <- rep(sequence, each=10)
ALL$roi <- "ALL"
#combine all three dataframes
P <-rbind(FOV,MID,PER,ALL)

## do the same for magno 
Mtrue <- read.csv("M_data_20200724.csv", sep = ",", header = TRUE)
#remove useless columns
Mtrue<-Mtrue[ , !(names(Mtrue) %in% "mask")]

#transformdata from wide to long for stats
temp <- melt(Mtrue, id=c("subj","voxels"))
#temp<- melt(MPtrue)

#rename variables
names(temp) <- c("subject","voxels", "roi", "signal")
#filter data by eccentricity
mFOV<- filter(temp,roi=="FOV_depth0" | roi=="FOV_depth2" | roi=="FOV_depth4"  | 
               roi=="FOV_depth6"  | roi=="FOV_depth8"  | roi=="FOV_depth10")
mMID <- filter(temp,roi=="MID_depth0" | roi=="MID_depth2" | roi=="MID_depth4"  | 
                roi=="MID_depth6"  | roi=="MID_depth8"  | roi=="MID_depth10")
mPER <-filter(temp,roi=="PER_depth0" | roi=="PER_depth2" | roi=="PER_depth4"  | 
               roi=="PER_depth6"  | roi=="PER_depth8"  | roi=="PER_depth10")
mALL <-filter(temp,roi=="ALL_depth0" | roi=="ALL_depth2" | roi=="ALL_depth4"  | 
               roi=="ALL_depth6"  | roi=="ALL_depth8"  | roi=="ALL_depth10")
#create a 5th column for depths & rename roi column in each datafram
sequence<- seq(0, 10, by=2)
mFOV$depth<-rep(sequence, each=10)
mFOV$roi <- "FOV"
mMID$depth<-rep(sequence, each=10)
mMID$roi <- "MID"
mPER$depth<-rep(sequence, each=10)
mPER$roi <- "PER"
mALL$depth <- rep(sequence, each=10)
mALL$roi <- "ALL"
#combine all three dataframes
M <-rbind(mFOV,mMID,mPER,mALL)
#add a column with M and P labels and combine both data frames so they can be plotted 

#group data and compute mean and sd for veins mask for P and M activation separately 
sumFOVp<- FOV %>% 
  group_by(depth) %>% 
  dplyr::summarise(average = mean(signal),sd=sd(signal))
sumMIDp<-MID %>% 
  group_by(depth) %>% 
  dplyr::summarise(average = mean(signal),sd=sd(signal))
sumPERp <- PER %>% 
  group_by(depth) %>% 
  dplyr::summarise(average = mean(signal),sd=sd(signal))
sumALLp <- ALL %>% 
  group_by(depth) %>% 
  dplyr::summarise(average = mean(signal),sd=sd(signal))
sumSignalP <- rbind(sumFOVp, sumMIDp, sumPERp)
roinew <- c("FOV","FOV","FOV","FOV","FOV","FOV",
            "MID","MID","MID","MID","MID","MID",
            "PER","PER","PER","PER","PER","PER")
#standardize depths
sumSignalP$depth<-sumSignalP$depth/max(sumSignalP$depth)
sumSignalP$cond<-rep("P", each=18)
#standard error
se <- sumSignalP$sd/sqrt(length(P$subj))
#put it into 1 dataframe by depths
Pplot<-data.frame(sumSignalP,roinew,se)

#group data and compute mean and sd for veins mask for P and M activation separately 
sumFOVm<- mFOV %>% 
  group_by(depth) %>% 
  dplyr::summarise(average = mean(signal),sd=sd(signal))
sumMIDm<-mMID %>% 
  group_by(depth) %>% 
  dplyr::summarise(average = mean(signal),sd=sd(signal))
sumPERm <- mPER %>% 
  group_by(depth) %>% 
  dplyr::summarise(average = mean(signal),sd=sd(signal))
sumALLm <- mALL %>% 
  group_by(depth) %>% 
  dplyr::summarise(average = mean(signal),sd=sd(signal))
sumSignalM <- rbind(sumFOVm, sumMIDm, sumPERm)

#standardize depths
sumSignalM$depth<-sumSignalM$depth/max(sumSignalM$depth)
sumSignalM$cond<-rep("M", each=18)
#standard error
se <- sumSignalM$sd/sqrt(length(M$subj))
#put it into 1 dataframe by depths
Mplot<-data.frame(sumSignalM,roinew,se)
PMplot<- rbind(Pplot,Mplot)
PMplot$roinew <- factor(PMplot$roinew, levels = c("FOV", "MID","PER"),
                  labels = c("Parafovea","Middle","Periphery")
)
##Plot data across eccentricity depth vs % signal change
#plot data and switch default axis for ploting SE
gg1<- ggplot(PMplot[order(PMplot$depth),], 
             aes(x=depth,
                 ymin = average-se,ymax = average+se,
                 group=cond
             ))
gg1<-gg1+facet_wrap( ~ roinew, nrow = 1) 
#plot shaded SE and flip coordinate system
gg1<- gg1 + geom_ribbon(aes(fill=cond),alpha=.4) + 
  coord_flip() 
#plot lines and pick colors
gg1<- gg1+geom_point(aes(y = average,color = cond),size=.75) 
gg1<- gg1+ geom_line(aes(y = average,color=cond), size=.25) 
#change axis limits
gg1<- gg1+ scale_x_continuous(limits=c(0,1),breaks = c(0,.2,.4,.6,.8,1))
gg1<- gg1+ scale_y_continuous(limits=c(0.8,4),
                            breaks = c(1,1.5,2,2.5,3,3.5))
#change labels name, size, and font 
gg1<- gg1 + labs(y="% Signal change",
               x= "Relative (equivolume) distance from WM")
gg1<-gg1 + theme(panel.background = element_rect(fill = "#ffffff", colour = "grey50"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
               legend.title = element_blank(),
               axis.text = element_text(colour = "black", 
                                        size = 10),
               legend.text = element_text(size=10))
gg1<-gg1 + theme(axis.title.x = element_text(size=12)) #x title
gg1<-gg1 + theme(axis.title.y = element_text(size =12)) # y title
#change line and fill colors & legend labels
gg1<-gg1+scale_colour_manual(values = c("#253494","#d7da42")
                             #labels = c("M-targeted ","P-targeted")
)
gg1<- gg1+scale_fill_manual(values = c("#253494","#d7da42")
)
gg1
ggsave("depthLaminar_PvsM.pdf", device="pdf",width = 8, height = 4)

FOVtest<- filter(PMplot,roinew=="Parafovea")
MIDtest<- filter(PMplot,roinew=="Middle")
PERtest<-filter(PMplot,roinew=="Periphery")
ALLtest<- filter(PMplot,roinew=="V1")

par(mfrow=c(1,3))
gg2<- ggplot(FOVtest[order(FOVtest$depth),], 
             aes(x=depth,
                 ymin = average-se,ymax = average+se,group=cond
             ))
#plot shaded SE and flip coordinate system
gg2<- gg2 + geom_ribbon(aes(fill=cond),alpha=.4) + 
  coord_flip() 
#plot lines and pick colors
gg2<- gg2+geom_point(aes(y = average,color=cond),size=1) 
gg2<- gg2+ geom_line(aes(y = average,color=cond), size=0.5) 

gg3<- ggplot(MIDtest[order(MIDtest$depth),], 
             aes(x=depth,
                 ymin = average-se,ymax = average+se,group=cond
             ))
#plot shaded SE and flip coordinate system
gg3<- gg3 + geom_ribbon(aes(fill=cond),alpha=.4) + 
  coord_flip() 
#plot lines and pick colors
gg3<- gg3+geom_point(aes(y = average,color=cond),size=1) 
gg3<- gg3+ geom_line(aes(y = average,color=cond), size=0.5) 


gg4<- ggplot(PERtest[order(PERtest$depth),], 
             aes(x=depth,
                 ymin = average-se,ymax = average+se,group=cond
             ))
gg4<-gg4+facet_wrap( ~ roinew , nrow = 1) 
#plot shaded SE and flip coordinate system
gg4<- gg4 + geom_ribbon(aes(fill=cond),alpha=.4) + 
  coord_flip() 
#plot lines and pick colors
gg4<- gg4+geom_point(aes(y = average,color=cond),size=1) 
gg4<- gg4+ geom_line(aes(y = average,color=cond), size=0.5) 

gg5<- ggplot(PERtest[order(PERtest$depth),], 
             aes(x=depth,
                 ymin = average-se,ymax = average+se,group=cond
             ))
gg5<-gg5+facet_wrap( ~ roinew , nrow = 1) 
#plot shaded SE and flip coordinate system
gg5<- gg5 + geom_ribbon(aes(fill=cond),alpha=.4) + 
  coord_flip() 
#plot lines and pick colors
gg5<- gg5+geom_point(aes(y = average,color=cond),size=1) 
gg5<- gg5+ geom_line(aes(y = average,color=cond), size=0.5) 


gg2
gg3
gg4

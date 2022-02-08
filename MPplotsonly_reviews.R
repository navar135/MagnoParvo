## This code solely plots laminar profiles and accessory plots from tranformed MP data 
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


####Process P and M data separate####
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
pFOV<- filter(temp,roi=="FOV_depth0" | roi=="FOV_depth2" | roi=="FOV_depth4"  | 
               roi=="FOV_depth6"  | roi=="FOV_depth8"  | roi=="FOV_depth10")
pMID <- filter(temp,roi=="MID_depth0" | roi=="MID_depth2" | roi=="MID_depth4"  | 
                roi=="MID_depth6"  | roi=="MID_depth8"  | roi=="MID_depth10")
pPER <-filter(temp,roi=="PER_depth0" | roi=="PER_depth2" | roi=="PER_depth4"  | 
               roi=="PER_depth6"  | roi=="PER_depth8"  | roi=="PER_depth10")
pALL <-filter(temp,roi=="ALL_depth0" | roi=="ALL_depth2" | roi=="ALL_depth4"  | 
               roi=="ALL_depth6"  | roi=="ALL_depth8"  | roi=="ALL_depth10")
#create a 5th column for depths & rename roi column in each datafram
sequence<- seq(0, 10, by=2)
pFOV$depth<-rep(sequence, each=10)
pFOV$roi <- "FOV"
pMID$depth<-rep(sequence, each=10)
pMID$roi <- "MID"
pPER$depth<-rep(sequence, each=10)
pPER$roi <- "PER"
pALL$depth <- rep(sequence, each=10)
pALL$roi <- "ALL"
#combine all three dataframes
P <-rbind(pFOV,pMID,pPER,pALL)

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
sumFOVp<- pFOV %>% 
  group_by(depth) %>% 
  dplyr::summarise(average = mean(signal),sd=sd(signal))
sumMIDp<-pMID %>% 
  group_by(depth) %>% 
  dplyr::summarise(average = mean(signal),sd=sd(signal))
sumPERp <- pPER %>% 
  group_by(depth) %>% 
  dplyr::summarise(average = mean(signal),sd=sd(signal))
sumALLp <- pALL %>% 
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
se <- sumSignalP$sd/sqrt(length(Ptrue$subj))
#put it into 1 dataframe by depths
Pplot<-data.frame(sumSignalP,roinew,se)

#quick stats to verify significance of plot 2 
PAovDat<- filter(P,roi !="ALL")
Paov<-aov(signal~roi,data=PAovDat)
summary(Paov)
PttestDat<- filter(PAovDat,depth == 0 | 
                    # depth == 2 | 
                    # depth == 8 |
                     depth == 10 )
PttestDat<- filter(PttestDat,roi != "MID" )
t.test(signal~roi,data=PttestDat,paired=T)

MAovDat<- filter(M,roi !="ALL")
Maov<-aov(signal~roi,data=MAovDat)
summary(Maov)
MttestDat<- filter(MAovDat,depth == 0 | 
                     # depth == 2 | 
                     # depth == 8 |
                     depth == 10 )
MttestDat<- filter(MttestDat,roi != "MID")
t.test(signal~roi,data=MttestDat,paired=T)

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
se <- sumSignalM$sd/sqrt(length(Mtrue$subj))
#put it into 1 dataframe by depths
Mplot<-data.frame(sumSignalM,roinew,se)
PMplot<- rbind(Pplot,Mplot)
PMplot$roinew <- factor(PMplot$roinew, levels = c("FOV", "MID","PER"),
                        labels = c("parafovea","middle","periphery"))
PMplot$cond <- factor(PMplot$cond, levels = c("P","M"),
                        labels = c("chromatic","achromatic"))

##Plot data across eccentricity depth vs % signal change
#plot data and switch default axis for ploting SE
pp1<- ggplot(PMplot[order(PMplot$depth),], 
             aes(x=depth,
                 ymin = average-se,ymax = average+se,
                 group=cond
             ))
pp1<- pp1+facet_wrap( ~ roinew, nrow = 1) 
#plot shaded SE and flip coordinate system
pp1<- pp1 + geom_ribbon(aes(fill=cond),alpha=.4) + 
  coord_flip() 
#plot lines and pick colors
#pp1<- pp1+geom_point(aes(y = average,color = cond),size=.75) 
pp1<- pp1+ geom_line(aes(y = average,color=cond), size=.25) 
#change axis limits
pp1<- pp1+ scale_x_continuous(limits=c(0,1),breaks = c(0,.2,.4,.6,.8,1))
pp1<- pp1+ scale_y_continuous(limits=c(0,4),
                              breaks =  c(1,2,3,4))
#change labels name, size, and font 
pp1<- pp1 + labs(y="% Signal change",
                 x= "Relative (equivolume) distance from WM")
pp1<-pp1 + theme(panel.background = element_rect(fill = "#ffffff", colour = "grey50"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 legend.title = element_blank(),
                 legend.position = "right",
                 #legend.background = element_rect(fill = "transparent", size=0.25, linetype="solid"),
                 axis.text = element_text(colour = "black", 
                                          size = 7),
                 legend.text = element_text(size=5))
pp1<-pp1 + theme(axis.title.x = element_text(size=9)) #x title
pp1<-pp1 + theme(axis.title.y = element_text(size =9)) # y title
#change line and fill colors & legend labels
pp1<-pp1+scale_colour_manual(values = c("#253494","#addd8e"),
                             labels = c("chromatic ","achromatic")
)
pp1<-pp1+scale_fill_manual(values = c("#253494","#d7da42"),labels = c("chromatic ","achromatic")
)
pp1
ggsave("depthLaminar_byEcc.pdf", device="pdf",width = 4.5, height = 3)

#plot all P and M activation across eccentricity separately 
pp2<-ggplot(PMplot[order(PMplot$depth),], 
            aes(x=depth,
                ymin = average-se,ymax = average+se,
                group=roinew
            ))
#plot shaded SE and flip coordinate system
pp2<- pp2 + geom_ribbon(aes(fill=roinew),alpha=.4) + 
  coord_flip() 
#plot lines and pick colors
#pp2<- pp2+geom_point(aes(y = average,color = roinew),size=.75) 
pp2<- pp2+ geom_line(aes(y = average,color=roinew), size=0.5) 
#change axis limits
pp2<- pp2+ scale_x_continuous(limits=c(0,1),breaks = c(0,.2,.4,.6,.8,1))
pp2<- pp2+ scale_y_continuous(limits=c(0,4),
                              breaks = c(0,1,2,3,4))
#change labels name, size, and font 
pp2<- pp2 + labs(y="% Signal change",
                 x= "Relative (equivolume) distance from WM")
pp2<-pp2 + theme(panel.background = element_rect(fill = "#ffffff", colour = "grey50"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 legend.title = element_blank(),
                 legend.position='right',
                 #legend.background = element_rect(fill = "transparent", size=0.1, linetype="solid"),
                 axis.text = element_text(colour = "black", 
                                          size = 7),
                 legend.text = element_text(size=7))
pp2<-pp2 + theme(axis.title.x = element_text(size=9)) #x title
pp2<-pp2 + theme(axis.title.y = element_text(size =9)) # y title
#pp2<-pp2 +guides(shape = guide_legend(override.aes = list(size = 0.2)),)
#change line and fill colors & legend labels
pp2<-pp2+scale_colour_manual(values = c("#081d58" ,"#3f007d", "#d7301f")
)
pp2<-pp2+scale_fill_manual(values = c("#081d58" ,"#3f007d", "#d7301f"),
                                      guide = 
                                        guide_legend(override.aes = list(size = 0.1))
)
pp2<-pp2+facet_wrap( ~ cond, nrow = 1) 
pp2
ggsave("depthLaminar_byPvsM.pdf", device="pdf",width = 3.5, height = 3)

## By: Karen Navarro 
## 08/17/2019

## This code analyses Magno-Parvo statistics from the csv file called MP_data.csv

#set environment
library(ggplot2)
library(dplyr)
library(reshape)
library(tidyr)
library(car)
library(DescTools)
setwd("/Users/navar135/Documents/UMN_Research/Projects/MagnoParvo/Stimuli")

####Set & transform Dataframes w/o Veins####
MP <- read.csv("MP_data.csv", sep = ",", header = TRUE)
#choose included datasets
included <-c("pnr161_20190124","pnr161_20190404","pnr492_20190204",
             "pnr492_20190429","pnr521_20190321","pnr739_20190131",
             "pnr739_20190415","pnr521_20190808","pnr328_20190808")
MP<- MP[MP$subj %in% included, ]
GMsig_75pct <- subset(MP,mask=="GMsig_75pct.nii")
#remove useless columns
GMsig_75pct<-GMsig_75pct[ , !(names(GMsig_75pct) %in% "mask")]
temp<- GMsig_75pct[ , !(names(GMsig_75pct) %in% "voxels")]
#transformdata from wide to long for stats
temp<- melt(temp)
#rename variables
names(temp) <- c("subject", "roi", "signal")
#filter data by eccentricity
FOV<- filter(temp,roi=="FOV_depth0" | roi=="FOV_depth2" | roi=="FOV_depth4"  | 
               roi=="FOV_depth6"  | roi=="FOV_depth8"  | roi=="FOV_depth10")
MID <- filter(temp,roi=="MID_depth0" | roi=="MID_depth2" | roi=="MID_depth4"  | 
                roi=="MID_depth6"  | roi=="MID_depth8"  | roi=="MID_depth10")
PER <-filter(temp,roi=="PER_depth0" | roi=="PER_depth2" | roi=="PER_depth4"  | 
               roi=="PER_depth6"  | roi=="PER_depth8"  | roi=="PER_depth10")
#crete a 4th column for depths & rename roi column in each datafram
sequence<- seq(0, 10, by=2)
FOV$depth<-rep(sequence, each=9)
FOV$roi <- "FOV"
MID$depth<-rep(sequence, each=9)
MID$roi <- "MID"
PER$depth<-rep(sequence, each=9)
PER$roi <- "PER"
#combine all three dataframes
GMsig_75pct <-rbind(FOV,MID,PER)
## set all variables as factors and set levels for stats
## there are two factors (roi and depth) and 3 levels
## for roi and 6 levels for depth 
GMsig_75pct$roi <- as.factor(GMsig_75pct$roi)
GMsig_75pct$depth <- as.factor(GMsig_75pct$depth)
levels(GMsig_75pct$roi)<- list("FOV"=1,"MID"=2, "PER"=3)
#repeat for data w/ veins
GMsigveins_75pct <- subset(MP, mask == "GMsigveins_75pct.nii")
GMsigveins_75pct <- GMsigveins_75pct[ , !(names(GMsigveins_75pct) %in% "mask")]
temp1<- GMsigveins_75pct[ , !(names(GMsigveins_75pct) %in% "voxels")]
#transform data from wide to long for stats
temp1<- melt(temp1)
#rename variables
names(temp1) <- c("subject", "roi", "signal")
#filter data by eccentricity
FOV1<- filter(temp1,roi=="FOV_depth0" | roi=="FOV_depth2" | roi=="FOV_depth4"  | 
                roi=="FOV_depth6"  | roi=="FOV_depth8"  | roi=="FOV_depth10")
MID1 <- filter(temp1,roi=="MID_depth0" | roi=="MID_depth2" | roi=="MID_depth4"  | 
                 roi=="MID_depth6"  | roi=="MID_depth8"  | roi=="MID_depth10")
PER1 <-filter(temp1,roi=="PER_depth0" | roi=="PER_depth2" | roi=="PER_depth4"  | 
                roi=="PER_depth6"  | roi=="PER_depth8"  | roi=="PER_depth10")
#crete a 4th column for depths & rename roi column in each datafram
FOV1$depth<-rep(sequence, each=9)
FOV1$roi <- "FOV"
MID1$depth<-rep(sequence, each=9)
MID1$roi <- "MID"
PER1$depth<-rep(sequence, each=9)
PER1$roi <- "PER"
#combine all three dataframes
GMsigveins_75pct <-rbind(FOV1,MID1,PER1)
## set all variables as factors and set levels for stats
## there are two factors (roi and depth) and 3 levels
## for roi and 6 levels for depth 
GMsigveins_75pct$roi <- as.factor(GMsigveins_75pct$roi)
GMsigveins_75pct$depth <- as.factor(GMsigveins_75pct$depth)
levels(GMsigveins_75pct$roi)<- list("FOV"=1,"MID"=2, "PER"=3)

####compute statistics####
#glance at means and sd of each level
#depth
(tapply(GMsigveins_75pct$signal,GMsigveins_75pct$depth,mean))
(tapply(GMsigveins_75pct$signal,GMsigveins_75pct$depth,sd))
#roi
(tapply(GMsigveins_75pct$signal,GMsigveins_75pct$roi,mean))
(tapply(GMsigveins_75pct$signal,GMsigveins_75pct$roi,sd))

#check for sphericity violations
# (sphericity<-ezANOVA(data=GMsigveins_75pct,
#         within=.(roi, depth),
#         wid=.(subject),
#         dv=.(signal)))
#perform 2way ANOVA with details for main effects & interaction
(twoAnova<-with(GMsigveins_75pct,aov(signal ~ roi * depth +
                                       Error(subject / (roi * depth)))))
summary(twoAnova)
#perform some posthoc pairwise comparisons against other conditions
justAnova<-aov(signal ~ roi * depth,data=GMsigveins_75pct)
(TukeyHSD(justAnova, "roi"))
(TukeyHSD(justAnova, "depth"))
#perform some posthoc pairwise comparisons against 0
#perform ttest for fovea and periphery
t.test(FOV1$signal,mu=0)
t.test(PER1$signal,mu=0)
#filter eccentricity by selected depths
FOVbyDeep<- filter(FOV1,depth==0)
FOVbyMid <- filter(FOV1, depth == 4)
FOVbySup <- filter(FOV1, depth == 10)
PERbyDeep<- filter(PER1,depth==0)
PERbyMid <- filter(PER1, depth == 4)
PERbySup <- filter(PER1, depth == 10) 
#perform ttest for each depth by ecc
t.test(FOVbyDeep$signal,mu=0)
t.test(FOVbyMid$signal,mu=0)
t.test(FOVbySup$signal,mu=0)
t.test(PERbyDeep$signal,mu=0)
t.test(PERbyMid$signal,mu=0)
t.test(PERbySup$signal,mu=0)

#ttest for each depth collapsed for ecc (not necessary)
GMsigveinsbyDeep<- filter(GMsigveins_75pct,depth==0)
GMsigveinsbyMid <- filter(GMsigveins_75pct, depth == 4)
GMsigveinsbySup <- filter(GMsigveins_75pct, depth == 10)
t.test(GMsigveinsbyDeep$signal,mu=0)
t.test(GMsigveinsbyMid$signal,mu=0)
t.test(GMsigveinsbySup$signal,mu=0)
##ttest for fovea-periphery 
Deep<-FOVbyDeep$signal-PERbyDeep$signal
t.test(Deep,mu=0)
Middle<- FOVbyMid$signal-PERbyMid$signal
t.test(Middle,mu=0)
Superficial<- FOVbySup$signal-PERbySup$signal
t.test(Superficial,mu=0)
###PLOT STUFF###
library(tidyverse)
#group data and compute mean and sd for veins mask
sumFOV<- FOV1 %>% 
  group_by(depth) %>% 
  summarise(average = mean(signal),sd=sd(signal))
sumMID<-MID1 %>% 
  group_by(depth) %>% 
  summarise(average = mean(signal),sd=sd(signal))
sumPER <- PER1 %>% 
  group_by(depth) %>% 
  summarise(average = mean(signal),sd=sd(signal))
sumSignal <- rbind(sumFOV, sumMID, sumPER)
roinew <- c("FOV","FOV","FOV","FOV","FOV","FOV",
         "MID","MID","MID","MID","MID","MID",
         "PER","PER","PER","PER","PER","PER")


#compute standard error
se <- sumSignal$sd/sqrt(length(included))
#standardize depths
sumSignal$depth<-sumSignal$depth/max(sumSignal$depth)
#put it into 1 dataframe by roi
laminar<-data.frame(sumSignal,roinew,se)

#
sumSignalV1<- GMsigveins_75pct %>%
  group_by(depth) %>% 
  summarise(average = mean(signal),sd=sd(signal))
seV1 <- sumSignalV1$sd/sqrt(length(included))
#standardize depths
sumSignalV1$depth<-c(0,2,4,6,8,10)
sumSignalV1$depth<-sumSignalV1$depth/max(sumSignalV1$depth)
#put it into 1 dataframe by depths
drawnV1<-data.frame(sumSignalV1,seV1)

####Make Plots####
library(ggplot2)
library(RColorBrewer)
library(plotly)
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
gg<- gg+ geom_line(aes(y = average,color=roinew), size=1) 
#change axis limits
gg<- gg+ scale_x_continuous(limits=c(0,1),breaks = c(0,.2,.4,.6,.8,1))
gg<- gg+ scale_y_continuous(limits=c(-.1,.4),
                            breaks = c(-.1,0,.1,.2,.3,.4))
#change labels name, size, and font 
gg<- gg + labs(y="% Signal change \n (P-M)/(P+M)",
          x= "Relative (equivolume) distance from WM")
gg<-gg + theme(panel.background = element_rect(fill = "white", colour = "grey50"),
               legend.title = element_blank(),
               axis.text = element_text(colour = "black", 
                                        size = 16),
               legend.text = element_text(size=16))
gg<-gg + theme(axis.title.x = element_text(size=20)) #x title
gg<-gg + theme(axis.title.y = element_text(size =20)) # y title
#change line and fill colors & legend labels
gg<-gg+scale_colour_manual(values = c("#045a8d" ,"#3690c0", "#a6bddb"), 
                          labels = c("Parafovea (PF)", "Middle (Mid)", "Periphery (Per)")
                          #values = c("0.3","0.6","0.9")
                          )
gg<- gg+scale_fill_manual(values = c("#045a8d" ,"#3690c0", "#a6bddb"), 
                          labels = c("Parafovea (PF)", "Middle (Mid)", "Periphery (Per)")
                          #values = c("0.3","0.6","0.9")
                          )
gg
ggsave("depthLaminar.pdf", device="pdf",width = 6, height = 8)
## plot data collapsed across eccentricity
gg1<- ggplot(drawnV1[order(drawnV1$depth),], 
             aes(x=depth,
                 ymin = average-seV1,ymax = average+seV1, 
                 lty = '          V1 \n (PF + Mid + Per)'))
#plot shaded SE and flip coordinate system
gg1<- gg1 + geom_ribbon(alpha=.4,fill="#2171b5") + 
  coord_flip() 
#plot line and points and set size
gg1<- gg1+geom_point(aes(y = average),size=2,color="#2171b5") 
gg1<- gg1+ geom_line(aes(y = average), size=1, color="#2171b5")
#change axis limits
gg1<- gg1+ scale_x_continuous(limits=c(0,1),breaks = c(0,.2,.4,.6,.8,1))
gg1<- gg1+ scale_y_continuous(limits=c(0,.3),
                            breaks = c(0,.1,.2,.3))
#change labels name, size, and font 
gg1<- gg1 + labs(y="% Signal change \n (P-M)/(P+M)",
               x= "Relative (equivolume) distance from WM")
gg1<-gg1 + theme(panel.background = element_rect(fill = "white", colour = "grey50"),
               legend.title = element_blank(),
               axis.text = element_text(colour = "black", 
                                        size = 16),
               legend.text = element_text(size=16))
gg1<-gg1 + theme(axis.title.x = element_text(size=20)) #x title
gg1<-gg1 + theme(axis.title.y = element_text(size =20)) # y title
gg1
ggsave("depthLaminarV1.pdf", device="pdf",width = 6, height = 8)
#
gg2<- ggplot() + 
  geom_path(data=sumFOV[order(sumFOV$depth),], aes(x=average, y = depth),color="darkred")+
  geom_path(data=sumMID, aes(x=average, y = depth),color="darkgreen")+
  geom_path(data=sumPER, aes(x=average, y = depth),color="darkblue")
gg2+ geom_ribbon(data=sumFOV,aes(x=average, y=depth,ymin = average-sd,ymax = average+sd, fill="salmon1"))
gg2

gg2<- gg2 + geom_ribbon(aes(fill=roinew)) + 
  coord_flip() 
gg2<- gg2+geom_point(aes(y = average,color = roinew)) 
gg2<- gg2+ geom_line(aes(y = average,color=roinew)) + 
  scale_color_manual(values=c('darkred','darkgreen','darkblue')) +
  scale_fill_manual(values=c("salmon1", "lightgreen","skyblue1"), name="fill")
gg2 + xlim(-.2, 0.6)
gg2
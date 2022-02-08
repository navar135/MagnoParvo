## This code analyses Magno-Parvo statistics from the csv file called MP_data.csv
## and creates laminar profiles for both MP and OD data
## By: Karen Navarro 
## 03/07/2020

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
install.packages("sjstats")
library(sjstats)
#clear out current environment
rm(list = ls(all.names = TRUE))
setwd("/Users/navar135/Documents/UMN_Research/Projects/MagnoParvo/Stimuli")

####Process MP data####
#Set & transform Dataframes w/o Veins
MPtrue <- read.csv("PminusMnorm_data_20200724.csv", sep = ",", header = TRUE)
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

#create a 4th column for depths & rename roi column in each datafram
sequence<- seq(0, 10, by=2)
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

####compute statistics####
#glance at means and sd of each level
#depth
(tapply(MP$signal,MP$depth,mean))
(tapply(MP$signal,MP$depth,sd))
#roi
(tapply(MP$signal,MP$roi,mean))
(tapply(MP$signal,MP$roi,sd))

#check for sphericity violations
(sphericity<-ezANOVA(data=MP,
        within=.(roi, depth),
        wid=.(subject),
        dv=.(signal)))
#perform 2way ANOVA with details for main effects & interaction
(twoAnova<-with(MP,aov(signal ~ roi * depth +
                                       Error(subject / (roi * depth)))))
summary(twoAnova)

#perform 2way ANOVA with normal error
justAnova<-aov(signal ~ roi * depth,data=MP)
summary(justAnova)
summary.lm(justAnova)
#different way to compute anova
model <- lm(signal ~ roi + depth +roi:depth,data=MP)
Anova(model,
      type = "II")
#different way to compute tukeysHSD (gives t-ratio)
(int<-lsmeans(model,
        pairwise ~ roi:depth,
        adjust="tukey"))
(tukROi<-lsmeans(model,
        pairwise ~ roi,
        ))
#perform some posthoc pairwise comparisons against other conditions
(TukeyHSD(justAnova, "roi"))
(TukeyHSD(justAnova, "depth"))
TukeyHSD(justAnova,"roi:depth")

#compute one way anova for each roi and each depth
FOV1<- filter(MP,roi=="FOV")
fovAov<- lm(signal~depth,data=FOV1)
Anova(fovAov)
(tukFov<-lsmeans(fovAov,
        pairwise ~ depth,
        adjust="tukey"))


pairwise.t.test(FOV1$signal, FOV1$depth,
                p.adjust.method = "BH")
oneway.test(signal~depth,data=FOV1)

MID1<-filter(MP,roi=="MID")
midAov<- lm(signal~depth,data=MID1)
Anova(midAov)
(tukFov<-lsmeans(midAov,
                 pairwise ~ depth,
                 adjust="tukey"))

PER1<-filter(MP,roi=="PER")
perAov<-lm(signal~depth,data=PER1)
Anova(perAov)
(tukFov<-lsmeans(perAov,
                 pairwise ~ depth,
                 adjust="tukey"))

depth0<- filter(MP,depth=="0")
depth0Aov<- lm(signal~roi,data=depth0)
Anova(depth0Aov)
(lsmeans(depth0Aov,
                 pairwise ~ roi,
                 adjust="tukey"))
mean(depth0$voxels)
sd(depth0$voxels)

depth2<- filter(MP,depth=="2")
depth2Aov<- lm(signal~roi,data=depth2)
Anova(depth2Aov)
(lsmeans(depth2Aov,
                 pairwise ~ roi,
                 adjust="tukey"))


depth4<- filter(MP,depth=="4")
depth4Aov<- lm(signal~roi,data=depth4)
Anova(depth4Aov)
(lsmeans(depth4Aov,
         pairwise ~ roi,
         adjust="tukey"))
depth6<- filter(MP,depth=="6")
depth6Aov<- lm(signal~roi,data=depth6)
Anova(depth6Aov)
(lsmeans(depth6Aov,
         pairwise ~ roi,
         adjust="tukey"))
depth8<- filter(MP,depth=="8")
depth8Aov<- lm(signal~roi,data=depth8)
Anova(depth8Aov)
(lsmeans(depth8Aov,
         pairwise ~ roi,
         adjust="tukey"))
depth10<- filter(MP,depth=="10")
depth10Aov<- lm(signal~roi,data=depth10)
Anova(depth10Aov)
(lsmeans(depth10Aov,
         pairwise ~ roi,
         adjust="tukey"))
mean(depth10$voxels)
sd(depth10$voxels)

#perform some posthoc pairwise comparisons against 0
#perform ttest for fovea and periphery
t.test(FOV$signal,mu=0)
t.test(PER$signal,mu=0)
#compare significance difference between fovea and periphery
t.test(FOV$signal,PER$signal)
#filter eccentricity by selected depths
FOVbyDeep<- filter(FOV,depth==0)
FOVbyMid <- filter(FOV, depth == 4)
FOVbySup <- filter(FOV, depth == 10)
PERbyDeep<- filter(PER,depth==0)
PERbyMid <- filter(PER, depth == 4)
PERbySup <- filter(PER, depth == 10) 
#perform ttest for each depth by ecc
t.test(FOVbyDeep$signal,mu=0)
t.test(FOVbyMid$signal,mu=0)
t.test(FOVbySup$signal,mu=0)
t.test(PERbyDeep$signal,mu=0)
t.test(PERbyMid$signal,mu=0)
t.test(PERbySup$signal,mu=0)

(onewayAov<-aov(signal ~ roi, data = MP))
summary(onewayAov)
summary.lm(onewayAov)
TukeyHSD(onewayAov)
#ttest for each depth collapsed for ecc (not necessary)
GMsigveinsbyDeep<- filter(MP,depth==0)
GMsigveinsbyMid <- filter(MP, depth == 4)
GMsigveinsbySup <- filter(MP, depth == 10)
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

#ttest for 2 deepest layers vs two superficial layers
deep<- filter(MP, depth == 0|depth ==2)
sup<- filter(MP, depth == 8|depth ==10)
t.test(deep$signal,sup$signal)
deepest<-filter(MP, depth == 0)
supMost<-filter(MP, depth == 10)
t.test(deepest$signal,supMost$signal)

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

#
sumSignalV1<- MP %>%
  group_by(depth) %>% 
  dplyr::summarise(average = mean(signal),sd=sd(signal))
seV1 <- sumSignalV1$sd/sqrt(length(MPtrue$subj))
#standardize depths
sumSignalV1$depth<-c(0,2,4,6,8,10)
sumSignalV1$depth<-sumSignalV1$depth/max(sumSignalV1$depth)
#put it into 1 dataframe by depths
drawnV1<-data.frame(sumSignalV1,seV1)

####process OD data ####
#Set & transform Dataframes w/o Veins
# OD <- read.csv("OD_data.csv", sep = ",", header = TRUE)
# 
# #remove useless columns
# ODveins<-OD[ , !(names(OD) %in% "mask")]
# ODveins<- ODveins[ , !(names(ODveins) %in% "voxels")]
# #transformdata from wide to long for stats
# ODveins<- melt(ODveins)
# #rename variables
# names(ODveins) <- c("subject","roi", "signal")
# #filter data by eccentricity
# FOV<- filter(ODveins,roi=="FOV_depth0" | roi=="FOV_depth2" | roi=="FOV_depth4"  | 
#                roi=="FOV_depth6"  | roi=="FOV_depth8"  | roi=="FOV_depth10")
# MID <- filter(ODveins,roi=="MID_depth0" | roi=="MID_depth2" | roi=="MID_depth4"  | 
#                 roi=="MID_depth6"  | roi=="MID_depth8"  | roi=="MID_depth10")
# PER <-filter(ODveins,roi=="PER_depth0" | roi=="PER_depth2" | roi=="PER_depth4"  | 
#                roi=="PER_depth6"  | roi=="PER_depth8"  | roi=="PER_depth10")
# #crete a 4th column for depths & rename roi column in each datafram
# sequence<- seq(0, 10, by=2)
# FOV$depth<-rep(sequence, each=10)
# FOV$roi <- "FOV"
# MID$depth<-rep(sequence, each=10)
# MID$roi <- "MID"
# PER$depth<-rep(sequence, each=10)
# PER$roi <- "PER"
# #combine all three dataframes
# ODveins <-rbind(FOV,MID,PER)
# 
# #
# sumV1<- ODveins %>%
#   group_by(depth) %>% 
#   dplyr::summarise(average = mean(signal),sd=sd(signal))
# seV1 <- sumV1$sd/sqrt(length(OD))
# #standardize depths
# sumV1$depth<-sumV1$depth/max(sumV1$depth)
# #put it into 1 dataframe by depths
# drawnRL_V1<-data.frame(sumV1,seV1)

#### Selected depths statistics ####
MPsub<- subset(MP, depth %in% c("0", "10"))

#check for sphericity violations
(sphericity<-ezANOVA(data=MPsub,
                     within=.(roi, depth),
                     wid=.(subject),
                     dv=.(signal)))

#perform 2way ANOVA with details for main effects & interaction with calculated error
#same as anova from sphericity
(twoAnova<-aov(signal ~ roi * depth +
                         Error(subject / (roi * depth)),data=MPsub))
summary(twoAnova)
eta_sq(twoAnova)
anova_stats(twoAnova)
#perform 2way ANOVA with normal error 
justAnova1<-aov(signal ~ roi * depth,data=MPsub) #same results as model1
eta_sq(justAnova1)
summary(justAnova1)
anova_stats(justAnova1)

summary.lm(justAnova1)
#different way to compute anova
model <- lm(signal ~ roi + depth +roi:depth,data=MPsub) #same results as model1
Anova(model,
      type = "II")
anova_stats(model)
#same as justAnova1
model1 <- lm(signal ~ roi*depth,data=MPsub) #this is what we are using in the paper! 
Anova(model1,type = "II")
anova_stats(model1)

#different way to compute tukeysHSD (gives t-ratio)
(int<-lsmeans(model,
              pairwise ~ roi:depth,
              adjust="tukey"))
(tukROi<-lsmeans(model1,
                 pairwise ~ roi,
                 adjust="tukey"))
(tukdepth<-lsmeans(model1,
                   pairwise ~ depth,
                   adjust="tukey"))
#perform a bonferroni 
(bonRoi<-lsmeans(model1,
                   pairwise ~ roi,
                   adjust="bonferroni"))
(bonDepth<-lsmeans(model1,
                 pairwise ~ depth,
                 adjust="bonferroni"))
(pairwise.t.test(MPsub$signal, MPsub$roi, p.adj = "bonferroni"))
(pairwise.t.test(MPsub$signal, MPsub$depth, p.adj = "bonferroni"))
#perform some posthoc pairwise comparisons against other conditions
(TukeyHSD(justAnova1, "roi"))
(TukeyHSD(justAnova1, "depth"))
TukeyHSD(justAnova1,"roi:depth")

#compute one way anova for each roi and each depth
FOV2<- filter(MPsub,roi=="FOV")
fovAov<- lm(signal~depth,data=FOV1)
Anova(fovAov)
(tukFov<-lsmeans(fovAov,
                 pairwise ~ depth,
                 adjust="tukey"))

pairwise.t.test(FOV1$signal, FOV1$depth,
                p.adjust.method = "BH")
oneway.test(signal~depth,data=FOV1)

MID2<-filter(MPsub,roi=="MID")
midAov<- lm(signal~depth,data=MID1)
Anova(midAov)
(tukFov<-lsmeans(midAov,
                 pairwise ~ depth,
                 adjust="tukey"))

PER2<-filter(MPsub,roi=="PER")
perAov<-lm(signal~depth,data=PER1)
Anova(perAov)
(tukFov<-lsmeans(perAov,
                 pairwise ~ depth,
                 adjust="tukey"))

depth0<- filter(MPsub,depth=="0")
depth0Aov<- lm(signal~roi,data=depth0)
Anova(depth0Aov)
(lsmeans(depth0Aov,
         pairwise ~ roi,
         adjust="tukey"))
mean(depth0$voxels)
sd(depth0$voxels)

depth10<- filter(MPsub,depth=="10")
depth10Aov<- lm(signal~roi,data=depth10)
Anova(depth10Aov)
(lsmeans(depth10Aov,
         pairwise ~ roi,
         adjust="tukey"))
mean(depth10$voxels)
sd(depth10$voxels)

#perform some posthoc pairwise comparisons against 0
#perform ttest for fovea and periphery
t.test(FOV2$signal,mu=0)
t.test(PER2$signal,mu=0)
#compare significance difference between fovea and periphery
t.test(FOV2$signal,PER$signal)

####Compute Means and SEM across sessions for bands measurements ####
dorsal<-read.csv("dorsal_spacing.csv", header = F)
names(dorsal) <- c("subject", "day1", "day2")
ventral<- read.csv("ventral_spacing.csv", header = F)
names(ventral) <- c("subject", "day1", "day2")
dorsal<- melt(dorsal)
dorsal<- na.omit(dorsal)
ventral<-melt(ventral)
ventral<- na.omit(ventral)
names(dorsal) <- c("subject", "day", "meanSpace")
names(ventral) <- c("subject", "day", "meanSpace")
dorsal$totMSpace<-dorsal$meanSpace*2 
ventral$totMSpace<-ventral$meanSpace*2 

#find mean and SEM
mean(dorsal$totMSpace)
sd(dorsal$totMSpace)/sqrt(length(dorsal$subject))
mean(ventral$totMSpace)
sd(ventral$totMSpace)/sqrt(length(ventral$subject))
#find total mean
(sum(dorsal$totMSpace)+sum(ventral$totMSpace))/(8+6)

#create a dataframe with the same labels and columns as "laminar" so we can 
#plot single subject points on the laminar profiles
MPplot<- MP
names(MPplot)<- c("subject","voxels","roinew", "average","depth") 
#standardize depths
MPplot$depth<-as.numeric(paste(MPplot$depth))
MPplot$depth<-MPplot$depth/max(MPplot$depth)
#here we are creating the sd of the average signal created 
#by each subject collapsed across roi 
#so it is averages the three means at one depth and finding the sd in that 
sdMP<- MPplot %>%
  group_by(subject,depth) %>% 
  dplyr::summarise(sd=sd(average))
#order data by depth then subject so we can duplicate the sd for each subject
#and add it to the data.frame
MPplot <- MPplot[order(MPplot$depth),] 
MPplot <- MPplot[order(MPplot$subject),] 
sdCol<- rep(sdMP$sd, each = 3)
MPplot$sd<- sdCol
#remove all middle depth 
MPplot<- dplyr::filter(MPplot, depth != 0.4 & 
                         depth != 0.6) 
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

#
sumSignalV1<- MP %>%
  group_by(depth) %>% 
  dplyr::summarise(average = mean(signal),sd=sd(signal))
seV1 <- sumSignalV1$sd/sqrt(length(MP$subj))
#standardize depths
sumSignalV1$depth<-c(0,2,4,6,8,10)
sumSignalV1$depth<-sumSignalV1$depth/max(sumSignalV1$depth)
#put it into 1 dataframe by depths
drawnV1<-data.frame(sumSignalV1,seV1)

##Plot data across eccentricity depth vs % signal change
#plot data and switch default axis for ploting SE
gg<- ggplot(laminar[order(laminar$depth),], 
            aes(x=depth,ymin = average-se,ymax = average+se, 
                group=roinew))
#plot shaded SE and flip coordinate system
gg<- gg + geom_ribbon(aes(fill=roinew),alpha=.4) + 
  coord_flip() 
#plot lines and pick colors
#gg<- gg+geom_point(aes(y = average,color = roinew),size=2) 
gg<- gg+ geom_line(aes(y = average,color=roinew), size=1) 
#change axis limits
gg<- gg+ scale_x_continuous(limits=c(0,1),breaks = c(0,.2,.4,.6,.8,1))
gg<- gg+ scale_y_continuous(limits=c(-.2,.45),
                            breaks = c(-.1, 0,.1,.2,.3,0.4))
#change labels name, size, and font 
gg<- gg + labs(y="% Signal change \n (C-A)/(C+A)",
          x= "Relative (equivolume) distance from WM")
gg<-gg + theme(panel.background = element_rect(fill = "white", colour = "grey50"),
               legend.title = element_blank(),
               axis.text = element_text(colour = "black", 
                                        size = 10),
               legend.text = element_text(size=10))
gg<-gg + theme(axis.title.x = element_text(size=12)) #x title
gg<-gg + theme(axis.title.y = element_text(size =12)) # y title
#change line and fill colors & legend labels
gg<-gg+scale_colour_manual(values = c("#014636" ,"#02818a", "#67a9cf"), 
                          labels = c("parafovea", "middle", "periphery")
                          #values = c("0.3","0.6","0.9")
                          )
gg<- gg+scale_fill_manual(values = c("#014636" ,"#02818a", "#67a9cf"), 
                          labels = c("parafovea", "middle", "periphery")
                          #values = c("0.3","0.6","0.9")
                          )

gg
ggsave("depthLaminar_reviews.pdf", device="pdf",width = 4, height = 5)

## try to plot single subjects using depth as a category
singPlot <- MPsub
names(singPlot)<- c("subject","voxels","roinew", "average","depth") 
#this will offset the grouped data so it s easier to visualization
dodge <- position_dodge(width=0.5)  

bb2<- ggplot(singPlot, aes(y= average, x = depth, colour=roinew))+ 
  geom_dotplot(binaxis='y', stackdir='center',aes(fill=roinew), dotsize = 0.5, 
               stackgroups = TRUE, binpositions = "all", position=dodge) 
bb2<- bb2+stat_summary(mapping = aes(group = roinew, color = roinew,fill=roinew), 
                       #fun.data="mean_sdl", 
                       fun.args = list(mult=1),geom="crossbar", size =1,width=0.2, alpha = 0.5,  position=dodge)
bb2<- bb2+ coord_flip()
bb2<-bb2+scale_colour_manual(values = c("#014636" ,"#02818a", "#67a9cf"), 
                             labels = c("parafovea", "middle", "periphery"))
bb2<- bb2+scale_fill_manual(values = c("#014636" ,"#02818a", "#67a9cf"), 
                            labels = c("parafovea", "middle", "periphery"))
#change axis limits
bb2<- bb2+ scale_y_continuous(limits=c(-0.3,1),
                              breaks = c(-.3,-.2,-.1,0,.1,.2,.3,.4,
                                         0.5,0.6,0.7,0.8,0.9,1.0))
bb2<- bb2+ scale_x_discrete(labels = c("deep","superficial"))
#change labels name, size, and font 
bb2<- bb2 + labs(y="% Signal change \n (C-A)/(C+A)",
                 x= "Relative distance from WM")
bb2<-bb2 + theme(panel.background = element_rect(fill = "white", colour = "grey50"),
                 legend.title = element_blank(),
                 axis.text = element_text(colour = "black", 
                                          size = 10),
                 legend.text = element_text(size=10))
bb2<-bb2 + theme(axis.title.x = element_text(size=12)) #x title
bb2<-bb2 + theme(axis.title.y = element_text(size =12)) # y title

bb2
ggsave("depthLaminar_singleSub_reviews.pdf", device="pdf",width = 6.5, height = 5)

#dev.off()
## plot data collapsed across eccentricity
gg1<- ggplot(drawnV1[order(drawnV1$depth),], 
             aes(x=depth,
                 ymin = average-seV1,ymax = average+seV1 
                 #lty = 'V1'
                 ))
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
                                        size = 10),
               legend.text = element_text(size=10))
gg1<-gg1 + theme(axis.title.x = element_text(size=12)) #x title
gg1<-gg1 + theme(axis.title.y = element_text(size =12)) # y title
gg1
ggsave("depthLaminarV1.pdf", device="pdf",width = 3, height = 5)
#dev.off()

##plot OD stuff 
## plot data collapsed across eccentricity
gg3<- ggplot(drawnRL_V1[order(drawnRL_V1$depth),], 
             aes(x=depth,
                 ymin = average-seV1,ymax = average+seV1 
                 #lty = 'V1'
                 ))
#plot shaded SE and flip coordinate system
gg3<- gg3 + geom_ribbon(alpha=.4,fill="#41ae76") + 
  coord_flip() 
#plot line and points and set size
gg3<- gg3+geom_point(aes(y = average),size=2,color="#41ae76") 
gg3<- gg3+ geom_line(aes(y = average), size=1, color="#41ae76")
#change axis limits
gg3<- gg3+ scale_x_continuous(limits=c(0,1),breaks = c(0,.2,.4,.6,.8,1))
gg3<- gg3+ scale_y_continuous(limits=c(0,.3),
                              breaks = c(0,.1,.2,.3))
#change labels name, size, and font 
gg3<- gg3 + labs(y="% Signal change \n |L-R|/(L+R)",
                 x= "Relative (equivolume) distance from WM")
gg3<-gg3 + theme(panel.background = element_rect(fill = "white", colour = "grey50"),
                 legend.title = element_blank(),
                 axis.text = element_text(colour = "black", 
                                          size = 10),
                 legend.text = element_text(size=10))
gg3<-gg3 + theme(axis.title.x = element_text(size=12)) #x title
gg3<-gg3 + theme(axis.title.y = element_text(size =12)) # y title
ggsave("depthLaminarV1_OD.pdf", device="pdf",width = 3, height = 5)
gg3
#grid.arrange(gg1, gg3, ncol=2)
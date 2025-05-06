
#title: "Beta Diversity 1969 vs. 2024 understory"
#author: "Melissa DeSiervo" "mhdesiervo@gmail.com, desiervm@union.edu"
#date: "4/7/2025"


## load packages## 

library(plyr)
library(dplyr)
library(ggplot2)
library(devtools)
library(tidyr)
library(vegan) ##must be the most recent version from github##
library(viridis)
library(tidyverse)
library(ggpubr)
library(ecotraj)
library(here)


###  read it cleaned understory data, 1969 and 2024###

forb_long <- read.csv(here("SAWYER/rawdata/4.6.25_Cleaned_Plant_List.csv"), header=T)  ## updated in Jan 2025##

## a few manual changes##

## drop plot B4##

forb_long2<-subset(forb_long, Plot!="B04")

forb_long2<-subset(forb_long2, Taxon!= "Unknown herb")



### Group elevations into 6 bins ### <1600 m, 1600-1800, 1800-2000, 2000-2200, 2200-2400, 2400+ ##

### CHANGED Group elevations into 5 bins ### <1600 m, 1600-1800, 1800-2000, 2000-2200, 2200+ ##

forb_long3 <- forb_long2 %>% mutate(Elev_bin=cut(Elevation, breaks=c(-Inf, 1600, 1800, 2000, 2200, Inf), labels=c("A", "B", "C", "D", "E")))


#### plot data ###

plotdata<-forb_long3 %>% select(Plot, Year, Elevation, Elevation, Elev_bin) %>% distinct() %>% mutate(PlotYear = paste(Plot, Year, sep='_')) %>% mutate(ElevYear=paste(Elev_bin, Year, sep='_'))%>% select(PlotYear, Plot, Year, Elevation, Elevation, Elev_bin, ElevYear)

## veg data ###

## identify doubles##

doubles<-forb_long3 %>% mutate(PlotYear = paste(Plot, Year, sep='_')) %>% dplyr::group_by(PlotYear, Taxon) %>%  dplyr::summarise(n = dplyr::n(), .groups = "drop") 

### for multiple entries (doubles etc) keep the one with the highest cover class ##

forb_long4<-forb_long3 %>% mutate(PlotYear = paste(Plot, Year, sep='_')) %>% dplyr::group_by(PlotYear, Taxon) %>%  slice_max(Class, with_ties = FALSE)

## check that there are no more doubles##

doubles_check<-forb_long4 %>% dplyr::group_by(PlotYear, Taxon) %>%  dplyr::summarise(n = dplyr::n(), .groups = "drop") 
## good..worked###

#### Convert class to cover values### ### Median cover: 1=0.2, 2=0.5, 3=5.5, 4=17.5, 5=37.5, 6=62.5, 7=87.5.##

forb_long5<-forb_long4 %>% mutate(Cover2=cut(Class, breaks=c(0,1,2,3,4,5,6,7), labels=c("0.2", "0.5", "5.5", "17.5", "37.5", "62.5", "87.5")))

forb_long5$Cover2<-as.numeric(as.character(forb_long5$Cover2))

#### Also make the OTV transformed version of data###

forb_long5<-forb_long5 %>% mutate(OTV = 1.415*log(Cover2) + 2)
  
## creates wide version of the data using raw COVER data###

forb_wide<-forb_long5 %>% mutate(PlotYear = paste(Plot, Year, sep='_')) %>% select(PlotYear, Taxon, Cover2) %>% pivot_wider(names_from = Taxon, values_from = Cover2)

# fill with 0s ###
forb_wide[is.na(forb_wide)] <- 0

dim(forb_wide)


## creates wide version of the data using OTV transformed cover data###

forb_wideOTV<-forb_long5 %>% mutate(PlotYear = paste(Plot, Year, sep='_')) %>% select(PlotYear, Taxon, OTV) %>% pivot_wider(names_from = Taxon, values_from = OTV)

# fill with 0s ###
forb_wideOTV[is.na(forb_wideOTV)] <- 0

dim(forb_wideOTV)


###calculate the alpha richness, evenness, shannon diversity by plot###

##species richness###
Richness<-specnumber(forb_wide[,2:222])##calculates species richness for all sites
ShannonD<-diversity(forb_wide[,2:222])
Evenness<- diversity(forb_wide[,2:222])/log(specnumber(forb_wide[,2:222])) 

alphaplotdata<-cbind(plotdata, Richness, ShannonD, Evenness)


alphadiversityavg<- alphaplotdata %>%  group_by(ElevYear) %>% dplyr::summarise(meanShannon=mean(ShannonD, na.rm=TRUE),sdShannon=sd(ShannonD, na.rm=TRUE),seShannon=sdShannon/sqrt(n()), n=n())

alpha_elev_df<-separate(data = alphadiversityavg, col = ElevYear, into = c("Elev_bin", "Year"))



## plot diversity##

alpha_elev_plot<-ggplot(aes(x = Elev_bin, y = meanShannon ,group=Year, color=Year),data = alpha_elev_df) +
  geom_point() +
  geom_errorbar(aes(ymin=meanShannon-seShannon, ymax=meanShannon+seShannon, width=0.5))+
  ylab("Shannon Diversity") +
  xlab("Elevation (m)")+ scale_x_discrete(labels=c('< 1600m', '1600-1800', '1800-2000', '2000-2200', '2200+'))+
  mytheme+
  theme(legend.position="right")


### BETA DIVERSITY####

##using the apply function to calculation the dispersion (beta disp) for each elevation group for each year#

beta_elev<-betadisper(vegdist(forb_wide[,2:222]),plotdata$ElevYear, type = "centroid", bias.adjust = TRUE)    ####Note that you have to specify the centroid ##  Bias adjust = true should adjust for differences in sample sizes##

### pull out the distances between each point and centroid###

distances<-as.data.frame(beta_elev$distances, row.names = rownames(plotdata$ElevYear))

distances2<-cbind(plotdata, distances)

colnames(distances2)[colnames(distances2) == 'beta_elev$distances'] <- 'dist'

## calculated mean and se around group distances###

meanSEdist<- distances2 %>%  group_by(ElevYear) %>% dplyr::summarise(meandist=mean(dist, na.rm=TRUE),sddist=sd(dist, na.rm=TRUE),sedist=sddist/sqrt(n()), n=n())
## good..this matches the output from betadisper and now we have SEs for each group###


beta_elev_df<-separate(data = meanSEdist, col = ElevYear, into = c("Elev_bin", "Year"))


### make a figure##

###<1600 m, 1600-1800, 1800-2000, 2000-2200, 2200+ ##

mytheme<- theme_bw()+ theme(axis.line.x= element_line(colour = "black", size=0.3))+theme(axis.line.y= element_line(colour = "black", size=0.3))+theme(axis.text.x=element_text(size=8, colour = "black"))+theme(axis.text.y=element_text(size=8, colour = "black"))+theme(axis.title=element_text(size=12))+theme(plot.title=element_text(size=6) +theme(plot.title = element_text(hjust = 0.5)))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(legend.position = "none")+ theme(strip.background = element_blank()) + theme(strip.text = element_text(face = "italic"))


beta_elev_plot<-ggplot(aes(x = Elev_bin, y = meandist ,group=Year, color=Year),data = beta_elev_df) +
  geom_point() +
  geom_errorbar(aes(ymin=meandist-sedist, ymax=meandist+sedist, width=0.5))+
  ylab("Beta Diversity") +
  xlab("Elevation (m)")+ scale_x_discrete(labels=c('< 1600m', '1600-1800', '1800-2000', '2000-2200', '2200+'))+
  mytheme+
  theme(legend.position="right")
  


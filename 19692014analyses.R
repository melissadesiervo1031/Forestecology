---
  title: "Cover_elevation_stats_fig"
author: "Melissa DeSiervo"
date: "4/7/2021"
output: html_document
---
  
  library(vegan)
library(indicspecies)
library(tidyverse)
library(data.table)
library(ggplot2)
library(gplots)
library(rcompanion)
library(betapart)
library(simba)
library(weights)
library(stringr)
library(viridis)



```{r setup, include=FALSE}

setwd("C:/Users/Mitzi/Desktop/Collaborations/Sugar Creek Aug 2020/'")

###import files##


##read in table with species data
sp_imp<-read.csv("C:/Users/Mitzi/Desktop/Academic Collaborations/Sugar Creek Aug 2020/russian_species.csv", header=TRUE)

##read in table with env data
env_imp<-read.csv("C:/Users/Mitzi/Desktop/Academic Collaborations/Sugar Creek Aug 2020/sawyer_topo_data_09212020.csv", header=TRUE)
#make year a factor
env_imp$Year<-as.factor(env_imp$Year)


###merge env and sp data 
all_dat<-merge(sp_imp, env_imp, by="Plot")


##read in ABCO results from ERIK###
ABCOcan<-read.csv("C:/Users/Mitzi/Desktop/Academic Collaborations/Sugar Creek Aug 2020/ABCO_CAN.csv", header=TRUE)


##read in wide plot-level BCM (downscaled climate data) ##
plotbcm_wide<-read.csv("C:/Users/Mitzi/Desktop/Academic Collaborations/Sugar Creek Aug 2020/Plot_BCM_WaterYear_Monthly_wide.csv", header=TRUE)

#covert to long##
plotbcm_long <- gather(plotbcm_wide, ID, measurement, A01:C05, factor_key=TRUE)

##then make weather variables wide###

plotbcm2 <- spread(plotbcm_long, Var, measurement)

######merge the plot bcm with the env attributes of plots#####

env_imp2<-env_imp %>% select(ID, Easting, Northing, ELEV, HLI, TPI_25)

bcm_plot<-merge(plotbcm2, env_imp2, by="ID")

##filter to the years of interest###


bcm_plot$Year

bcm_plot19592014<-subset(bcm_plot, Year > 1957)

##remove duplicates##

bcm_plot19592014_2<-unique(bcm_plot19592014)

bcm_plot19592014_2$Month <- factor(bcm_plot19592014_2$Month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))

##add column water year## (so that Oct - Sept = one precipitation season)

bcm_plot19592014_3<-bcm_plot19592014_2 %>% mutate(wyear = ifelse(Month=="Oct"|Month=="Nov"|Month=="Dec", Year, Year-1))%>%select(ID, Year, wyear, Month, ELEV, HLI, Easting, Northing, aet, cwd, pck, pet, ppt, tmn, tmx)




```


```{r species names}

###species names###
dist_dat2<-dist_dat%>%mutate(SPECIES2 = SPECIES)%>%mutate(LAYER = SPECIES)

dist_dat2$SPECIES2<-str_remove(dist_dat2$SPECIES2, "CAN1969")
dist_dat2$SPECIES2<-str_remove(dist_dat2$SPECIES2, "SAP1969")
dist_dat2$SPECIES2<-str_remove(dist_dat2$SPECIES2, "SEED1969")

dist_dat2$LAYER<-str_remove(dist_dat2$LAYER, "1969")
dist_dat2$LAYER<-str_remove(dist_dat2$LAYER, "ABCO")
dist_dat2$LAYER<-str_remove(dist_dat2$LAYER, "ABMA")
dist_dat2$LAYER<-str_remove(dist_dat2$LAYER, "ABLAS")
dist_dat2$LAYER<-str_remove(dist_dat2$LAYER, "CADE")
dist_dat2$LAYER<-str_remove(dist_dat2$LAYER, "PICBRE")
dist_dat2$LAYER<-str_remove(dist_dat2$LAYER, "PINCON")
dist_dat2$LAYER<-str_remove(dist_dat2$LAYER, "PICENG")
dist_dat2$LAYER<-str_remove(dist_dat2$LAYER, "PINLAM")
dist_dat2$LAYER<-str_remove(dist_dat2$LAYER, "PINMON")
dist_dat2$LAYER<-str_remove(dist_dat2$LAYER, "PINPON")
dist_dat2$LAYER<-str_remove(dist_dat2$LAYER, "PSME")
dist_dat2$LAYER<-str_remove(dist_dat2$LAYER, "TSMER")




dist_dat2$SPECIES2<-plyr::mapvalues(dist_dat2$SPECIES2, from=c("ABCO", "ABLAS", "ABMA", "CADE", "PICBRE", "PINCON", "PICENG", "PINLAM", "PINMON", "PINPON", "PSME", "TSMER"), to=c("Abies concolor", "Abies lasiocarpa", "Abies magnifica", "Calocedrus decurrens", "Picea breweriana", "Picea engelmannii", "Pinus contorta","Pinus lambertiana","Pinus monticola", "Pinus ponderosa", "Pseudotsuga menziesii", "Tsuga mertensiana"))
dist_dat2$LAYER<-plyr::mapvalues(dist_dat2$LAYER, from=c("CAN","SAP", "SEED"), to=c("Canopy", "Sapling", "Seedling"))





```


```{r subset by year for elevation, echo=FALSE}



###subset dataframe by species/ layer##

elev1969<-dist_dat2%>%filter(Year=="1969")%>%filter(measurement>0)
elev2015<-dist_dat2%>%filter(Year=="2015")%>%filter(measurement>0)

## split the dataframe by species / layer##

dist_dat_split<-split.data.frame(x = dist_dat,f = dist_dat$SPECIES)


dist_dat_split1969<-split.data.frame(x = elev1969,f = elev1969$SPECIES)
dist_dat_split2015<-split.data.frame(x = elev2015,f = elev2015$SPECIES)


```


```{r sum up layers for total cover, echo=FALSE}

########sum strata by species
##abco
sp_imp$ABCO<-sp_imp$ABCOCAN1969 +sp_imp$ABCOSAP1969+sp_imp$ABCOSEED1969
##ablas
sp_imp$ABLAS<-sp_imp$ABLASCAN1969 +sp_imp$ABLASSAP1969+sp_imp$ABLASSEED1969
##abma
sp_imp$ABMA<-sp_imp$ABMACAN1969+sp_imp$ABMASAP1969+sp_imp$ABMASEED1969
##cade
sp_imp$CADE<-sp_imp$CADECAN1969+sp_imp$CADESAP1969+sp_imp$CADESEED1969
##picbre
sp_imp$PIBR<-sp_imp$PICBRECAN1969+sp_imp$PICBRESAP1969+sp_imp$PICBRESEED1969
##piceng
sp_imp$PIEN<-sp_imp$PICENGCAN1969+sp_imp$PICENGSAP1969+sp_imp$PICENGSEED1969
##pincon
sp_imp$PICO<-sp_imp$PINCONCAN1969+sp_imp$PINCONSAP1969+sp_imp$PINCONSEED1969
##pinlam
sp_imp$PILA<-sp_imp$PINLAMCAN1969+sp_imp$PINLAMSAP1969+sp_imp$PINLAMSEED1969
##pinmon
sp_imp$PIMO<-sp_imp$PINMONCAN1969+sp_imp$PINMONSAP1969+sp_imp$PINMONSEED1969
##pinpon
sp_imp$PIPO<-sp_imp$PINPONCAN1969+sp_imp$PINPONSAP1969+sp_imp$PINPONSEED1969
##psme
sp_imp$PSME<-sp_imp$PSMECAN1969+sp_imp$PSMESAP1969+sp_imp$PSMESEED1969
##tsme
sp_imp$TSME<-sp_imp$TSMERCAN1969+sp_imp$TSMERSAP1969+sp_imp$TSMERSEED1969



```


```{r total canopy cover, echo=FALSE}


###merge env and sp data 
all_dat<-merge(sp_imp, env_imp, by="Plot")


try<-gather(all_dat, SPECIES,measurement, ABCO:TSME)

totalcover<-try[,38:46] #### keep in the zeroes##

avgcovertotal<-totalcover %>% group_by(Year, SPECIES) %>% dplyr::summarise(meancover=mean(measurement, na.rm=TRUE),sdcover=sd(measurement, na.rm=TRUE),secover=sd(measurement, na.rm=TRUE)/sqrt(n()))

###
totalcover1969<-subset(totalcover, Year==1969)
totalcover2015<-subset(totalcover, Year==2015)

hist(totalcover1969$measurement)

total_dist_dat_split1969<-split.data.frame(x = totalcover1969,f = totalcover1969$SPECIES)
total_dist_dat_split2015<-split.data.frame(x = totalcover2015,f = totalcover2015$SPECIES)



#####

```



```{r weighted t tests elevation by layer}


## weighted t test for elevation#####
weightedttests<-as.data.frame(stack(mapply(function(x, y) wtd.t.test(x$ELEV,y$ELEV, weight=as.numeric(x$measurement), weighty=as.numeric(y$measurement))$coefficients, dist_dat_split1969, dist_dat_split2015)))

### weighted elevation averages###
weightedmean1969<-stack(lapply(X = dist_dat_split1969,function(t){
  means1969<-weighted.mean(t$ELEV, t$measurement)
  return(means1969)
}
))

weightedmean2015<-stack(lapply(X = dist_dat_split2015,function(t){
  means2015<-weighted.mean(t$ELEV, t$measurement)
  return(means2015)
}
))

elev1969model <- lm(values ~ 1, data = weightedmean1969)
elev2015model <- lm(values ~ 1, data = weightedmean2015)

par(mfrow = c(2, 2))
plot(elev1969model)

##testing for normality##



```


```{r average canopy cover by layer}

avgcoverbylayer<-dist_dat2 %>% group_by(Year, SPECIES, SPECIES2, LAYER) %>% dplyr::summarise(meancover=mean(measurement, na.rm=TRUE),sdcover=sd(measurement, na.rm=TRUE),secover=sd(measurement, na.rm=TRUE)/sqrt(n()))

avgcoverbylayerdf<-as.data.frame(avgcoverbylayer)

#####

```


```{r Wilcoxon tests canopy cover by layer}


##wilcox signed rank tests by layer##

try<-gather(all_dat, SPECIES,measurement, ABCOCAN1969:TSMERSEED1969) ##keep the zeroes##
try$SPECIES<-as.factor(try$SPECIES)

###subset dataframe by species/ layer##

elev1969withzeros<-try%>%filter(Year=="1969")
elev2015withzeros<-try%>%filter(Year=="2015")

dist_dat_split1969_withzeros<-split.data.frame(x = elev1969withzeros,f = elev1969withzeros$SPECIES)
dist_dat_split2015_withzeros<-split.data.frame(x = elev2015withzeros,f = elev2015withzeros$SPECIES)


## wilcoxon signed rank tests for layer #####
wilcoxbylayerteststat<-as.data.frame(stack(mapply(function(x, y) wilcox.test(x$measurement,y$measurement, paired=TRUE)$statistic, dist_dat_split1969_withzeros, dist_dat_split2015_withzeros)))

wilcoxbylayertestpval<-as.data.frame(stack(mapply(function(x, y) wilcox.test(x$measurement,y$measurement, paired=TRUE)$p.value, dist_dat_split1969_withzeros, dist_dat_split2015_withzeros)))

tests<-mapply(function(x, y) wilcox.test(x$measurement,y$measurement, paired=TRUE), dist_dat_split1969_withzeros, dist_dat_split2015_withzeros)

wilcoxbylayertestpval$values<-round(wilcoxbylayertestpval$values,6)

#####we've run four tests (seedlings, saplings, canopy, and sum total), so let's use the adjustment of 0.05/4 = 0.0125 for significance (two asterisks) and 0.1/4 = 0.025 for marginal significance (one asterisk). 

wilcoxbylayertestpval <- data.table(wilcoxbylayertestpval)
wilcoxbylayertestpval$signif <- ifelse(wilcoxbylayertestpval$values < 0.0125,"significant","not significant")

```


```{r Wilcoxon tests total canopy cover}


##wilcox signed rank tests by layer##


## wilcoxon signed rank tests for layer #####
wilcoxteststat<-as.data.frame(stack(mapply(function(x, y) wilcox.test(x$measurement,y$measurement, paired=TRUE)$statistic, total_dist_dat_split1969, total_dist_dat_split2015)))
wilcoxtestpval<-as.data.frame(stack(mapply(function(x, y) wilcox.test(x$measurement,y$measurement, paired=TRUE)$p.value, total_dist_dat_split1969, total_dist_dat_split2015)))

wilcoxtestpval$values<-round(wilcoxtestpval$values,4)




```


```{r elevation violin plots by layer}

mytheme<- theme_bw()+ theme(axis.line.x= element_line(colour = "black", size=0.3))+theme(axis.line.y= element_line(colour = "black", size=0.3))+theme(axis.text.x=element_text(size=8, colour = "black"))+theme(axis.text.y=element_text(size=8, colour = "black"))+theme(axis.title=element_text(size=12))+theme(plot.title=element_text(size=6) +theme(plot.title = element_text(hjust = 0.5)))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(legend.position = "none")+ theme(strip.background = element_blank()) + theme(strip.text = element_text(face = "italic"))


dodge <- position_dodge(width = 0.8)

Elevbylayer<-ggplot(aes(x = LAYER, y = ELEV),data = dist_dat2)+   geom_violin(position=dodge, trim=TRUE, alpha=0.7, aes(fill=Year)) + geom_boxplot(width=0.2,lwd=0.3, position=dodge,outlier.color=NA, aes(color=Year))+
  facet_wrap(~SPECIES2,nrow = 4, ncol = 3) + 
  mytheme +
  ylab("Elevation (m)") +
  xlab("")+
  ylim(1000,2500)+
  theme(strip.text = element_text(size=9.5))+
  theme(panel.grid.minor = element_line(colour = "grey", size = 0.5))+ 
  theme(panel.grid.major.y  = element_line(colour = "grey", size = 0.5))+
  scale_y_continuous(limits=c(1200, 2500), breaks=seq(1200,2500,300))+
  scale_fill_manual(values=c("gray30", "green4"))+
  scale_color_manual(values=c("black", "black"))+
  theme(legend.position="right")+
  theme(legend.title = element_blank())


### to export...add in file path where you want to save ##

jpeg("C:/Users/Mitzi/Desktop/Collaborations/Sugar Creek Aug 2020/elev_violin_layer.jpeg",width = 6, height = 5,units = 'in', res = 600)
Elevbylayer
dev.off()


```



```{r HLI violin plots by layer}

mytheme<- theme_bw()+ theme(axis.line.x= element_line(colour = "black", size=0.3))+theme(axis.line.y= element_line(colour = "black", size=0.3))+theme(axis.text.x=element_text(size=8, colour = "black"))+theme(axis.text.y=element_text(size=8, colour = "black"))+theme(axis.title=element_text(size=12))+theme(plot.title=element_text(size=6) +theme(plot.title = element_text(hjust = 0.5)))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(legend.position = "none")+ theme(strip.background = element_blank()) + theme(strip.text = element_text(face = "italic"))


HLIbylayer<-ggplot(aes(x = LAYER, y = HLI,fill = Year,),data = dist_dat2) +
  geom_violin(position=position_dodge(),trim=TRUE) +
  facet_wrap(~SPECIES2,nrow = 4, ncol = 3) + 
  mytheme +
  ylab("Heat load index") +
  xlab("")+
  scale_fill_manual(values=alpha(c("gray30", "green4"), 0.7))+
  theme(strip.text = element_text(size=9.5))+
  theme(panel.grid.minor = element_line(colour = "grey", size = 0.5))+ 
  theme(panel.grid.major.y  = element_line(colour = "grey", size = 0.5))+
  theme(legend.position="right")+
  theme(legend.title = element_blank())

HLIbylayer<-ggplot(aes(x = LAYER, y =HLI),data = dist_dat2)+   geom_violin(position=dodge, trim=TRUE, alpha=0.7, aes(fill=Year)) + geom_boxplot(width=0.2,lwd=0.3, position=dodge,outlier.color=NA, aes(color=Year))+
  facet_wrap(~SPECIES2,nrow = 4, ncol = 3) + 
  mytheme +
  ylab("Heat load index") +
  xlab("")+
  theme(strip.text = element_text(size=9.5))+
  theme(panel.grid.minor = element_line(colour = "grey", size = 0.5))+ 
  theme(panel.grid.major.y  = element_line(colour = "grey", size = 0.5))+
  scale_fill_manual(values=c("gray30", "green4"))+
  scale_color_manual(values=c("black", "black"))+
  ylim(c(0.25,1))+
  theme(legend.position="right")+
  theme(legend.title = element_blank())



### to export...add in file path where you want to save ##

jpeg("C:/Users/Mitzi/Desktop/Collaborations/Sugar Creek Aug 2020/HLI_violin_layer.jpeg",width = 6, height = 5,units = 'in', res = 600)
HLIbylayer
dev.off()


```




```{r plot canopy cover by layer}

mytheme<- theme_bw()+ theme(axis.line.x= element_line(colour = "black", size=0.3))+theme(axis.line.y= element_line(colour = "black", size=0.3))+theme(axis.text.x=element_text(size=8, colour = "black"))+theme(axis.text.y=element_text(size=8, colour = "black"))+theme(axis.title=element_text(size=12))+theme(plot.title=element_text(size=6) +theme(plot.title = element_text(hjust = 0.5)))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(legend.position = "none")+ theme(strip.background = element_blank()) + theme(strip.text = element_text(face = "italic"))


##add stars for significance in paint##

##

coverbylayer<-ggplot(aes(x = LAYER, y = meancover,fill = Year,),data = avgcoverbylayerdf) +
  geom_bar(stat="identity",position=position_dodge()) +
  geom_errorbar(aes(ymin=meancover-secover, ymax=meancover+secover), width=.2,position=position_dodge(.9))+
  facet_wrap(~SPECIES2,nrow = 3, ncol = 4)+
  ylab("Cover (%)") +
  xlab("")+
  scale_fill_manual(values=alpha(c("gray30", "green4"), 1))+
  theme(strip.text = element_text(size=9.5))+
  theme(panel.grid.minor = element_line(colour = "grey", size = 0.5))+ 
  theme(panel.grid.major.y  = element_line(colour = "grey", size = 0.5))+
  scale_y_continuous(limits=c(0, 30), breaks=seq(0,30,10))+
  theme(legend.position="right")+
  theme(legend.title = element_blank())+
  mytheme


##

jpeg("C:/Users/Mitzi/Desktop/Collaborations/Sugar Creek Aug 2020/coverby_layer.jpeg",width = 6, height = 5,units = 'in', res = 600)
coverbylayer
dev.off()





```



```{r plot canopy cover total}

mytheme<- theme_bw()+ theme(axis.line.x= element_line(colour = "black", size=0.3))+theme(axis.line.y= element_line(colour = "black", size=0.3))+theme(axis.text.x=element_text(size=8, colour = "black"))+theme(axis.text.y=element_text(size=8, colour = "black"))+theme(axis.title=element_text(size=12))+theme(plot.title=element_text(size=6) +theme(plot.title = element_text(hjust = 0.5)))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(legend.position = "none")+ theme(strip.background = element_blank()) + theme(strip.text = element_text(face = "italic"))


###species names###
avgcovertotaldf2<-avgcovertotaldf%>%mutate(SPECIES2 = SPECIES)
avgcovertotaldf2$SPECIES2<-plyr::mapvalues(avgcovertotaldf2$SPECIES2, from=c("ABCO", "ABLAS", "ABMA", "CADE", "PIBR", "PICO", "PIEN", "PILA", "PIMO", "PIPO", "PSME", "TSME"), to=c("Abies concolor", "Abies lasiocarpa", "Abies magnifica", "Calocedrus decurrens", "Picea breweriana", "Picea engelmannii", "Pinus contorta","Pinus lambertiana","Pinus monticola", "Pinus ponderosa", "Pseudotsuga menziesii", "Tsuga mertensiana"))


covertotal<-ggplot(aes(x = SPECIES2, y = meancover,fill = Year,),data = avgcovertotaldf2) +
  geom_bar(stat="identity",position=position_dodge()) +
  geom_errorbar(aes(ymin=meancover-secover, ymax=meancover+secover), width=.2,position=position_dodge(.9))+
  mytheme +
  ylab("Total Cover (%)") +
  xlab("")+
  scale_fill_manual(values=alpha(c("gray30", "green4"), 1))+
  theme(strip.text = element_text(size=9.5))+
  theme(panel.grid.minor = element_line(colour = "grey", size = 0.5))+ 
  theme(panel.grid.major.y  = element_line(colour = "grey", size = 0.5))+
  scale_y_continuous(limits=c(0, 30), breaks=seq(0,30,10))+
  theme(legend.position="right")+
  theme(legend.title = element_blank())+
  theme(legend.key.size = unit(0.5, 'cm'))+
  theme(axis.text.x=element_text(angle=50, face="italic",hjust = 1))


jpeg("C:/Users/Mitzi/Desktop/Collaborations/Sugar Creek Aug 2020/totalcover.jpeg",width = 5, height = 3,units = 'in', res = 600)
covertotal
dev.off()

```




```{r plot ABCO variability bands}

mytheme<- theme_bw()+ theme(axis.line.x= element_line(colour = "black", size=0.3))+theme(axis.line.y= element_line(colour = "black", size=0.3))+theme(axis.text.x=element_text(size=8, colour = "black"))+theme(axis.text.y=element_text(size=8, colour = "black"))+theme(axis.title=element_text(size=12))+theme(plot.title=element_text(size=14) +theme(plot.title = element_text(hjust = 0.5)))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(legend.position = "none")+ theme(strip.background = element_blank()) + theme(strip.text = element_text(face = "italic"))

head(ABCOcan)

ABCOcan$Year<-as.factor(as.character(ABCOcan$Year))


ABCOcoverelev<-ggplot(aes(x = as.numeric(Elevation2015m), y = Estim,color= Year, fill = Year),data = ABCOcan) +
  geom_point() +
  geom_ribbon(aes(ymin=LowBand, ymax=HighBand), color=NA, alpha=0.2)+
  mytheme +
  ylab("Canopy Cover (%)") +
  xlab("Elevation (m)")+
  annotate(geom="text", x=2300, y=55, label="Abies concolor", size=4, fontface = "italic")+
  scale_color_manual(values=alpha(c("gray30", "green4"), 1))+
  scale_fill_manual(values=alpha(c("gray30", "green4"), 1))+
  theme(legend.position=c(0.76, 0.75))+
  theme(legend.title = element_blank())+
  theme(legend.key.size = unit(0.7, 'cm'))


jpeg("C:/Users/Mitzi/Desktop/Academic Collaborations/Sugar Creek Aug 2020/ABCOvariability.jpeg",width = 4, height = 3,units = 'in', res = 600)
ABCOcoverelev
dev.off()

```


```{r climate variables}

head(bcm_plot19592014_3)

##count months per plot per year##

countmonths<-bcm_plot19592014_3 %>% count(ID, wyear)

####

precipdata<-bcm_plot19592014_3 %>% group_by(ID, wyear, ELEV, HLI) %>% dplyr::summarise(yearlyppt=sum(ppt, na.rm=FALSE))

tmndata<-bcm_plot19592014_3 %>% group_by(ID, Year) %>% dplyr::summarise(tempmin=min(tmn, na.rm=FALSE))

growingseasondata<-bcm_plot19592014_3 %>% filter(Month=="Jun"|Month=="Jul"|Month=="Aug")%>% group_by(ID, Year) %>% dplyr::summarise(gsmax=max(tmx, na.rm=FALSE))

cwddata<-bcm_plot19592014_3 %>% group_by(ID, Year) %>% dplyr::summarise(yearlycwd=sum(cwd, na.rm=FALSE))

####

allclimate<-cbind(precipdata,tempmin=tmndata$tempmin, gsmax=growingseasondata$gsmax, cwd=cwddata$yearlycwd)

###1960s avg###


climate1960s_1<-allclimate %>% filter(wyear>1959 & wyear < 1970)%>% group_by(ID, ELEV, HLI) %>% dplyr::summarise(precip1960s=mean(yearlyppt, na.rm=FALSE), tempmin1960s=mean(tempmin, na.rm=FALSE), gsmax1960s=mean(gsmax, na.rm=FALSE), cwd1960s=mean(cwd, na.rm=FALSE), n=n())

climate1960s_2<-allclimate %>% filter(wyear>1959 & wyear < 1970)%>% group_by(ID, ELEV, HLI) %>% dplyr::summarise(precip1960s=mean(yearlyppt, na.rm=TRUE), tempmin1960s=mean(tempmin, na.rm=TRUE), gsmax1960s=mean(gsmax, na.rm=TRUE), cwd1960s=mean(cwd, na.rm=TRUE), n=n())

##2000s average##

climate2000s<-allclimate %>% filter(wyear>2004 & wyear < 2015)%>% group_by(ID, ELEV, HLI) %>% dplyr::summarise(precip2000s=mean(yearlyppt, na.rm=TRUE), tempmin2000s=mean(tempmin, na.rm=TRUE), gsmax2000s=mean(gsmax, na.rm=TRUE), cwd2000s=mean(cwd, na.rm=TRUE), n=n())

climate2000s_2<-as.data.frame(climate2000s[,4:7])

climate1960s_2<-as.data.frame(climate1960s)
####

climate1960s2000s_2<-cbind(climate1960s_2, climate2000s_2)


###
climate1960s2000s_2$diffprecip=climate1960s2000s_2$precip2000s-climate1960s2000s_2$precip1960s

climate1960s2000s_2$difftempmin=climate1960s2000s_2$tempmin2000s-climate1960s2000s_2$tempmin1960s

climate1960s2000s_2$diffgsmax=climate1960s2000s_2$gsmax2000s-climate1960s2000s_2$gsmax1960s

climate1960s2000s_2$diffcwd=climate1960s2000s_2$cwd2000s-climate1960s2000s_2$cwd1960s


##relative diff##

climate1960s2000s_2$reldiffprecip=(climate1960s2000s_2$diffprecip/climate1960s2000s_2$precip1960s)*100

climate1960s2000s_2$reldifftempmin=(climate1960s2000s_2$difftempmin/climate1960s2000s_2$tempmin1960s)*100

climate1960s2000s_2$reldiffgsmax=(climate1960s2000s_2$diffgsmax/climate1960s2000s_2$gsmax1960s*100)

climate1960s2000s_2$reldiffcwd=(climate1960s2000s_2$diffcwd/climate1960s2000s_2$cwd1960s)*100

###is regression significant?##

modelprecip<-lm(diffprecip~ELEV+HLI, data=climate1960s2000s_2)
modelGStemp<-lm(diffgsmax~ELEV+HLI, data=climate1960s2000s_2)
modelmintemp<-lm(difftempmin~ELEV+HLI, data=climate1960s2000s_2)
modelCWD<-lm(diffcwd~ELEV+HLI, data=climate1960s2000s_2)

####

mytheme<- theme_bw()+ theme(axis.line.x= element_line(colour = "black", size=0.3))+theme(axis.line.y= element_line(colour = "black", size=0.3))+theme(axis.text.x=element_text(size=10, colour = "black"))+theme(axis.text.y=element_text(size=10, colour = "black"))+theme(axis.title=element_text(size=12))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



diffprecipplot<-ggplot(aes(x = ELEV, y = diffprecip, color=HLI),data = climate1960s2000s_2) +
  geom_point() +
  mytheme +
  geom_smooth(method='lm', formula= y~x, color="black")+
  ylim(-150, 150)+
  xlim(1200, 2600)+
  ylab(expression(atop("Change in precipitation (mm)", paste("2000s - 1960s"))))+
  xlab("Elevation (m)")+
  scale_colour_viridis()




diffwinterminplot<-ggplot(aes(x = ELEV, y = difftempmin, color=HLI),data = climate1960s2000s_2) +
  geom_point() +
  mytheme +
  geom_smooth(method='lm', formula= y~x, color="black")+
  xlim(1200, 2600)+
  ylim(1.0, 2.6)+
  ylab(expression(atop("Change in min temp(C)", paste("2000s - 1960s"))))+
  xlab("Elevation (m)")+
  scale_colour_viridis()


diffgsmaxplot<-ggplot(aes(x = ELEV, y = diffgsmax, color=HLI),data = climate1960s2000s_2) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x, color="black")+
  mytheme +
  xlim(1200, 2600)+
  ylim(-0.5, 1.7)+
  ylab(expression(atop("Change in max GS temp (C)", paste("2000s - 1960s"))))+
  xlab("Elevation (m)")+
  scale_colour_viridis()

diffCWDplot<-ggplot(aes(x = ELEV, y = diffcwd, color=HLI),data = climate1960s2000s_2) +
  geom_point() +
  xlim(1200, 2600)+
  ylim(-100, 10)+
  geom_smooth(method='lm', formula= y~x, color="black")+
  mytheme +
  ylab(expression(atop("Change in CWD (mm H2O)", paste("2000s - 1960s"))))+
  xlab("Elevation (m)")+
  scale_colour_viridis()

##combine plots###

multiplot(diffprecipplot, diffwinterminplot, diffgsmaxplot,  diffCWDplot, cols=2)

jpeg("C:/Users/Mitzi/Desktop/Academic Collaborations/Sugar Creek Aug 2020/climatediffraw.jpeg",width = 8, height = 6,units = 'in', res = 600)
multiplot(diffprecipplot, diffwinterminplot, diffgsmaxplot,  diffCWDplot, cols=2)
dev.off()


####


```


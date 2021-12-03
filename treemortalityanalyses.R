##upload libraries##
library(lme4)
library(reshape2)
library(ggplot2)
library(here)
library(dplyr)
library(rms)
library(ROCR)


##aic function###
aictable<-function(X,m){
  #
  #  This function will create a full model selection table based on AICc.
  #  
  #  Inputs:
  #
  #  X is the AIC table produced by R by using AIC(model1,model2, ...)
  #  m is the sample size
  #
  #
  rnames<-row.names(X)
  AICc<-X$AIC+2*X$df*(X$df+1)/(m-X$df-1)     #small-sample correction
  logL<-X$df-X$AIC/2                         #Log-likelihood
  tab<-data.frame(X[,1],logL,AICc)           #Remove AIC column; add logL and AICc
  colnames(tab)[1]<-c("Params")              #Rename "df" column   
  row.names(tab)<-rnames
  tab<-tab[order(tab$AICc),]                 #Sort by ascending AICc value
  deltaAICc<-tab$AICc-min(tab$AICc)          #Delta AICc
  weight<-exp(-deltaAICc/2)/sum(exp(-deltaAICc/2))  #Weights
  cumwt<-weight                              #Column for cumulative weight
  for(i in 2:dim(X)[1]){                  
    cumwt[i]<-cumwt[i-1]+cumwt[i]              #Accumulate weight from the top
  }
  tab<-data.frame(tab,deltaAICc,weight,cumwt)
  tab<-round(tab,4)
  tab
}


##upload Klamath tree data ##
trees<- read.csv(here("fulltreedata.csv"), header=T, check.names=FALSE) 


##Mortality analyses done for three species: ABMA, PINCON, ABLAS### This code only shows the example for ABMA##

ABMAonly<-subset(trees, Species == "ABMA")

str(ABMAonly)

##first scale all the quantitative variables###

##scaleallvariables##
ABMAonly$DBH<-scale(ABMAonly$DBH)
ABMAonly$Elevation<-scale(ABMAonly$Elevation)
ABMAonly$TRMI<-scale(ABMAonly$TRMI)
ABMAonly$BAABMA<-scale(ABMAonly$BAABMA)
ABMAonly$BAotherABMA<-scale(ABMAonly$BAotherABMA)
ABMAonly$QMDABMA<-scale(ABMAonly$QMDABMA)
ABMAonly$Change5yrtempmin<-scale(ABMAonly$Change5yrtempmin)
ABMAonly$Change5yrGstempmax<-scale(ABMAonly$Change5yrGstempmax)
ABMAonly$Change10yrtempmin<-scale(ABMAonly$Change10yrtempmin)
ABMAonly$Change10yrGstempmax<-scale(ABMAonly$Change10yrGstempmax)
ABMAonly$raw5yrtempmin<-scale(ABMAonly$raw5yrtempmin)
ABMAonly$raw5yrppt<-scale(ABMAonly$raw5yrppt)
ABMAonly$raw5yrGstempmax<-scale(ABMAonly$raw5yrGstempmax)
ABMAonly$raw5yrGSCWD<-scale(ABMAonly$raw5yrGSCWDraw)

#quick way to visualize how all covariates in the model will vary by live and dead trees## 

tmp<-melt(ABMAonly[,c("Status4","BB","Mtoe","DBH","Elevation","TRMI","BAABMA","QMDABMA", "BAotherABMA")],
          id.vars="Status4")
ggplot(tmp,aes(factor(Status4),y=value,fill=factor(Status4)))+
  geom_boxplot()+
  facet_wrap(~variable,scales="free_y")+theme_classic()


##BEFORE CREATING MODEL, WE WANT TO LOOK AT VARIANCE INFLATION FACTOR OF QUANTITATIVE VARIABLES###

###FUNCTIONS FOR VIF CORRELATIONS#####

#corvif#
##VIF FUNCTION.
#To use:  corvif(YourDataFile)
corvif <- function(dataz) {
  dataz <- as.data.frame(dataz)
  #correlation part
  #cat("Correlations of the variables\n\n")
  #tmp_cor <- cor(dataz,use="complete.obs")
  #print(tmp_cor)
  
  #vif part
  form    <- formula(paste("fooy ~ ",paste(strsplit(names(dataz)," "),collapse=" + ")))
  dataz   <- data.frame(fooy=1,dataz)
  lm_mod  <- lm(form,dataz)
  
  cat("\n\nVariance inflation factors\n\n")
  print(myvif(lm_mod))
}


#Support function for corvif. Will not be called by the user
myvif <- function(mod) {
  v <- vcov(mod)
  assign <- attributes(model.matrix(mod))$assign
  if (names(coefficients(mod)[1]) == "(Intercept)") {
    v <- v[-1, -1]
    assign <- assign[-1]
  } else warning("No intercept: vifs may not be sensible.")
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2) stop("The model contains fewer than 2 terms")
  if (length(assign) > dim(v)[1] ) {
    diag(tmp_cor)<-0
    if (any(tmp_cor==1.0)){
      return("Sample size is too small, 100% collinearity is present")
    } else {
      return("Sample size is too small")
    }
  }
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1)) {
    result <- data.frame(GVIF=result[, 1])
  } else {
    result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  }
  invisible(result)
}
#END VIF FUNCTIONS


# calculate the vif of all quantitative covariates
ABMAonly$DBH<-scale(ABMAonly$DBH)
ABMAonly$Elevation<-scale(ABMAonly$Elevation)
ABMAonly$TRMI<-scale(ABMAonly$TRMI)
ABMAonly$BAABMA<-scale(ABMAonly$BAABMA)
ABMAonly$BAotherABMA<-scale(ABMAonly$BAotherABMA)
ABMAonly$Change5yrtempmin<-scale(ABMAonly$Change5yrtempmin)
ABMAonly$Change5yrGstempmax<-scale(ABMAonly$Change5yrGstempmax)


ppnew<-ABMAonly[,c("DBH","Elevation","TRMI","BAABMA","BAotherABMA",
                   "Change5yrtempmin","Change5yrGstempmax")]
ppnew<-as.data.frame(ppnew)
corvif(ppnew)

####Nothing too troubling about the VIF of covariates##

##MAKING THE MODELS###

##FIRSTAICnoclimate##
fullABMA <- glmer(Status4~DBH + BB + Mtoe + Elevation + TRMI + BAABMA + BAotherABMA + QMDABMA+ BB*DBH + Mtoe*DBH + BB*BAABMA + Mtoe*BAABMA + (1|Plot), family=binomial(link="logit"), data=ABMAonly, control=glmerControl(optimizer="bobyqa"))
m11<- glmer(Status4~DBH + BB + Mtoe + Elevation + BAABMA + QMDABMA +  BB*DBH + Mtoe*DBH + BB*BAABMA + Mtoe*BAABMA  + (1|Plot), family=binomial(link="logit"), data=ABMAonly, control=glmerControl(optimizer="bobyqa"))
m33<- glmer(Status4~DBH + BB + Mtoe + Elevation + BAABMA + BB*BAABMA + Mtoe*BAABMA + (1|Plot), family=binomial(link="logit"), data=ABMAonly, control=glmerControl(optimizer="bobyqa"))
m44<- glmer(Status4~DBH + BB + Mtoe + Elevation + BAABMA +  BB*DBH + Mtoe*DBH + BB*BAABMA + Mtoe*BAABMA + (1|Plot), family=binomial(link="logit"), data=ABMAonly, control=glmerControl(optimizer="bobyqa"))
m1<- glmer(Status4~DBH + BB + Mtoe + Elevation + TRMI + QMDABMA + BAABMA + BAotherABMA + (1|Plot), family=binomial(link="logit"), data=ABMAonly, control=glmerControl(optimizer="bobyqa"))
m2<- glmer(Status4~DBH + BB + Mtoe + BAABMA + BAotherABMA + BB*DBH + Mtoe*DBH + BB*BAABMA + Mtoe*BAABMA + (1|Plot), family=binomial(link="logit"), data=ABMAonly, control=glmerControl(optimizer="bobyqa"))
m3<- glmer(Status4~DBH + BB + Mtoe + BAABMA + QMDABMA +  BB*DBH + Mtoe*DBH + BB*BAABMA + Mtoe*BAABMA + (1|Plot), family=binomial(link="logit"), data=ABMAonly, control=glmerControl(optimizer="bobyqa"))
m4<- glmer(Status4~DBH + BB + Mtoe + BAABMA + QMDABMA + BB*DBH + Mtoe*DBH + (1|Plot), family=binomial(link="logit"), data=ABMAonly, control=glmerControl(optimizer="bobyqa"))
m5<- glmer(Status4~DBH + BB + Mtoe + BAABMA + QMDABMA + BB*BAABMA + Mtoe*BAABMA + (1|Plot), family=binomial(link="logit"), data=ABMAonly, control=glmerControl(optimizer="bobyqa"))
m6<- glmer(Status4~BB + Mtoe + BAABMA + QMDABMA + BB*DBH + Mtoe*DBH + BB*BAABMA + Mtoe*BAABMA + (1|Plot), family=binomial(link="logit"), data=ABMAonly, control=glmerControl(optimizer="bobyqa"))
m7<- glmer(Status4~BB + Mtoe + BAABMA + QMDABMA + (1|Plot), family=binomial(link="logit"), data=ABMAonly, control=glmerControl(optimizer="bobyqa"))
m8<- glmer(Status4~BB + Mtoe + BAABMA + QMDABMA + Elevation + TRMI+ (1|Plot), family=binomial(link="logit"), data=ABMAonly, control=glmerControl(optimizer="bobyqa"))
m9<- glmer(Status4~DBH + BB + BAABMA + QMDABMA + BB*DBH + BB*BAABMA + (1|Plot), family=binomial(link="logit"), data=ABMAonly, control=glmerControl(optimizer="bobyqa"))
nullABMA <- glmer(Status4~(1|Plot), family=binomial(link="logit"), data=ABMAonly, control=glmerControl(optimizer="bobyqa"))

rawaic<-AIC(fullABMA, m11, m33, m44, m1, m2, m3, m4, m5, m6, m7, m8, m9, nullABMA)
nR<-dim(ABMAonly)[1]  #Sample size 
aictable(rawaic,nR)

#mm11bestmodel#

##including climatic variables####
##AICallclimate##
m111 <- glmer(Status4~DBH + BB +  Mtoe + QMDABMA + BAABMA + BB*BAABMA +  Mtoe*BAABMA + raw5yrtotalCWD + (1|Plot), family=binomial(link="logit"), data=ABMAonly, control=glmerControl(optimizer="bobyqa"))
m222 <- glmer(Status4~DBH + BB +  Mtoe + QMDABMA + BAABMA + BB*BAABMA +  Mtoe*BAABMA + raw5yrtotalCWD +  Change5yrtempmin + Change5yrGstempmax + (1|Plot), family=binomial(link="logit"), data=ABMAonly, control=glmerControl(optimizer="bobyqa"))
m333 <- glmer(Status4~DBH + BB +  Mtoe + QMDABMA + BAABMA + BB*BAABMA +  Mtoe*BAABMA + raw5yrtotalCWD + Change5yrGstempmax + (1|Plot), family=binomial(link="logit"), data=ABMAonly, control=glmerControl(optimizer="bobyqa"))
m444 <- glmer(Status4~DBH + BB +  Mtoe + QMDABMA + BAABMA + BB*BAABMA +  Mtoe*BAABMA + raw5yrtotalCWD + Change5yrtempmin + (1|Plot), family=binomial(link="logit"), data=ABMAonly, control=glmerControl(optimizer="bobyqa"))
m555 <- glmer(Status4~DBH + BB + Mtoe + QMDABMA + BAABMA + BB*BAABMA +  Mtoe*BAABMA + raw5yrtotalCWD +  Change10yrtempmin + Change10yrGstempmax + (1|Plot), family=binomial(link="logit"), data=ABMAonly, control=glmerControl(optimizer="bobyqa"))
m666 <- glmer(Status4~DBH + BB +  Mtoe + QMDABMA + BAABMA + BB*BAABMA +  Mtoe*BAABMA + raw5yrtotalCWD + Change10yrGstempmax + (1|Plot), family=binomial(link="logit"), data=ABMAonly, control=glmerControl(optimizer="bobyqa"))
m777 <- glmer(Status4~DBH + BB +  Mtoe + QMDABMA + BAABMA + BB*BAABMA +  Mtoe*BAABMA + raw5yrtotalCWD + Change10yrtempmin + (1|Plot), family=binomial(link="logit"), data=ABMAonly, control=glmerControl(optimizer="bobyqa"))
m888 <- glmer(Status4~DBH + BB +  Mtoe + QMDABMA + BAABMA + BB*BAABMA +  Mtoe*BAABMA +  Change5yrtempmin + Change5yrGstempmax + (1|Plot), family=binomial(link="logit"), data=ABMAonly, control=glmerControl(optimizer="bobyqa"))
m999 <- glmer(Status4~DBH + BB +  Mtoe + QMDABMA + BAABMA + BB*BAABMA +  Mtoe*BAABMA + Change5yrGstempmax + (1|Plot), family=binomial(link="logit"), data=ABMAonly, control=glmerControl(optimizer="bobyqa"))
m1000 <- glmer(Status4~DBH + BB +  Mtoe + QMDABMA + BAABMA + BB*BAABMA +  Mtoe*BAABMA + Change5yrtempmin + (1|Plot), family=binomial(link="logit"), data=ABMAonly, control=glmerControl(optimizer="bobyqa"))
m1100 <- glmer(Status4~DBH + BB +  Mtoe + QMDABMA + BAABMA + BB*BAABMA +  Mtoe*BAABMA +  Change10yrtempmin + Change10yrGstempmax + (1|Plot), family=binomial(link="logit"), data=ABMAonly, control=glmerControl(optimizer="bobyqa"))
m1200 <- glmer(Status4~DBH + BB +  Mtoe + QMDABMA + BAABMA + BB*BAABMA +  Mtoe*BAABMA + Change10yrGstempmax + (1|Plot), family=binomial(link="logit"), data=ABMAonly, control=glmerControl(optimizer="bobyqa"))
m1300 <- glmer(Status4~DBH + BB +  Mtoe + QMDABMA + BAABMA + BB*BAABMA +  Mtoe*BAABMA + Change10yrtempmin + (1|Plot), family=binomial(link="logit"), data=ABMAonly, control=glmerControl(optimizer="bobyqa"))
nullABMA3<- glmer(Status4~DBH + BB +  Mtoe + Elevation + BAABMA + BB*BAABMA + DMR2*BAABMA + (1|Plot), family=binomial(link="logit"), data=ABMAonly, control=glmerControl(optimizer="bobyqa"))

rawaic<-AIC(m111, m222, m333, m444, m555, m666, m777, m888, m999, m1000, m1100, m1200, m1300, nullABMA2)
nR<-dim(ABMAonly)[1]  #Sample size 
aictable(rawaic,nR)

#m444 is the best model###


## a way to graph the interaction effects##

# graph interactions
ABMAonly$BB<-as.factor(ABMAonly$BB)
m1<-ggplot(ABMAonly,aes(x=BAABMA,y=Status4,colour=BB))+geom_point()+geom_smooth(method="lm",alpha=0.3)+theme_classic()
m1

ABMAonly$Mtoe<-as.factor(ABMAonly$Mtoe)
m2<-ggplot(ABMAonly,aes(x=BAABMA,y=Status4,colour=Mtoe))+geom_point()+geom_smooth(method="lm",alpha=0.3)+theme_classic()
m2

###Final model in Table 2 ###

modelABMAfinal <- glmer(Status4~DBH + BB +Elevation + QMDABMA+  Change10yrtempmin+ (1|Plot), family=binomial(link="logit"), data=ABMAonly, control=glmerControl(optimizer="bobyqa"))


summary(modelABMAfinal) ##these estimates approximate what is in table 2##

##some other things to report with this model, calculate the AUC##

#Error rate (Gelman & Hill, p 99)
ABMAonly2<-ABMAonly %>% select(Plot, Elevation, Species, Status, Status4, DBH,  BAABMA, BAotherABMA, QMDABMA, Mtoe, BB, Change10yrtempmin)

ABMAonly3<- ABMAonly2[complete.cases(ABMAonly2), ]

y.pred = predict(modelABMAfinal, type ="response") 

ABMAonly4<-cbind(ABMAonly3, y.pred)

error.rate <- mean ((ABMAonly4$y.pred>0.5 & ABMAonly4$Status4==0 |  ABMAonly4$y.pred<0.5 & ABMAonly4$Status4==1))

error.rate ##tell you how often your model predicts the wrong outcome###

vp <- val.prob(ABMAonly4$y.pred, ABMAonly4$Status4, logit)

pred <- prediction(ABMAonly4$y.pred, ABMAonly4$Status4 )

perf <- performance(pred,"tpr","fpr")
perf <- performance(pred,"tpr","fpr")
plot(perf, colorize=T)

# calculating AUC
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
auc

r.auc <- round(auc, digits = 4)
auct <- paste(c("AUC = "),r.auc,sep="")
legend(0.25, 0.25, auct ,border="white",cex=1,box.col = "white")


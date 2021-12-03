##upload libraries##
library(lme4)
library(reshape2)
library(ggplot2)
library(here)

##upload Klamath tree data ##
trees<- read.csv(here("fulltreedata.csv"), header=T, check.names=FALSE) 


##Mortality analyses done for three species: ABMA, PINCON, ABLAS###

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
ABMAonly$'5yrSpringsnowpack'<-scale(ABMAonly$'5yrSpringsnowpack')
ABMAonly$'5yrtempmin'<-scale(ABMAonly$'5yrtempmin')
ABMAonly$'5yrppt'<-scale(ABMAonly$'5yrppt')
ABMAonly$'5yrGstempmax'<-scale(ABMAonly$'5yrGstempmax')
ABMAonly$'5yrGSCWD'<-scale(ABMAonly$'5yrGSCWD')
ABMAonly$'5yrtotalCWD'<-scale(ABMAonly$'5yrtotalCWD')

#quick way to visualize how all covariates in the model will vary by live and dead trees## 

tmp<-melt(ABMAonly[,c("Status4","BB","Mtoe","DMR2","DBH","Elevation","TRMI","BAABMA","BAotherABMA","Change5yrtempmin","Change5yrGstempmax")],
          id.vars="Status4")

ggplot(tmp,aes(factor(Status4),y=value,fill=factor(Status4)))+
  geom_boxplot()+
  facet_wrap(~variable,scales="free_y")+theme_classic()


##BEFORE CREATING MODEL, WE WANT TO LOOK AT VARIANCE INFLATION FACTOR OF VARIABLES###

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


# models
null<-glmer(Status4~DBH + BB + Mtoe + DMR2 + Elevation + TRMI + BAABMA + BAother + Change5yrtempmin + Change5yrGstempmax +(1|Plot), family=binomial(link="logit"), data=ABMAonly)
null2<-glmer(Status4~DBH + Elevation + TRMI + BAABMA + BAother + Change5yrtempmin + Change5yrGstempmax +(1|Plot), family=binomial(link="logit"), data=ABMAonly)


interactions<-glmer(Status4~DBH + BB + Mtoe + BB*Mtoe+ BAABMA + BAother + BB*BAABMA+ Mtoe*BAABMA + DBH*BB+ DBH*Mtoe+ Change5yrtempmin + Change5yrGstempmax +(1|Plot),family=binomial(link="logit"),data=ABMAonly) 

ABMAonly$BB<-as.factor(ABMAonly$BB)
m1<-ggplot(ABMAonly,aes(x=BAABMA,y=Status4,colour=BB))+geom_point()+geom_smooth(method="lm",alpha=0.3)+theme_classic()
m1

ABMAonly$Mtoe<-as.factor(ABMAonly$Mtoe)
m2<-ggplot(ABMAonly,aes(x=BAABMA,y=Status4,colour=Mtoe))+geom_point()+geom_smooth(method="lm",alpha=0.3)+theme_classic()
m2

#plotlevelanalysis?##
PlotABMA<- drop_read_csv("/Red fir paper/Rfiles/PlotABMA_12_9_15.csv")


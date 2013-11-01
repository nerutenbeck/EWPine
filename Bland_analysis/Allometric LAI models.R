library(RODBC);library(nlme);library(car);library(plyr);library(doBy)
temp=odbcConnectAccess("C:/Users/Frank/Desktop/Demeritt WP- Bland.mdb")
LAI=sqlFetch(temp,"Plot_litterfall LAI_all years");Plots=sqlFetch(temp,"Plots")
SBA11=sqlFetch(temp,"2011 Sapwood Area Final");SBA08=sqlFetch(temp,"Sapwood Basal Area 2008");SBA01=sqlFetch(temp,"SapwoodCores_2001_adj4RadialGrowth")
close(temp)
# Renames Sapwood areas from all years to "SBA"
names(SBA08)[names(SBA08)=="SAP_BA_08"] = "SBA";names(SBA01)[names(SBA01)=="MeanSBA"] = "SBA";names(SBA11)[names(SBA11)=="SapAreaBH"] = "SBA"
Plots=subset(Plots,select=c(Stand,PlotID,Plot,AreaExpansion)) # Subsets relevant plot attribute information
DM=read.csv("C:/Users/Frank/Desktop/R Output/DM.csv")
DM=subset(DM,select=c(Year,PlotID,Plot,TOPHT,SI,SDI))
############################################################################################
#########################  Litterfall LAI   ##################################
# Fits a spline to LAI and adjusts year to appropriate LAI
# This removes Blue-Spring b/c there is only 1 measurement- it is added back later
LAI.BS=subset(LAI,Plot=="Blue-Spring");LAI.BS$LAI.spline=LAI.BS$Litterfall_LAI
LAI=subset(LAI,Plot!="Blue-Spring")
LAImod <- function(df) lm(Litterfall_LAI~bs(Year),df,na.action=na.omit)
#this dlply will fit the above function to each Stand, PlotID, and Plot in the dataframe and return it as a list
LAIList <- dlply(LAI, .(Stand,PlotID,Plot), LAImod)
#this dlply will separate the data by each Stand, PlotID, and Plot in the dataframe and return it as a list
LAIData <- dlply(LAI, .(Stand,PlotID,Plot))
#(make sure they are the same length so you know you did it right)
preds.LAI <- mdply(cbind(mod = LAIList, df = LAIData), function(mod, df) {
  mutate(df, LAI.spline = predict(mod, df))# mutate adds a new 'pred' column with predictions to the old dataframe called 'preds.LAI'
})
# Replaces original LAI dataframe with one above with LAI.spline values 
LAI<-preds.LAI;LAI$X1<-NULL # Removes extra column created when fitting spline
# Merges Blue-Spring plot 
LAI=merge(LAI,LAI.BS,all=T)
################################################################################
# This uses litterfall collected 2 years after inventory year 
ifelse(LAI$Year<2010,LAI$Year-2,LAI$Year)
head(LAI)
##############################################################################################
###############################################################################################
# Loads Tree summary statistics- Metric units
TSS=read.csv("C:/Users/Frank/Desktop/R output/TSS.csv")
TSS=subset(TSS,select=c(Year,Plot,TreeNo,Stand,PlotID,thinType,DBH,TotHt,TPH,SI,CL_LLB,mLCR,TreeBAperHa,TreeBAperAc))
# If mLCR is larger than 1 then mLCR=1
TSS$mLCR<-ifelse(TSS$mLCR>1,1,TSS$mLCR)
names(TSS)[names(TSS)=="thinType"] = "Treatment";names(TSS)[names(TSS)=="CL_LLB"] = "CL"

PSS=ddply(TSS,.(Year,Stand,PlotID,Plot,Treatment),colwise(sum,c('DBH','TotHt','TPH','CL','mLCR','TreeBAperHa','TreeBAperAc')),.progress="text")
names(PSS)[names(PSS)=="DBH"] = "Plot.DBH";names(PSS)[names(PSS)=="TotHt"] = "Plot.Ht";names(PSS)[names(PSS)=="DBH"] = "Plot.DBH";
names(PSS)[names(PSS)=="TreeBAperHa"]= "BAperHa";names(PSS)[names(PSS)=="TreeBAperAc"]= "BAperAc";names(PSS)[names(PSS)=="CL"]= "Plot.CL"
names(PSS)[names(PSS)=="mLCR"]= "Plot.mLCR"
head(PSS)
# Merge Plot attribute data with PSS- This removes 2x2 out etc.
PSS=merge(PSS,Plots,by=c("Stand", "PlotID","Plot"))
# converts summed DBH, Ht, CL to per Hectare basis
PSS$DBHperHa=PSS$Plot.DBH*PSS$AreaExpansion;PSS$HtperHa=PSS$Plot.Ht*PSS$AreaExpansion
PSS$CLperHa=PSS$Plot.CL*PSS$AreaExpansion;PSS$mLCRperHa=PSS$Plot.mLCR*PSS$AreaExpansion
PSS=subset(PSS,select=c(Year,Stand,Plot,PlotID,Treatment,TPH,BAperHa,BAperAc,AreaExpansion,DBHperHa,HtperHa,CLperHa,mLCRperHa))
# Calculates stand averages
PSS$MeanCL=PSS$CLperHa/PSS$AreaExpansion
PSS$MeanmLCR=PSS$mLCRperHa/PSS$TPH
PSS$MeanDBH=PSS$DBHperHa/PSS$AreaExpansion
PSS$MeanBA=PSS$BAperHa/PSS$AreaExpansion
###############################################################
# Merges Litterfall LAI with per hectare attribute information
PSS=merge(LAI,PSS,by=c("Year","Stand","Plot","PlotID","Treatment"))
################################################
# Takes sapwood area data, sums to stand level, combines them and merges with previous data.  
SBA11=subset(SBA11,select=c(Year,Plot,PlotID,TreeNo,SBA))
SBA11=ddply(SBA11,.(Year,PlotID,Plot),colwise(sum,c('SBA')),.progress="text")
SBA08=subset(SBA08,select=c(Year,Plot,PlotID,TreeNo,SBA))
SBA08=ddply(SBA08,.(Year,PlotID,Plot),colwise(sum,c('SBA')),.progress="text");SBA08=subset(SBA08,PlotID!="NA")
SBA01=subset(SBA01,select=c(Year,Plot,PlotID,TreeNo,SBA))
SBA01=ddply(SBA01,.(Year,PlotID,Plot),colwise(sum,c('SBA')),.progress="text");SBA01=subset(SBA01,PlotID!="NA")
SBA1=rbind(SBA01,SBA08)
SBA=rbind(SBA1,SBA11)
SBA=merge(SBA,Plots,by=c("PlotID","Plot"))
SBA$SBAperHa=SBA$SBA*SBA$AreaExpansion
SBA=subset(SBA,select=c(Year,Stand,Plot,PlotID,SBAperHa))
# Everything works until now
test=merge(SBA,PSS,all=T)
#test=merge(SBA,PSS,by=c("Year","stand","Plot","PlotID","AreaExpansion"))
PSS=subset(test,Litterfall_LAI!="NA")
# Merge with DM
PSS=merge(DM,PSS,by=c("Year","Plot","PlotID"))
##################    Add Thinned      ##########################
PSS$Thinned<-ifelse(PSS$Treatment=="NoThin",0,1)
# Calculates SBA/CL ratio for SAP.MK and SMAG.MK
PSS$CL.SBA=(PSS$CLperHa/PSS$SBAperHa)
# Calculates DBH/Ht ratio for SMAG models
PSS$DBH.Ht=(PSS$DBHperHa/PSS$HtperHa)
# Calculates BA/mLCR ratio
PSS$BAmLCR=PSS$BAperHa*PSS$mLCRperHa
PSS$BA.mLCR=PSS$BAperHa*PSS$MeanmLCR
head(PSS)
######################################################################
plot(LAI.spline~DBH.Ht,data=PSS)
plot(LAI.spline~BA.mLCR,data=PSS)
plot(LAI.spline~CL.SBA,data=PSS)
######################################################################################
par(mar=c(5,5,4,2))
# Plots LAI over Sapwood area per hectare
plot(LAI.spline~SBAperHa,data=PSS,Ylab="LAI",xlab="Sapwood Area per Hectare",type="n")
# Unthinned plots (1990, 1970, Thinning Control, HD, MM, Blue-Spring, Nutting)
points(PSS[PSS$Treatment=="NoThin",]$SBAperHa,PSS[PSS$Treatment=="NoThin",]$LAI.spline,col="black",pch=20)
points(PSS[PSS$Treatment=="PCT",]$SBAperHa,PSS[PSS$Treatment=="PCT",]$LAI.spline,col="yellow",pch=17)
points(PSS[PSS$PlotID=="17",]$SBAperHa,PSS[PSS$PlotID=="17",]$LAI.spline,col="blue",pch=3)
points(PSS[PSS$PlotID=="18",]$SBAperHa,PSS[PSS$PlotID=="18",]$LAI.spline,col="blue",pch=3)
points(PSS[PSS$PlotID=="19",]$SBAperHa,PSS[PSS$PlotID=="19",]$LAI.spline,col="blue",pch=3)
points(PSS[PSS$PlotID=="19",]$SBAperHa,PSS[PSS$PlotID=="19",]$LAI.spline,col="blue",pch=3)
points(PSS[PSS$PlotID=="20",]$SBAperHa,PSS[PSS$PlotID=="20",]$LAI.spline,col="blue",pch=3)
points(PSS[PSS$PlotID=="21",]$SBAperHa,PSS[PSS$PlotID=="21",]$LAI.spline,col="blue",pch=3)
points(PSS[PSS$Treatment=="LD",]$SBAperHa,PSS[PSS$Treatment=="LD",]$LAI.spline,col="green",pch=15)
points(PSS[PSS$Treatment=="B",]$SBAperHa,PSS[PSS$Treatment=="B",]$LAI.spline,col="purple",pch=18)
legend("bottomright", c("LD","B","NT","Spacing","1990"),cex=1, col=c("green","purple","black","yellow","blue"),pch=c(15,18,20,17,3))

################################################################################################
# Plots LAI over BAperHa
plot(LAI.spline~BAperHa,data=PSS,Ylab="LAI",xlab="Basal Area per Hectare",type="n")
points(PSS[PSS$Treatment=="NoThin",]$BAperHa,PSS[PSS$Treatment=="NoThin",]$LAI.spline,col="black",pch=20)
points(PSS[PSS$Treatment=="PCT",]$BAperHa,PSS[PSS$Treatment=="PCT",]$LAI.spline,col="yellow",pch=17)
points(PSS[PSS$PlotID=="17",]$BAperHa,PSS[PSS$PlotID=="17",]$LAI.spline,col="blue",pch=3)
points(PSS[PSS$PlotID=="18",]$BAperHa,PSS[PSS$PlotID=="18",]$LAI.spline,col="blue",pch=3)
points(PSS[PSS$PlotID=="19",]$BAperHa,PSS[PSS$PlotID=="19",]$LAI.spline,col="blue",pch=3)
points(PSS[PSS$PlotID=="19",]$BAperHa,PSS[PSS$PlotID=="19",]$LAI.spline,col="blue",pch=3)
points(PSS[PSS$PlotID=="20",]$BAperHa,PSS[PSS$PlotID=="20",]$LAI.spline,col="blue",pch=3)
points(PSS[PSS$PlotID=="21",]$BAperHa,PSS[PSS$PlotID=="21",]$LAI.spline,col="blue",pch=3)
points(PSS[PSS$Treatment=="LD",]$BAperHa,PSS[PSS$Treatment=="LD",]$LAI.spline,col="green",pch=15)
points(PSS[PSS$Treatment=="B",]$BAperHa,PSS[PSS$Treatment=="B",]$LAI.spline,col="purple",pch=18)
legend("bottomright", c("LD","B","NT","Spacing","1990"),cex=0.85, col=c("green","purple","black","yellow","blue"),pch=c(15,18,20,17,3))

# Plots LAI over CLperHa
plot(LAI.spline~CLperHa,data=PSS,Ylab="LAI",ylim=c(1,8),xlab="Crown length per Hectare (m)",type="n")
points(PSS[PSS$Treatment=="NoThin",]$CLperHa,PSS[PSS$Treatment=="NoThin",]$LAI.spline,col="black",pch=20)
points(PSS[PSS$Treatment=="PCT",]$CLperHa,PSS[PSS$Treatment=="PCT",]$LAI.spline,col="yellow",pch=17)
points(PSS[PSS$PlotID=="17",]$CLperHa,PSS[PSS$PlotID=="17",]$LAI.spline,col="blue",pch=3)
points(PSS[PSS$PlotID=="18",]$CLperHa,PSS[PSS$PlotID=="18",]$LAI.spline,col="blue",pch=3)
points(PSS[PSS$PlotID=="19",]$CLperHa,PSS[PSS$PlotID=="19",]$LAI.spline,col="blue",pch=3)
points(PSS[PSS$PlotID=="19",]$CLperHa,PSS[PSS$PlotID=="19",]$LAI.spline,col="blue",pch=3)
points(PSS[PSS$PlotID=="20",]$CLperHa,PSS[PSS$PlotID=="20",]$LAI.spline,col="blue",pch=3)
points(PSS[PSS$PlotID=="21",]$CLperHa,PSS[PSS$PlotID=="21",]$LAI.spline,col="blue",pch=3)
points(PSS[PSS$Treatment=="LD",]$CLperHa,PSS[PSS$Treatment=="LD",]$LAI.spline,col="green",pch=15)
points(PSS[PSS$Treatment=="B",]$CLperHa,PSS[PSS$Treatment=="B",]$LAI.spline,col="purple",pch=18)
legend("bottomright", c("LD","B","NT","Spacing","1990"),cex=0.9, col=c("green","purple","black","yellow","blue"),pch=c(15,18,20,17,3))
#
plot(LAI.spline~TOPHT,data=PSS)
plot(LAI.spline~DBHperHa,data=PSS)
plot(LAI.spline~HtperHa,data=PSS)
plot(LAI.spline)
###############    Equations    #########################
################################################################################
source("C:/Users/Frank/Desktop/Tinn R Files/Furnival.r")
rmse<-function(obs,pred) {
  sqrt(sum((obs-pred)^2,na.rm=T)/length(obs))
}
mab<-function(obs,pred) {
  (sum(abs(obs-pred),na.rm=T))/length(obs)
}
SAPNLS=gnls(LAI.spline~(a1*SBAperHa^a2),data=PSS,na.action=na.omit,start=c(a1=1,a2=.26),control=nlmeControl(minScale=1e-10,returnObject=T));summary(SAPNLS)
PSS$LAI.SAPNLS.pred=(coef(SAPNLS)[1]*(PSS$SBAperHa^(coef(SAPNLS)[2])))
rmse.LAI.SAPNLS<-rmse(PSS$LAI.spline,PSS$LAI.SAPNLS.pred);rmse.LAI.SAPNLS
mab.LAI.SAPNLS<-mab(PSS$LAI.spline,PSS$LAI.SAPNLS.pred);mab.LAI.SAPNLS
furnival(SAPNLS)
# SAP- AIC=68.70254- P-value=0.05
SAPNLS.TOPHT=gnls((LAI.spline)~(a0*TOPHT)+(a1*log10(TOPHT)+(SBAperHa^a2)),data=PSS,na.action=na.omit,start=c(a0=0.05,a1=-5,a2=.26),control=nlmeControl(minScale=1e-10,returnObject=T));summary(SAPNLS.TOPHT)
PSS$LAI.SAPNLS.TOPHT.pred=(coef(SAPNLS.TOPHT)[1]*(PSS$TOPHT)+(coef(SAPNLS.TOPHT)[2]*(log10(PSS$TOPHT)))+(PSS$SBAperHa^(coef(SAPNLS.TOPHT)[3])))
rmse.LAI.SAPNLS.TOPHT<-rmse(PSS$LAI.spline,PSS$LAI.SAPNLS.TOPHT.pred);rmse.LAI.SAPNLS.TOPHT
mab.LAI.SAPNLS.TOPHT<-mab(PSS$LAI.spline,PSS$LAI.SAPNLS.TOPHT.pred);mab.LAI.SAPNLS.TOPHT
furnival(SAPNLS.TOPHT)
# Random effect on plot- ###################  Try ranef on b1+b2  ################
SAPNLMER.TOPHT=nlme((LAI.spline)~(a0*TOPHT)+(a1*log10(TOPHT)+(SBAperHa^a2)),data=PSS,na.action=na.omit,verbose=T,fixed=a0+a1+a2~1,random=a2~1|Plot,start=c(a0=5,a1=9.7,a2=0.4),control=nlmeControl(minScale=1e-10,returnObject=T));summary(SAPNLMER.TOPHT)
PSS$LAI.SAPNLMER.TOPHT.pred=(fixef(SAPNLMER.TOPHT)[1]*(PSS$TOPHT)+(fixef(SAPNLMER.TOPHT)[2]*(log10(PSS$TOPHT)))+(PSS$SBAperHa^(fixef(SAPNLMER.TOPHT)[3])))
rmse.LAI.SAPNLMER.TOPHT<-rmse(PSS$LAI.spline,PSS$LAI.SAPNLMER.TOPHT.pred);rmse.LAI.SAPNLMER.TOPHT
mab.LAI.SAPNLMER.TOPHT<-mab(PSS$LAI.spline,PSS$LAI.SAPNLMER.TOPHT.pred);mab.LAI.SAPNLMER.TOPHT
furnival(SAPNLMER.TOPHT)
###########################################################################################################
# AIC=54.9, a1- P-value=0.0793
SMAGNLS=gnls((LAI.spline)~(a1*SBAperHa^a2)*exp(a3*(DBH.Ht)),data=PSS, na.action=na.omit,start=c(a1=0.1,a2=1.3,a3=-0.1),control=nlmeControl(minScale=1e-10,returnObject=T));summary(SMAGNLS)  
PSS$LAI.SMAGNLS.pred=(coef(SMAGNLS)[1]*(PSS$SBAperHa^(coef(SMAGNLS)[2]))*(exp(coef(SMAGNLS)[3]*(PSS$DBHperHa/PSS$HtperHa))))
rmse.LAI.SMAGNLS<-rmse(PSS$LAI.spline,PSS$LAI.SMAGNLS.pred);rmse.LAI.SMAGNLS
mab.LAI.SMAGNLS<-mab(PSS$LAI.spline,PSS$LAI.SMAGNLS.pred);mab.LAI.SMAGNLS
furnival(SMAGNLS)
# AIC =45.12
SMAGNLS.TOPHT=gnls((LAI.spline)~(a0*TOPHT)+(a1*log10(TOPHT)+(SBAperHa^a2)*exp(a3*(DBH.Ht))),data=PSS,na.action=na.omit,start=c(a0=1,a1=0.1,a2=1.3,a3=-0.1),control=nlmeControl(minScale=1e-10,returnObject=T));summary(SMAGNLS.TOPHT)
PSS$LAI.SMAGNLS.TOPHT.pred=(coef(SMAGNLS.TOPHT)[1]*(PSS$TOPHT)+(coef(SMAGNLS.TOPHT)[2]*(log10(PSS$TOPHT)))+(PSS$SBAperHa^(coef(SMAGNLS.TOPHT)[3]))*(exp(coef(SMAGNLS.TOPHT)[4]*(PSS$DBHperHa/PSS$HtperHa))))
rmse.LAI.SMAGNLS.TOPHT<-rmse(PSS$LAI.spline,PSS$LAI.SMAGNLS.TOPHT.pred);rmse.LAI.SMAGNLS.TOPHT
mab.LAI.SMAGNLS.TOPHT<-mab(PSS$LAI.spline,PSS$LAI.SMAGNLS.TOPHT.pred);mab.LAI.SMAGNLS.TOPHT
furnival(SMAGNLS.TOPHT)
#AIC= 42.31-- Best Model ###########
SMAGNLMER.TOPHT=nlme((LAI.spline)~(a0*TOPHT)+(a1*log10(TOPHT))+((SBAperHa^a2)*exp(a3*(DBH.Ht))),na.action=na.omit,data=PSS,fixed=a0+a1+a2+a3~1,random=a3~1|Plot,start=c(a0=1,a1=0.1,a2=1.3,a3=-0.1),control=nlmeControl(minScale=1e-10,returnObject=T));summary(SMAGNLMER.TOPHT)
PSS$LAI.SMAGNLS.TOPHT.pred=(fixef(SMAGNLMER.TOPHT)[1]*(PSS$TOPHT)+(fixef(SMAGNLMER.TOPHT)[2]*(log10(PSS$TOPHT)))+(PSS$SBAperHa^(fixef(SMAGNLMER.TOPHT)[3]))*(exp(fixef(SMAGNLMER.TOPHT)[4]*(PSS$DBHperHa/PSS$HtperHa))))
rmse.LAI.SMAGNLMER.TOPHT<-rmse(PSS$LAI.spline,PSS$LAI.SMAGNLS.TOPHT.pred);rmse.LAI.SMAGNLMER.TOPHT
mab.LAI.SMAGNLMER.TOPHT<-mab(PSS$LAI.spline,PSS$LAI.SMAGNLS.TOPHT.pred);mab.LAI.SMAGNLMER.TOPHT
furnival(SMAGNLMER.TOPHT)
# plot residuals
PSS.subset=subset(PSS,SBAperHa!="NA")
PSS.subset$SMAGNLMER.resid=resid(SMAGNLMER.TOPHT,na.action=na.omit)
plot(PSS.subset$LAI.spline,PSS.subset$SMAGNLMER.resid)
plot(PSS.subset$TOPHT,PSS.subset$SMAGNLMER.resid)
############################################################################################################################
################################################################################################################################
#AIC= 255.9139, P-value= 0.0944
DCLNLS=gnls((LAI.spline)~(a1*DBHperHa^a2*(CLperHa^a3)),data=PSS,na.action=na.omit,start=c(a1=0.3970,a2=0.6149,a3=1.6450),control=nlmeControl(minScale=1e-10,returnObject=T));summary(DCLNLS)

DCLNLMER=nlme((LAI.spline)~(a1*(DBHperHa^a2)*(CLperHa^a3)),data=PSS,na.action=na.omit,weights=varPower(0.5,form="DBHperHa"),fixed=a1+a2+a3~1,random=a3~1|Treatment,start=c(a1=0.06,a2=0.88,a3=1.97),control=nlmeControl(minScale=1e-10,returnObject=T))
summary(DCLNLMER) 
ranef(DCLNLMER)
fixef(DCLNLMER)
plot(DCLNLMER)
###########################################################################################################################
# AIC= 75.72799- a1 and a3 non-significant
SCLNLS=gnls((LAI.spline)~(a1*SBAperHa^a2*CLperHa^a3),data=PSS,na.action=na.omit,start=c(a1=0.15,a2=0.5,a3=1.6),control=nlmeControl(minScale=1e-10,returnObject=T));summary(SCLNLS)  
PSS$LAI.SCLNLS.pred=(coef(SCLNLS)[1]*(PSS$SBAperHa^(coef(SCLNLS)[2]))*(PSS$CLperHa)^(coef(SCLNLS)[3]))
rmse.LAI.SCLNLS<-rmse(PSS$LAI.spline,PSS$LAI.SCLNLS.pred);rmse.LAI.SCLNLS
mab.LAI.SCLNLS<-mab(PSS$LAI.spline,PSS$LAI.SCLNLS.pred);mab.LAI.SCLNLS
furnival(SCLNLS)
# AIC= 49.86204- All parameters significant P=0.05
SCLNLS.TOPHT=gnls((LAI.spline)~((a0*TOPHT)+(a1*log10(TOPHT))+SBAperHa^a2*CLperHa^a3),data=PSS,na.action=na.omit,start=c(a0=1,a1=0.15,a2=0.5,a3=1.6),control=nlmeControl(minScale=1e-10,returnObject=T));summary(SCLNLS.TOPHT)
PSS$LAI.SCLNLS.TOPHT.pred=(coef(SCLNLS.TOPHT)[1]*(PSS$TOPHT)+(coef(SCLNLS.TOPHT)[2]*(log10(PSS$TOPHT))+(PSS$SBAperHa^(coef(SCLNLS.TOPHT)[3]))*(PSS$CLperHa)^(coef(SCLNLS.TOPHT)[4])))
rmse.LAI.SCLNLS.TOPHT<-rmse(PSS$LAI.spline,PSS$LAI.SCLNLS.TOPHT.pred);rmse.LAI.SCLNLS.TOPHT
mab.LAI.SCLNLS.TOPHT<-mab(PSS$LAI.spline,PSS$LAI.SCLNLS.TOPHT.pred);mab.LAI.SCLNLS.TOPHT
furnival(SCLNLS.TOPHT)
# AIC= 49.22822###################  Try ranef on b1+b2  ################
SCLNLMER.TOPHT=nlme((LAI.spline)~((a0*TOPHT)+(a1*log10(TOPHT))+SBAperHa^a2*CLperHa^a3),data=PSS,na.action=na.omit,fixed=a0+a1+a2+a3~1,random=a2~1|Plot,start=c(a0=1,a1=0.15,a2=0.5,a3=1.6),control=nlmeControl(minScale=1e-10,returnObject=T));summary(SCLNLMER.TOPHT)
PSS$LAI.SCLNLMER.TOPHT.pred=(fixef(SCLNLMER.TOPHT)[1]*(PSS$TOPHT)+(fixef(SCLNLMER.TOPHT)[2]*(log10(PSS$TOPHT))+(PSS$SBAperHa^(fixef(SCLNLMER.TOPHT)[3]))*(PSS$CLperHa)^(fixef(SCLNLMER.TOPHT)[4])))
rmse.LAI.SCLNLMER.TOPHT<-rmse(PSS$LAI.spline,PSS$LAI.SCLNLMER.TOPHT.pred);rmse.LAI.SCLNLMER.TOPHT
mab.LAI.SCLNLMER.TOPHT<-mab(PSS$LAI.spline,PSS$LAI.SCLNLMER.TOPHT.pred);mab.LAI.SCLNLMER.TOPHT
furnival(SCLNLMER.TOPHT)


points(PSS$LAI.SCLNLMER.TOPHT.pred~PSS$LAI.spline,pch=2,col="blue")
legend("bottomright",c("TreePLA","Allometric LAI"),pch=c(3,2),col=c("red","blue"),cex=0.9)
legend("bottomright",c("LME","TreePLA","Allometric LAI"),pch=c(20,3,2),col=c("black","red","blue"),cex=0.9)
###########################################################################################################################
# AIC= AIC= 261.8
VALNLS=gnls((LAI.spline)~(a1*((BA.mLCR)^a2)),data=PSS, na.action=na.omit,start=c(a1=1.2,a2=0.84),control=nlmeControl(minScale=1e-10,returnObject=T));summary(VALNLS)                            
# AIC= 222.8
VALNLS.TOPHT=gnls((LAI.spline)~((a1*TOPHT)+(a2*log10(TOPHT))+((BA.mLCR)^a3)),data=PSS, na.action=na.omit,start=c(a1=1.2,a2=0.84,a3=1),control=nlmeControl(minScale=1e-10,returnObject=T)) 
summary(VALNLS.TOPHT)
# AIC= 185.8
VALNLMER=nlme((LAI.spline)~(a1*TOPHT)+(a2*log10(TOPHT))+((BA.mLCR)^a3),data=PSS,na.action=na.omit,start=c(a1=0.9054,a2=0.8374,a3=1),fixed=a1+a2+a3~1,random=a2~1|Plot,control=nlmeControl(minScale=1e-10,returnObject=T),weights=varPower(0.5,form="BAperHa"))                                                                                                                                                              
summary(VALNLMER)
plot(VALWNLS)
plot(resid$VALWNLS
     
#AIC= ???????????????????????????????????????????????
VALNLMER=nlme((LAI.spline)~(a1*(BAperHa*mLCRperHa)^a2),data=PSS,na.action=na.omit,weights=varPower(0.5,form="BAperHa"),      
fixed=a1+a2~1,random=a2~1|Location/Thinned,start=c(a1=0.01,a2=1.6),control=nlmeControl(minScale=1e-10,returnObject=T)) 
summary(VALNLMER)

#############  SAP.MK  ###############################3
# AIC= 75.72799, a1 and a3 nonsignificant!
SAPNLS.MK=gnls(((LAI.spline)~(a1*(SBAperHa^a2)*(CL.SBA)^a3)),data=PSS,na.action=na.omit,start=c(a1=0.3,a2=.1,a3=2),control=nlmeControl(minScale=1e-10,returnObject=T));summary(SAPNLS.MK)
PSS$LAI.SAPNLSMK.pred=(coef(SAPNLS.MK)[1]*(PSS$SBAperHa^(coef(SAPNLS.MK)[2]))*(PSS$CL.SBA)^(coef(SAPNLS.MK)[3]))
rmse.LAI.SAPNLS.MK<-rmse(PSS$LAI.spline,PSS$LAI.SAPNLSMK.pred);rmse.LAI.SAPNLS.MK
mab.LAI.SAPNLS.MK<-mab(PSS$LAI.spline,PSS$LAI.SAPNLSMK.pred);mab.LAI.SAPNLS.MK
furnival(SAPNLS.MK)
#AIC= 49.86204 a3 P-value=0.001
SAPNLS.MK.TOPHT=gnls((LAI.spline)~(a1*TOPHT)+(a2*log10(TOPHT))+(SBAperHa^a3)*((CL.SBA)^a4),data=PSS,na.action=na.omit,start=c(a1=0.3,a2=-10,a3=2,a4=.5),control=nlmeControl(minScale=1e-10,returnObject=T));summary(SAPNLS.MK.TOPHT)
PSS$LAI.SAPNLSMK.TOPHT.pred=(coef(SAPNLS.MK.TOPHT)[1]*(PSS$TOPHT)+(coef(SAPNLS.MK.TOPHT)[2]*(log10(PSS$TOPHT))+(PSS$SBAperHa^(coef(SAPNLS.MK.TOPHT)[3]))*(PSS$CL.SBA)^(coef(SAPNLS.MK.TOPHT)[4])))
rmse.LAI.SAPNLS.MK.TOPHT<-rmse(PSS$LAI.spline,PSS$LAI.SAPNLSMK.TOPHT.pred);rmse.LAI.SAPNLS.MK.TOPHT
mab.LAI.SAPNLS.MK.TOPHT<-mab(PSS$LAI.spline,PSS$LAI.SAPNLSMK.TOPHT.pred);mab.LAI.SAPNLS.MK.TOPHT          
furnival(SAPNLS.MK.TOPHT)
# AIC 47.7664
SAPNLMER.MK.TOPHT=nlme((LAI.spline)~(a1*TOPHT)+(a2*log10(TOPHT))+(SBAperHa^a3)*((CL.SBA)^a4),data=PSS,na.action=na.omit,fixed=a1+a2+a3+a4~1,random=a4~1|Plot,start=c(a1=0.3,a2=.1,a3=2,a4=1),control=nlmeControl(minScale=1e-10,returnObject=T));summary(SAPNLMER.MK.TOPHT)
PSS$LAI.SAPNLSMK.TOPHT.pred=(fixef(SAPNLMER.MK.TOPHT)[1]*(PSS$TOPHT)+(fixef(SAPNLMER.MK.TOPHT)[2]*(log10(PSS$TOPHT))+(PSS$SBAperHa^(fixef(SAPNLMER.MK.TOPHT)[3]))*(PSS$CL.SBA)^(fixef(SAPNLMER.MK.TOPHT)[4])))
rmse.LAI.SAPNLS.MK.TOPHT<-rmse(PSS$LAI.spline,PSS$LAI.SAPNLSMK.TOPHT.pred);rmse.LAI.SAPNLS.MK.TOPHT
mab.LAI.SAPNLS.MK.TOPHT<-mab(PSS$LAI.spline,PSS$LAI.SAPNLSMK.TOPHT.pred);mab.LAI.SAPNLS.MK.TOPHT
furnival(SAPNLMER.MK.TOPHT)
##########################  SMAG.MK  #############################################################
# AIC= 51.75297
SMAGNLS.MK=gnls((LAI.spline)~(a1*(SBAperHa^a2)*((CL.SBA)^a3)*(exp(DBH.Ht)^a4)),data=PSS,na.action=na.omit,start=c(a1=0.14,a2=2.13,a3=1.6,a4=.13),control=nlmeControl(minScale=1e-10,returnObject=T));summary(SMAGNLS.MK)
PSS$LAI.SMAGNLS.MK.pred=(coef(SMAGNLS.MK)[1]*(PSS$SBAperHa^(coef(SMAGNLS.MK)[2])))*((PSS$CL.SBA)^(coef(SMAGNLS.MK)[3])*(exp(PSS$DBH.Ht^(coef(SMAGNLS.MK)[4]))))
rmse.LAI.SMAGNLS.MK<-rmse(PSS$LAI.spline,PSS$LAI.SMAGNLS.MK.pred);rmse.LAI.SMAGNLS.MK
mab.LAI.SMAGNLS.MK<-mab(PSS$LAI.spline,PSS$LAI.SMAGNLS.MK.pred);mab.LAI.SMAGNLS.MK
furnival(SMAGNLS.MK)

# AIC=  43.82223
SMAGNLS.MK.TOPHT=gnls((LAI.spline)~(a0*TOPHT)+(a1*log10(TOPHT))+(SBAperHa^a2)*((CL.SBA)^a3)*(exp(DBH.Ht)^a4),data=PSS,na.action=na.omit,start=c(a0=.1,a1=-4,a2=0.8,a3=0.4,a4=3),control=nlmeControl(minScale=1e-10,returnObject=T));summary(SMAGNLS.MK.TOPHT)
PSS$LAI.SMAGNLS.TOPHT.pred=(coef(SMAGNLS.MK.TOPHT)[1]*(PSS$TOPHT)+(coef(SMAGNLS.MK.TOPHT)[2]*(log10(PSS$TOPHT)))+(PSS$SBAperHa^(coef(SMAGNLS.MK.TOPHT)[3]))*(PSS$CL.SBA^(coef(SMAGNLS.MK.TOPHT)[4]))*(exp(PSS$DBH.Ht)^(coef(SMAGNLS.MK.TOPHT)[5])))
rmse.LAI.SMAGNLS.MK.TOPHT<-rmse(PSS$LAI.spline,PSS$LAI.SMAGNLS.TOPHT.pred);rmse.LAI.SMAGNLS.MK.TOPHT
mab.LAI.SMAGNLS.MK.TOPHT<-mab(PSS$LAI.spline,PSS$LAI.SMAGNLS.TOPHT.pred);mab.LAI.SMAGNLS.MK.TOPHT
furnival(SMAGNLS.MK.TOPHT)
                                                                                                                                                           
# AIC=43.36634, 
SMAGNLMER.MK.TOPHT=nlme((LAI.spline)~(a0*TOPHT)+(a1*log10(TOPHT))+(SBAperHa^a2)*((CL.SBA)^a3)*(exp(DBH.Ht)^a4),data=PSS,na.action=na.omit,fixed=a0+a1+a2+a3+a4~1,random=a2~1|Plot,start=c(a0=.1,a1=-4,a2=0.8,a3=0.4,a4=3),control=nlmeControl(minScale=1e-10,returnObject=T));summary(SMAGNLMER.MK)
PSS$LAI.SMAGNLS.TOPHT.pred=(fixef(SMAGNLMER.MK.TOPHT)[1]*(PSS$TOPHT)+(fixef(SMAGNLMER.MK.TOPHT)[2]*(log10(PSS$TOPHT)))+(PSS$SBAperHa^(fixef(SMAGNLMER.MK.TOPHT)[3]))*(PSS$CL.SBA^(fixef(SMAGNLMER.MK.TOPHT)[4]))*(exp(PSS$DBH.Ht)^(fixef(SMAGNLMER.MK.TOPHT)[5])))
rmse.LAI.SMAGNLMER.TOPHT<-rmse(PSS$LAI.spline,PSS$LAI.SMAGNLS.TOPHT.pred);rmse.LAI.SMAGNLMER.TOPHT
mab.LAI.SMAGNLMER.TOPHT<-mab(PSS$LAI.spline,PSS$LAI.SMAGNLS.TOPHT.pred);mab.LAI.SMAGNLMER.TOPHT
furnival(SMAGNLMER.TOPHT)

# Plots allometric LAI data on plot
head(PSS)
PSS1=subset(PSS,SBAperHa!="NA")
PSS2=subset(PSS,,select=c(Year,Plot,PlotID,LAI.spline,LAI.SCLNLMER.TOPHT.pred))
PSS2=subset(PSS2,LAI.SCLNLMER.TOPHT.pred!="NA")
str(PSS2)
PSS2
write.csv(PSS2, file = ("C:/Users/Frank/Desktop/Allometric LAI.csv"), row.names = FALSE)  
PSS2
points(PSS2$LAI.SCLNLMER.TOPHT.pred~PSS2$LAI.spline,pch=2,col="blue")
lines(lowess(PSS2$LAI.spline,PSS2$LAI.SCLNLMER.TOPHT.pred),lty=5)
     legend("bottomright",c(expression(paste("PLA"["Tree"],sep="")),"Allometric LAI"),pch=c(3,2),col=c("red","blue"),cex=1.0)
     legend("bottomright",c("LME",expression(paste("PLA"["Tree"],sep="")),"Allometric LAI "),pch=c(20,3,2),col=c("black","red","blue"),cex=1)
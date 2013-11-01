library(RODBC);library(plyr);library(doBy);library(TTR);library(beanplot)
temp=odbcConnectAccess("C:/Users/Frank/Desktop/Demeritt WP- Bland.mdb")
LAI=sqlFetch(temp,"Plot_litterfall LAI_all years");Plots=sqlFetch(temp,"Plots")
DM=read.csv("C:/Users/Frank/Desktop/R Output/DM.csv")
#########################################################################################
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
#############################################################################################
# This uses litterfall collected 2 years after inventory year 
ifelse(LAI$Year<2010,LAI$Year-2,LAI$Year)
##############################################################################################
# Loads parameter tree leaf area predicted from different models and fitting techniques
##############################################################################################
#SAP NLS, WNLS, NLME
SAP.NLS1=sqlFetch(temp,"Tree_SAP_NLS_LA_2001");SAP.NLS8=sqlFetch(temp,"Tree_SAP_NLS_LA_2008");SAP.NLS11=sqlFetch(temp,"Tree_SAP_NLS_2011")
SAP.NLS<-rbind(SAP.NLS1,SAP.NLS8,SAP.NLS11)
SAP.WNLS1=sqlFetch(temp,"Tree_SAP_WNLS_LA_2001");SAP.WNLS8=sqlFetch(temp,"Tree_SAP_WNLS_LA_2008");SAP.WNLS11=sqlFetch(temp,"Tree_SAP_WNLS_2011")
SAP.WNLS<-rbind(SAP.WNLS1,SAP.WNLS8,SAP.WNLS11)
SAP.NLME1=sqlFetch(temp,"Tree_SAP_NLME_LA_2001");SAP.NLME8=sqlFetch(temp,"Tree_SAP_NLME_LA_2008");SAP.NLME11=sqlFetch(temp,"Tree_SAP_NLME_2011")
SAP.NLME<-rbind(SAP.NLME1,SAP.NLME8,SAP.NLME11)
# SMAG NLS, WNLS, NLME
SMAG.NLS1=sqlFetch(temp,"Tree_SMAG_NLS_2001");SMAG.NLS8=sqlFetch(temp,"Tree_SMAG_NLS_2008");SMAG.NLS11=sqlFetch(temp,"Tree_SMAG_NLS_LA")
SMAG.NLS<-rbind(SMAG.NLS1,SMAG.NLS8,SMAG.NLS11)
SMAG.WNLS1=sqlFetch(temp,"Tree_SMAG_WNLS_2001");SMAG.WNLS8=sqlFetch(temp,"Tree_SMAG_WNLS_2008");SMAG.WNLS11=sqlFetch(temp,"Tree_SMAG_LA_2011")
SMAG.WNLS<-rbind(SMAG.WNLS1,SMAG.WNLS8,SMAG.WNLS11)
SMAG.NLME1=sqlFetch(temp,"Tree_SMAG_NLME_2001");SMAG.NLME8=sqlFetch(temp,"Tree_SMAG_NLME_2008");SMAG.NLME11=sqlFetch(temp,"Tree_SMAG_NLME_LA")
SMAG.NLME<-rbind(SMAG.NLME1,SMAG.NLME8,SMAG.NLME11)
#SCL NLS, WNLS, NLME 
SCL.NLS1=sqlFetch(temp,"Tree_SCL_NLS_2001");SCL.NLS8=sqlFetch(temp,"Tree_SCL_NLS_2008");SCL.NLS11=sqlFetch(temp,"Tree_SCL_NLS_2011")
SCL.NLS=rbind(SCL.NLS1,SCL.NLS8,SCL.NLS11)
SCL.WNLS1=sqlFetch(temp,"Tree_SCL_WNLS_2001");SCL.WNLS8=sqlFetch(temp,"Tree_SCL_WNLS_2008");SCL.WNLS11=sqlFetch(temp,"Tree_SCL_WNLS_2011")
SCL.WNLS=rbind(SCL.WNLS1,SCL.WNLS8,SCL.WNLS11)
SCL.NLME1=sqlFetch(temp,"Tree_SCL_NLME_2001");SCL.NLME8=sqlFetch(temp,"Tree_SCL_NLME_2008");SCL.NLME11=sqlFetch(temp,"Tree_SCL_NLME_2011")
SCL.NLME=rbind(SCL.NLME1,SCL.NLME8,SCL.NLME11)
# MAG
MAG.NLS=sqlFetch(temp,"Tree_MAG_NLS_LA");MAG.WNLS=sqlFetch(temp,"Tree_MAG_WNLS_LA");MAG.NLME=sqlFetch(temp,"Tree_MAG_LA_all years")
# VAL
VAL.NLS=sqlFetch(temp,"Tree_VAL_NLS_LA");VAL.WNLS=sqlFetch(temp,"Tree_VAL_LA_all years");VAL.NLME=sqlFetch(temp,"Tree_VAL_NLME_LA")
# DCL
DCL.NLS=sqlFetch(temp,"Tree_DCL_NLS_LA");DCL.WNLS=sqlFetch(temp,"Tree_DCL_LA_all years");DCL.NLME=sqlFetch(temp,"Tree_DCL_NLME_LA")
# SMAG NLS, WNLS, NLME
SMAGMK.NLS1=sqlFetch(temp,"Tree_SMAGMK_NLS_2001");SMAGMK.NLS8=sqlFetch(temp,"Tree_SMAGMK_NLS_2008");SMAGMK.NLS11=sqlFetch(temp,"Tree_SMAGMK_NLS_LA")
SMAGMK.NLS<-rbind(SMAGMK.NLS1,SMAGMK.NLS8,SMAGMK.NLS11)
SMAGMK.WNLS1=sqlFetch(temp,"Tree_SMAGMK_WNLS_2001");SMAGMK.WNLS8=sqlFetch(temp,"Tree_SMAGMK_WNLS_2008");SMAGMK.WNLS11=sqlFetch(temp,"Tree_SMAGMK_WNLS_LA")
SMAGMK.WNLS<-rbind(SMAGMK.WNLS1,SMAGMK.WNLS8,SMAGMK.WNLS11)
SMAGMK.NLME1=sqlFetch(temp,"Tree_SMAGMK_NLME_2001");SMAGMK.NLME8=sqlFetch(temp,"Tree_SMAGMK_NLME_2008");SMAGMK.NLME11=sqlFetch(temp,"Tree_SMAGMK_NLME_LA")
SMAGMK.NLME<-rbind(SMAGMK.NLME1,SMAGMK.NLME8,SMAGMK.NLME11)
# SAP.MK
SAPMK.NLS1=sqlFetch(temp,"Tree_SAPMK_NLS_2001");SAPMK.NLS8=sqlFetch(temp,"Tree_SAPMK_NLS_2008");SAPMK.NLS11=sqlFetch(temp,"Tree_SAPMK_NLS_LA")
SAPMK.NLS<-rbind(SAPMK.NLS1,SAPMK.NLS8,SAPMK.NLS11)
SAPMK.WNLS1=sqlFetch(temp,"Tree_SAPMK_WNLS_2001");SAPMK.WNLS8=sqlFetch(temp,"Tree_SAPMK_WNLS_2008");SAPMK.WNLS11=sqlFetch(temp,"Tree_SAPMK_WNLS_LA")
SAPMK.WNLS<-rbind(SAPMK.WNLS1,SAPMK.WNLS8,SAPMK.WNLS11)
SAPMK.NLME1=sqlFetch(temp,"Tree_SAPMK_NLME_2001");SAPMK.NLME8=sqlFetch(temp,"Tree_SAPMK_NLME_2008");SAPMK.NLME11=sqlFetch(temp,"Tree_SAPMK_NLME_LA")
SAPMK.NLME<-rbind(SAPMK.NLME1,SAPMK.NLME8,SAPMK.NLME11)
close(temp)
############################################################################################################
# Removes trees with no LA value- dead/too small to core/negative leaf area predicted from sapwood
SMAG.NLS=subset(SMAG.NLS,!is.na(SMAG_LA),select=,);SMAG.WNLS=subset(SMAG.WNLS,!is.na(SMAG_LA),select=,);SMAG.NLME=subset(SMAG.NLME,!is.na(SMAG_LA),select=,)                     
SAP.NLS=subset(SAP.NLS,!is.na(SAP_LA),select=,);SAP.WNLS=subset(SAP.WNLS,!is.na(SAP_LA),select=,);SAP.NLME=subset(SAP.NLME,!is.na(SAP_LA),select=,)                     
MAG.NLS=subset(MAG.NLS,!is.na(MAG_LA),select=,);MAG.WNLS=subset(MAG.WNLS,!is.na(MAG_LA),select=,);MAG.NLME=subset(MAG.NLME,!is.na(MAG_LA),select=,)                    
SCL.NLS=subset(SCL.NLS,!is.na(SCL_LA),select=,);SCL.WNLS=subset(SCL.WNLS,!is.na(SCL_LA),select=,);SCL.NLME=subset(SCL.NLME,!is.na(SCL_LA),select=,)                     
VAL.NLS=subset(VAL.NLS,!is.na(VAL_LA),select=,);VAL.WNLS=subset(VAL.WNLS,!is.na(VAL_LA),select=,);VAL.NLME=subset(VAL.NLME,!is.na(VAL_LA),select=,)
DCL.NLS=subset(DCL.NLS,!is.na(DCL_LA),select=,);DCL.WNLS=subset(DCL.WNLS,!is.na(DCL_LA),select=,);DCL.NLME=subset(DCL.NLME,!is.na(DCL_LA),select=,)
SMAGMK.NLS=subset(SMAGMK.NLS,!is.na(SMAGMK_LA),select=,);SMAGMK.WNLS=subset(SMAGMK.WNLS,!is.na(SMAGMK_LA),select=,);SMAGMK.NLME=subset(SMAGMK.NLME,!is.na(SMAGMK_LA),select=,)
SAPMK.NLS=subset(SAPMK.NLS,!is.na(SAPMK_LA),select=,);SAPMK.WNLS=subset(SAPMK.WNLS,!is.na(SAPMK_LA),select=,);SAPMK.NLME=subset(SAPMK.NLME,!is.na(SAPMK_LA),select=,)

#SMAG.LA=ddply(SMAG,.(Year,PlotID,Plot),summarise,sum('SMAG_LA')) -didn't work
# Sums tree leaf area to the plot-level
SMAG.NLS.LA=summaryBy(SMAG_LA~Year+PlotID, data=SMAG.NLS, FUN=sum);SMAG.WNLS.LA=summaryBy(SMAG_LA~Year+PlotID, data=SMAG.WNLS, FUN=sum);SMAG.NLME.LA=summaryBy(SMAG_LA~Year+PlotID, data=SMAG.NLME, FUN=sum)
SAP.NLS.LA=summaryBy(SAP_LA~Year+PlotID, data=SAP.NLS, FUN=sum);SAP.WNLS.LA=summaryBy(SAP_LA~Year+PlotID, data=SAP.WNLS, FUN=sum);SAP.NLME.LA=summaryBy(SAP_LA~Year+PlotID, data=SAP.NLME, FUN=sum)                     
MAG.NLS.LA=summaryBy(MAG_LA~Year+PlotID, data=MAG.NLS, FUN=sum);MAG.WNLS.LA=summaryBy(MAG_LA~Year+PlotID, data=MAG.WNLS, FUN=sum);MAG.NLME.LA=summaryBy(MAG_LA~Year+PlotID, data=MAG.NLME, FUN=sum)
SCL.NLS.LA=summaryBy(SCL_LA~Year+PlotID, data=SCL.NLS, FUN=sum);SCL.WNLS.LA=summaryBy(SCL_LA~Year+PlotID, data=SCL.WNLS, FUN=sum);SCL.NLME.LA=summaryBy(SCL_LA~Year+PlotID, data=SCL.NLME, FUN=sum)                     
VAL.NLS.LA=summaryBy(VAL_LA~Year+PlotID, data=VAL.NLS, FUN=sum);VAL.WNLS.LA=summaryBy(VAL_LA~Year+PlotID, data=VAL.WNLS, FUN=sum);VAL.NLME.LA=summaryBy(VAL_LA~Year+PlotID, data=VAL.NLME, FUN=sum)                     
DCL.NLS.LA=summaryBy(DCL_LA~Year+PlotID, data=DCL.NLS, FUN=sum);DCL.WNLS.LA=summaryBy(DCL_LA~Year+PlotID, data=DCL.WNLS, FUN=sum);DCL.NLME.LA=summaryBy(DCL_LA~Year+PlotID, data=DCL.NLME, FUN=sum)
SMAGMK.NLS.LA=summaryBy(SMAGMK_LA~Year+PlotID, data=SMAGMK.NLS, FUN=sum);SMAGMK.WNLS.LA=summaryBy(SMAGMK_LA~Year+PlotID, data=SMAGMK.WNLS, FUN=sum);SMAGMK.NLME.LA=summaryBy(SMAGMK_LA~Year+PlotID, data=SMAGMK.NLME, FUN=sum)
SAPMK.NLS.LA=summaryBy(SAPMK_LA~Year+PlotID, data=SAPMK.NLS, FUN=sum);SAPMK.WNLS.LA=summaryBy(SAPMK_LA~Year+PlotID, data=SAPMK.WNLS, FUN=sum);SAPMK.NLME.LA=summaryBy(SAPMK_LA~Year+PlotID, data=SAPMK.NLME, FUN=sum)

#Merges allometric LA with plot attribute table to calculate LAI 
SMAGLA.NLS=merge(SMAG.NLS.LA,Plots,by='PlotID',all=T);SMAGLA.WNLS=merge(SMAG.WNLS.LA,Plots,by='PlotID',all=T);SMAGLA.NLME=merge(SMAG.NLME.LA,Plots,by='PlotID',all=T)                     
SAPLA.NLS=merge(SAP.NLS.LA,Plots,by='PlotID',all=T);SAPLA.WNLS=merge(SAP.WNLS.LA,Plots,by='PlotID',all=T);SAPLA.NLME=merge(SAP.NLME.LA,Plots,by='PlotID',all=T)                     
MAGLA.NLS=merge(MAG.NLS.LA,Plots,by='PlotID',all=T);MAGLA.WNLS=merge(MAG.WNLS.LA,Plots,by='PlotID',all=T);MAGLA.NLME=merge(MAG.NLME.LA,Plots,by='PlotID',all=T)                     
SCLLA.NLS=merge(SCL.NLS.LA,Plots,by='PlotID',all=T);SCLLA.WNLS=merge(SCL.WNLS.LA,Plots,by='PlotID',all=T);SCLLA.NLME=merge(SCL.NLME.LA,Plots,by='PlotID',all=T)                    
VALLA.NLS=merge(VAL.NLS.LA,Plots,by='PlotID',all=T);VALLA.WNLS=merge(VAL.WNLS.LA,Plots,by='PlotID',all=T);VALLA.NLME=merge(VAL.NLME.LA,Plots,by='PlotID',all=T)                     
DCLLA.NLS=merge(DCL.NLS.LA,Plots,by='PlotID',all=T);DCLLA.WNLS=merge(DCL.WNLS.LA,Plots,by='PlotID',all=T);DCLLA.NLME=merge(DCL.NLME.LA,Plots,by='PlotID',all=T)                     
SAPMKLA.NLS=merge(SAPMK.NLS.LA,Plots,by='PlotID',all=T);SAPMKLA.WNLS=merge(SAPMK.WNLS.LA,Plots,by='PlotID',all=T);SAPMKLA.NLME=merge(SAPMK.NLME.LA,Plots,by='PlotID',all=T)
SMAGMKLA.NLS=merge(SMAGMK.NLS.LA,Plots,by='PlotID',all=T);SMAGMKLA.WNLS=merge(SMAGMK.WNLS.LA,Plots,by='PlotID',all=T);SMAGMKLA.NLME=merge(SMAGMK.NLME.LA,Plots,by='PlotID',all=T)
# Calculate LAI- (Leaf Area per plot(m2))/Plot Area (m)
SMAGLA.NLS$LAI=SMAGLA.NLS$SMAG_LA.sum/SMAGLA.NLS$AreaMeters;SMAGLA.WNLS$LAI=SMAGLA.WNLS$SMAG_LA.sum/SMAGLA.WNLS$AreaMeters;SMAGLA.NLME$LAI=SMAGLA.NLME$SMAG_LA.sum/SMAGLA.NLME$AreaMeters                     
SAPLA.NLS$LAI=SAPLA.NLS$SAP_LA.sum/SAPLA.NLS$AreaMeters;SAPLA.WNLS$LAI=SAPLA.WNLS$SAP_LA.sum/SAPLA.WNLS$AreaMeters;SAPLA.NLME$LAI=SAPLA.NLME$SAP_LA.sum/SAPLA.NLME$AreaMeters
MAGLA.NLS$LAI=MAGLA.NLS$MAG_LA.sum/MAGLA.NLS$AreaMeters;MAGLA.WNLS$LAI=MAGLA.WNLS$MAG_LA.sum/MAGLA.WNLS$AreaMeters;MAGLA.NLME$LAI=MAGLA.NLME$MAG_LA.sum/MAGLA.NLME$AreaMeters                     
SCLLA.NLS$LAI=SCLLA.NLS$SCL_LA.sum/SCLLA.NLS$AreaMeters;SCLLA.WNLS$LAI=SCLLA.WNLS$SCL_LA.sum/SCLLA.WNLS$AreaMeters;SCLLA.NLME$LAI=SCLLA.NLME$SCL_LA.sum/SCLLA.NLME$AreaMeters                     
VALLA.NLS$LAI=VALLA.NLS$VAL_LA.sum/VALLA.NLS$AreaMeters;VALLA.WNLS$LAI=VALLA.WNLS$VAL_LA.sum/VALLA.WNLS$AreaMeters;VALLA.NLME$LAI=VALLA.NLME$VAL_LA.sum/VALLA.NLME$AreaMeters                     
DCLLA.NLS$LAI=DCLLA.NLS$DCL_LA.sum/DCLLA.NLS$AreaMeters;DCLLA.WNLS$LAI=DCLLA.WNLS$DCL_LA.sum/DCLLA.WNLS$AreaMeters;DCLLA.NLME$LAI=DCLLA.NLME$DCL_LA.sum/DCLLA.NLME$AreaMeters   
SAPMKLA.NLS$LAI=SAPMKLA.NLS$SAPMK_LA.sum/SAPMKLA.NLS$AreaMeters;SAPMKLA.WNLS$LAI=SAPMKLA.WNLS$SAPMK_LA.sum/SAPMKLA.WNLS$AreaMeters;SAPMKLA.NLME$LAI=SAPMKLA.NLME$SAPMK_LA.sum/SAPMKLA.NLME$AreaMeters
SMAGMKLA.NLS$LAI=SMAGMKLA.NLS$SMAGMK_LA.sum/SMAGMKLA.NLS$AreaMeters;SMAGMKLA.WNLS$LAI=SMAGMKLA.WNLS$SMAGMK_LA.sum/SMAGMKLA.WNLS$AreaMeters;SMAGMKLA.NLME$LAI=SMAGMKLA.NLME$SMAGMK_LA.sum/SMAGMKLA.NLME$AreaMeters                     

# This gives Allometric LAI by Year and Plot
SMAG.NLS.LAI=aggregate(x=SMAGLA.NLS$LAI,by=list(SMAGLA.NLS$Year,SMAGLA.NLS$Plot,SMAGLA.NLS$PlotID),FUN=sum)
SMAG.WNLS.LAI=aggregate(x=SMAGLA.WNLS$LAI,by=list(SMAGLA.WNLS$Year,SMAGLA.WNLS$Plot,SMAGLA.WNLS$PlotID),FUN=sum)                     
SMAG.NLME.LAI=aggregate(x=SMAGLA.NLME$LAI,by=list(SMAGLA.NLME$Year,SMAGLA.NLME$Plot,SMAGLA.NLME$PlotID),FUN=sum)                     
SAP.NLS.LAI=aggregate(x=SAPLA.NLS$LAI,by=list(SAPLA.NLS$Year,SAPLA.NLS$Plot,SAPLA.NLS$PlotID),FUN=sum)
SAP.WNLS.LAI=aggregate(x=SAPLA.WNLS$LAI,by=list(SAPLA.WNLS$Year,SAPLA.WNLS$Plot,SAPLA.WNLS$PlotID),FUN=sum)                     
SAP.NLME.LAI=aggregate(x=SAPLA.NLME$LAI,by=list(SAPLA.NLME$Year,SAPLA.NLME$Plot,SAPLA.NLME$PlotID),FUN=sum)                     
MAG.NLS.LAI=aggregate(x=MAGLA.NLS$LAI,by=list(MAGLA.NLS$Year,MAGLA.NLS$Plot,MAGLA.NLS$PlotID),FUN=sum)
MAG.WNLS.LAI=aggregate(x=MAGLA.WNLS$LAI,by=list(MAGLA.WNLS$Year,MAGLA.WNLS$Plot,MAGLA.WNLS$PlotID),FUN=sum)                     
MAG.NLME.LAI=aggregate(x=MAGLA.NLME$LAI,by=list(MAGLA.NLME$Year,MAGLA.NLME$Plot,MAGLA.NLME$PlotID),FUN=sum)                     
SCL.NLS.LAI=aggregate(x=SCLLA.NLS$LAI,by=list(SCLLA.NLS$Year,SCLLA.NLS$Plot,SCLLA.NLS$PlotID),FUN=sum)
SCL.WNLS.LAI=aggregate(x=SCLLA.WNLS$LAI,by=list(SCLLA.WNLS$Year,SCLLA.WNLS$Plot,SCLLA.WNLS$PlotID),FUN=sum)                     
SCL.NLME.LAI=aggregate(x=SCLLA.NLME$LAI,by=list(SCLLA.NLME$Year,SCLLA.NLME$Plot,SCLLA.NLME$PlotID),FUN=sum)                     
VAL.NLS.LAI=aggregate(x=VALLA.NLS$LAI,by=list(VALLA.NLS$Year,VALLA.NLS$Plot,VALLA.NLS$PlotID),FUN=sum)
VAL.WNLS.LAI=aggregate(x=VALLA.WNLS$LAI,by=list(VALLA.WNLS$Year,VALLA.WNLS$Plot,VALLA.WNLS$PlotID),FUN=sum)                     
VAL.NLME.LAI=aggregate(x=VALLA.NLME$LAI,by=list(VALLA.NLME$Year,VALLA.NLME$Plot,VALLA.NLME$PlotID),FUN=sum)                     
DCL.NLS.LAI=aggregate(x=DCLLA.NLS$LAI,by=list(DCLLA.NLS$Year,DCLLA.NLS$Plot,DCLLA.NLS$PlotID),FUN=sum)
DCL.WNLS.LAI=aggregate(x=DCLLA.WNLS$LAI,by=list(DCLLA.WNLS$Year,DCLLA.WNLS$Plot,DCLLA.WNLS$PlotID),FUN=sum)                     
DCL.NLME.LAI=aggregate(x=DCLLA.NLME$LAI,by=list(DCLLA.NLME$Year,DCLLA.NLME$Plot,DCLLA.NLME$PlotID),FUN=sum)     
SMAGMK.NLS.LAI=aggregate(x=SMAGMKLA.NLS$LAI,by=list(SMAGMKLA.NLS$Year,SMAGMKLA.NLS$Plot,SMAGMKLA.NLS$PlotID),FUN=sum)
SMAGMK.WNLS.LAI=aggregate(x=SMAGMKLA.WNLS$LAI,by=list(SMAGMKLA.WNLS$Year,SMAGMKLA.WNLS$Plot,SMAGMKLA.WNLS$PlotID),FUN=sum)                     
SMAGMK.NLME.LAI=aggregate(x=SMAGMKLA.NLME$LAI,by=list(SMAGMKLA.NLME$Year,SMAGMKLA.NLME$Plot,SMAGMKLA.NLME$PlotID),FUN=sum)                     
SAPMK.NLS.LAI=aggregate(x=SAPMKLA.NLS$LAI,by=list(SAPMKLA.NLS$Year,SAPMKLA.NLS$Plot,SAPMKLA.NLS$PlotID),FUN=sum)
SAPMK.WNLS.LAI=aggregate(x=SAPMKLA.WNLS$LAI,by=list(SAPMKLA.WNLS$Year,SAPMKLA.WNLS$Plot,SAPMKLA.WNLS$PlotID),FUN=sum)                     
SAPMK.NLME.LAI=aggregate(x=SAPMKLA.NLME$LAI,by=list(SAPMKLA.NLME$Year,SAPMKLA.NLME$Plot,SAPMKLA.NLME$PlotID),FUN=sum)

#Adds names to new dataframes
names(SMAG.NLS.LAI)[1:4]<-c("Year","Plot","PlotID","SMAG.NLS.LAI");names(SMAG.WNLS.LAI)[1:4]<-c("Year","Plot","PlotID","SMAG.WNLS.LAI");names(SMAG.NLME.LAI)[1:4] <- c("Year","Plot","PlotID","SMAG.NLME.LAI")                     
names(SAP.NLS.LAI)[1:4]<-c("Year","Plot","PlotID","SAP.NLS.LAI");names(SAP.WNLS.LAI)[1:4]<-c("Year","Plot","PlotID","SAP.WNLS.LAI");names(SAP.NLME.LAI)[1:4] <- c("Year","Plot","PlotID","SAP.NLME.LAI")                     
names(MAG.NLS.LAI)[1:4]<-c("Year","Plot","PlotID","MAG.NLS.LAI");names(MAG.WNLS.LAI)[1:4]<-c("Year","Plot","PlotID","MAG.WNLS.LAI");names(MAG.NLME.LAI)[1:4] <- c("Year","Plot","PlotID","MAG.NLME.LAI")                     
names(SCL.NLS.LAI)[1:4]<-c("Year","Plot","PlotID","SCL.NLS.LAI");names(SCL.WNLS.LAI)[1:4]<-c("Year","Plot","PlotID","SCL.WNLS.LAI");names(SCL.NLME.LAI)[1:4] <- c("Year","Plot","PlotID","SCL.NLME.LAI")                     
names(VAL.NLS.LAI)[1:4]<-c("Year","Plot","PlotID","VAL.NLS.LAI");names(VAL.WNLS.LAI)[1:4]<-c("Year","Plot","PlotID","VAL.WNLS.LAI");names(VAL.NLME.LAI)[1:4] <- c("Year","Plot","PlotID","VAL.NLME.LAI")                     
names(DCL.NLS.LAI)[1:4]<-c("Year","Plot","PlotID","DCL.NLS.LAI");names(DCL.WNLS.LAI)[1:4]<-c("Year","Plot","PlotID","DCL.WNLS.LAI");names(DCL.NLME.LAI)[1:4] <- c("Year","Plot","PlotID","DCL.NLME.LAI")  
names(SMAGMK.NLS.LAI)[1:4]<-c("Year","Plot","PlotID","SMAGMK.NLS.LAI");names(SMAGMK.WNLS.LAI)[1:4]<-c("Year","Plot","PlotID","SMAGMK.WNLS.LAI");names(SMAGMK.NLME.LAI)[1:4] <- c("Year","Plot","PlotID","SMAGMK.NLME.LAI")                     
names(SAPMK.NLS.LAI)[1:4]<-c("Year","Plot","PlotID","SAPMK.NLS.LAI");names(SAPMK.WNLS.LAI)[1:4]<-c("Year","Plot","PlotID","SAPMK.WNLS.LAI");names(SAPMK.NLME.LAI)[1:4] <- c("Year","Plot","PlotID","SAPMK.NLME.LAI")

###########################################################################################################
# Now we have allometric LAI for all models
#############################################################################################
##############################  Merge  #########################################
# Combine Litterfall LAI with Allometric- 
head(LAI)
NLS.LAI1=merge(SAP.NLS.LAI,SCL.NLS.LAI,by=c('PlotID','Plot','Year'),na.action=na.omit,all=T);NLS.LAI2=merge(NLS.LAI1,SMAG.NLS.LAI,by=c('PlotID','Plot','Year'),na.action=na.omit,all=T)
NLS.LAI3=merge(NLS.LAI2,VAL.NLS.LAI,by=c('PlotID','Plot','Year'),na.action=na.omit,all=T);NLS.LAI4=merge(NLS.LAI3,MAG.NLS.LAI,by=c('PlotID','Plot','Year'),na.action=na.omit,all=T)
NLS.LAI5=merge(NLS.LAI4,DCL.NLS.LAI,by=c('PlotID','Plot','Year'),na.action=na.omit,all=T);NLS.LAI6=merge(NLS.LAI5,SMAGMK.NLS.LAI,by=c('PlotID','Plot','Year'),na.action=na.omit,all=T)
NLS.LAIall=merge(NLS.LAI6,SAPMK.NLS.LAI,by=c('PlotID','Plot','Year'),na.action=na.omit,all=T)
head(NLS.LAIall)
# Checks to see if 2008 trees were included
test=subset(NLS.LAIall,,select=c(Year,Plot,PlotID,SAP.NLS.LAI))
factor(NLS.LAIall$Year)
factor(LAI$Year)
##############################################################################################################
NLS=merge(NLS.LAIall,LAI,by=c('PlotID','Plot','Year'),na.action=na.omit,all=T)

# Calculate absolute value of LAI Diff
NLS$SAP.Dif<-abs(NLS$SAP.NLS.LAI-NLS$LAI.spline);NLS$SCL.Dif<-abs(NLS$SCL.NLS.LAI-NLS$LAI.spline)
NLS$SMAG.Dif<-abs(NLS$SMAG.NLS.LAI-NLS$LAI.spline);NLS$MAG.Dif<-abs(NLS$MAG.NLS.LAI-NLS$LAI.spline)
NLS$VAL.Dif<-abs(NLS$VAL.NLS.LAI-NLS$LAI.spline);NLS$DCL.Dif<-abs(NLS$DCL.NLS.LAI-NLS$LAI.spline)
NLS$SMAGMK.Dif<-abs(NLS$SMAGMK.NLS.LAI-NLS$LAI.spline);NLS$SAPMK.Dif<-abs(NLS$SAPMK.NLS.LAI-NLS$LAI.spline)
NLS=subset(NLS,DCL.Dif!="NA")

# Sum total diff/count of observations
#AAD.SAP=aggregate(x=NLS$SAP.Dif,by=list(NLS$Tree,NLS$Dbh,NLS$HtCB,tree$TotHt),FUN=sum)

WNLS.LAI1=merge(SAP.WNLS.LAI,SCL.WNLS.LAI,by=c('PlotID','Plot','Year'),na.action=na.omit,all=T)
WNLS.LAI2=merge(WNLS.LAI1,SMAG.WNLS.LAI,by=c('PlotID','Plot','Year'),na.action=na.omit,all=T)
WNLS.LAI3=merge(WNLS.LAI2,VAL.WNLS.LAI,by=c('PlotID','Plot','Year'),na.action=na.omit,all=T)
WNLS.LAI4=merge(WNLS.LAI3,MAG.WNLS.LAI,by=c('PlotID','Plot','Year'),na.action=na.omit,all=T)
WNLS.LAI5=merge(WNLS.LAI4,DCL.WNLS.LAI,by=c('PlotID','Plot','Year'),na.action=na.omit,all=T)
WNLS.LAI6=merge(WNLS.LAI5,SMAGMK.WNLS.LAI,by=c('PlotID','Plot','Year'),na.action=na.omit,all=T)
WNLS.LAIall=merge(WNLS.LAI6,SAPMK.WNLS.LAI,by=c('PlotID','Plot','Year'),na.action=na.omit,all=T)
head(WNLS.LAIall)
WNLS=merge(WNLS.LAIall,LAI,by=c('PlotID','Plot','Year'),na.action=na.omit,all=T)
str(WNLS)
factor(WNLS$Year)
# Plot Observed-predicted
# plot(WNLS$LAI.spline~WNLS$SAP.WNLS.LAI)
# plot(WNLS$LAI.spline~WNLS$SCL.WNLS.LAI)
# plot(WNLS$LAI.spline~WNLS$VAL.WNLS.LAI)
# plot(WNLS$LAI.spline~WNLS$MAG.WNLS.LAI)
# plot(WNLS$LAI.spline~WNLS$DCL.WNLS.LAI)

WNLS$SAP.Dif<-abs(WNLS$SAP.WNLS.LAI-WNLS$LAI.spline);WNLS$SCL.Dif<-abs(WNLS$SCL.WNLS.LAI-WNLS$LAI.spline)
WNLS$SMAG.Dif<-abs(WNLS$SMAG.WNLS.LAI-WNLS$LAI.spline);WNLS$MAG.Dif<-abs(WNLS$MAG.WNLS.LAI-WNLS$LAI.spline)
WNLS$VAL.Dif<-abs(WNLS$VAL.WNLS.LAI-WNLS$LAI.spline);WNLS$DCL.Dif<-abs(WNLS$DCL.WNLS.LAI-WNLS$LAI.spline)
WNLS$SMAGMK.Dif<-abs(WNLS$SMAGMK.WNLS.LAI-WNLS$LAI.spline);WNLS$SAPMK.Dif<-abs(WNLS$SAPMK.WNLS.LAI-WNLS$LAI.spline)
#WNLS=subset(WNLS,DCL.Dif!="NA")

NLME.LAI1=merge(SAP.NLME.LAI,SCL.NLME.LAI,by=c('PlotID','Plot','Year'),na.action=na.omit,all=T)
NLME.LAI2=merge(NLME.LAI1,SMAG.NLME.LAI,by=c('PlotID','Plot','Year'),na.action=na.omit,all=T)
NLME.LAI3=merge(NLME.LAI2,VAL.NLME.LAI,by=c('PlotID','Plot','Year'),na.action=na.omit,all=T)
NLME.LAI4=merge(NLME.LAI3,MAG.NLME.LAI,by=c('PlotID','Plot','Year'),na.action=na.omit,all=T)
NLME.LAI5=merge(NLME.LAI4,DCL.NLME.LAI,by=c('PlotID','Plot','Year'),na.action=na.omit,all=T)
NLME.LAI6=merge(NLME.LAI5,SMAGMK.NLME.LAI,by=c('PlotID','Plot','Year'),na.action=na.omit,all=T)
NLME.LAIall=merge(NLME.LAI6,SAPMK.NLME.LAI,by=c('PlotID','Plot','Year'),na.action=na.omit,all=T)
NLME=merge(NLME.LAIall,LAI,by=c('PlotID','Plot','Year'),na.action=na.omit,all=T)
# Calculates AAD
NLME$SAP.Dif<-abs(NLME$SAP.NLME.LAI-NLME$LAI.spline);NLME$SCL.Dif<-abs(NLME$SCL.NLME.LAI-NLME$LAI.spline)
NLME$SMAG.Dif<-abs(NLME$SMAG.NLME.LAI-NLME$LAI.spline);NLME$MAG.Dif<-abs(NLME$MAG.NLME.LAI-NLME$LAI.spline)
NLME$VAL.Dif<-abs(NLME$VAL.NLME.LAI-NLME$LAI.spline);NLME$DCL.Dif<-abs(NLME$DCL.NLME.LAI-NLME$LAI.spline)
NLME$SMAGMK.Dif<-abs(NLME$SMAGMK.NLME.LAI-NLME$LAI.spline);NLME$SAPMK.Dif<-abs(NLME$SAPMK.NLME.LAI-NLME$LAI.spline)
#NLME=subset(NLME,DCL.Dif!="NA")

write.csv(NLS, file = ("C:/Users/Frank/Desktop/R output/AAD.NLS.csv"), row.names = FALSE)
write.csv(WNLS, file = ("C:/Users/Frank/Desktop/R output/AAD.WNLS.csv"), row.names = FALSE)
write.csv(NLME, file = ("C:/Users/Frank/Desktop/R output/AAD.NLME.csv"), row.names = FALSE)

################################################################################################
# Rename Treatments
levels(NLS$Treatment)=c("B","LD","LT","NT","PCT")
levels(WNLS$Treatment)=c("B","LD","LT","NT","PCT")
levels(NLME$Treatment)=c("B","LD","LT","NT","PCT")

# Renames "thinType" to "Treatment" to match DM
names(DM)[names(DM)=="thinType"] = "Treatment";levels(DM$Treatment)=c("B","LD","LT","NT","PCT","SW")
# Adds TOPHT from DM
DM=subset(DM,select=c(Year,PlotID,Stand,Plot,Treatment,SI,TOPHT))
NLS=merge(NLS,DM,by=c("Year","PlotID","Plot","Stand","Treatment"))
WNLS=merge(WNLS,DM,by=c("Year","PlotID","Plot","Stand","Treatment"))
NLME=merge(NLME,DM,by=c("Year","PlotID","Plot","Stand","Treatment"))
library(beanplot)

# Calculate NLS residuals (Obs-Pred)
NLS$SAP.Resid<-(NLS$LAI.spline-NLS$SAP.NLS.LAI);NLS$SCL.Resid<-(NLS$LAI.spline-NLS$SCL.NLS.LAI)
NLS$SMAG.Resid<-(NLS$LAI.spline-NLS$SMAG.NLS.LAI);NLS$MAG.Resid<-(NLS$LAI.spline-NLS$MAG.NLS.LAI)
NLS$VAL.Resid<-(NLS$LAI.spline-NLS$VAL.NLS.LAI);NLS$DCL.Resid<-(NLS$LAI.spline-NLS$DCL.NLS.LAI)
NLS$SMAGMK.Resid<-(NLS$LAI.spline-NLS$SMAGMK.NLS.LAI);NLS$SAPMK.Resid<-(NLS$LAI.spline-NLS$SAPMK.NLS.LAI)
####################################
par(mfrow=c(2,4))
par(mar=c(5,5,4,2))
beanplot(NLS$SAP.Resid~NLS$Treatment,main="SAP",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-6,4),xlab="Treatment")
beanplot(NLS$SCL.Resid~NLS$Treatment,main="SCL",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-6,4),xlab="Treatment")
beanplot(NLS$SMAG.Resid~NLS$Treatment,main="SMAG",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-6,4),xlab="Treatment")
beanplot(NLS$MAG.Resid~NLS$Treatment,main="MAG",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-6,4),xlab="Treatment")
beanplot(NLS$VAL.Resid~NLS$Treatment,main="VAL",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-6,4),xlab="Treatment")
beanplot(NLS$DCL.Resid~NLS$Treatment,main="DCL",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-6,4),xlab="Treatment")
beanplot(NLS$SMAGMK.Resid~NLS$Treatment,main="SMAGMK",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-6,4),xlab="Treatment")
beanplot(NLS$SAPMK.Resid~NLS$Treatment,main="SAP.MK",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-6,4),xlab="Treatment")

plot(NLS$TOPHT,NLS$SAP.Resid,main="SAP",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-4,3),xlab="TOPHT (m)")
plot(NLS$TOPHT,NLS$SCL.Resid,main="SCL",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-4,3),xlab="TOPHT (m)")
plot(NLS$TOPHT,NLS$SMAG.Resid,main="SMAG",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-4,3),xlab="TOPHT (m)")
plot(NLS$TOPHT,NLS$MAG.Resid,main="MAG",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-4,3),xlab="TOPHT (m)")
plot(NLS$TOPHT,NLS$VAL.Resid,main="VAL",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-4,3),xlab="TOPHT (m)")
plot(NLS$TOPHT,NLS$DCL.Resid,main="DCL",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-4,3),xlab="TOPHT (m)")
plot(NLS$TOPHT,NLS$SMAKMK.Resid,main="SMAGMK",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-4,3),xlab="TOPHT (m)")
plot(NLS$TOPHT,NLS$SAPMK.Resid,main="SAP.MK",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-4,3),xlab="TOPHT (m)")
#####################################################################################################

# WNLS
WNLS$SAP.Resid<-(WNLS$LAI.spline-WNLS$SAP.WNLS.LAI);WNLS$SCL.Resid<-(WNLS$LAI.spline-WNLS$SCL.WNLS.LAI)
WNLS$SMAG.Resid<-(WNLS$LAI.spline-WNLS$SMAG.WNLS.LAI);WNLS$MAG.Resid<-(WNLS$LAI.spline-WNLS$MAG.WNLS.LAI)
WNLS$VAL.Resid<-(WNLS$LAI.spline-WNLS$VAL.WNLS.LAI);WNLS$DCL.Resid<-(WNLS$LAI.spline-WNLS$DCL.WNLS.LAI)
WNLS$SMAGMK.Resid<-(WNLS$LAI.spline-WNLS$SMAGMK.WNLS.LAI);WNLS$SAPMK.Resid<-(WNLS$LAI.spline-WNLS$SAPMK.WNLS.LAI)

par(mfrow=c(2,4))
par(mar=c(5,5,4,2))
beanplot(WNLS$SAP.Resid~WNLS$Treatment,main="SAP",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-4,3),xlab="Treatment")
beanplot(WNLS$SCL.Resid~WNLS$Treatment,main="SCL",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-4,3),xlab="Treatment")
beanplot(WNLS$SMAG.Resid~WNLS$Treatment,main="SMAG",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-4,3),xlab="Treatment")
beanplot(WNLS$MAG.Resid~WNLS$Treatment,main="MAG",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-4,3),xlab="Treatment")
beanplot(WNLS$VAL.Resid~WNLS$Treatment,main="VAL",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-4,3),xlab="Treatment")
beanplot(WNLS$DCL.Resid~WNLS$Treatment,main="DCL",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-4,3),xlab="Treatment")
beanplot(WNLS$SAPMK.Resid~WNLS$Treatment,main="SAP.MK",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-4,3),xlab="Treatment")
beanplot(WNLS$SMAGMK.Resid~WNLS$Treatment,main="SMAGMK",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-4,3),xlab="Treatment")


plot(WNLS$TOPHT,WNLS$SAP.Resid,main="SAP",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-8,4),xlab="TOPHT (m)")
plot(WNLS$TOPHT,WNLS$SCL.Resid,main="SCL",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-8,4),xlab="TOPHT (m)")
plot(WNLS$TOPHT,WNLS$SMAG.Resid,main="SMAG",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-8,4),xlab="TOPHT (m)")
plot(WNLS$TOPHT,WNLS$MAG.Resid,main="MAG",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-8,4),xlab="TOPHT (m)")
plot(WNLS$TOPHT,WNLS$VAL.Resid,main="VAL",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-8,4),xlab="TOPHT (m)")
plot(WNLS$TOPHT,WNLS$DCL.Resid,main="DCL",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-8,4),xlab="TOPHT (m)")
plot(WNLS$TOPHT,WNLS$SAPMK.Resid,main="SAP.MK",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-8,4),xlab="TOPHT (m)")
plot(WNLS$TOPHT,WNLS$SMAGMK.Resid,main="SMAG.MK",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-8,4),xlab="TOPHT (m)")

# NLME
NLME$SAP.Resid<-(NLS$LAI.spline-NLS$SAP.NLS.LAI);NLME$SCL.Resid<-(NLS$LAI.spline-NLS$SCL.NLS.LAI)
NLME$SMAG.Resid<-(NLS$LAI.spline-NLS$SMAG.NLS.LAI);NLME$MAG.Resid<-(NLS$LAI.spline-NLS$MAG.NLS.LAI)
NLME$VAL.Resid<-(NLS$LAI.spline-NLS$VAL.NLS.LAI);NLME$DCL.Resid<-(NLS$LAI.spline-NLS$DCL.NLS.LAI)
NLME$SMAGMK.Resid<-(NLS$LAI.spline-NLS$SMAGMK.NLS.LAI);NLME$SAPMK.Resid<-(NLS$LAI.spline-NLS$SAPMK.NLS.LAI)

par(mfrow=c(2,4))
par(mar=c(5,5,4,2))
beanplot(NLME$SAP.Resid~NLME$Treatment,main="SAP",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-4,4),xlab="Treatment")
beanplot(NLME$SCL.Resid~NLME$Treatment,main="SCL",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-4,4),xlab="Treatment")
beanplot(NLME$SMAG.Resid~NLME$Treatment,main="SMAG",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-4,4),xlab="Treatment")
beanplot(NLME$MAG.Resid~NLME$Treatment,main="MAG",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-4,4),xlab="Treatment")
beanplot(NLME$VAL.Resid~NLME$Treatment,main="VAL",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-4,4),xlab="Treatment")
beanplot(NLME$DCL.Resid~NLME$Treatment,main="DCL",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-4,4),xlab="Treatment")
beanplot(NLME$SAPMK.Resid~NLME$Treatment,main="SAP.MK",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-4,4),xlab="Treatment")
beanplot(NLME$SMAGMK.Resid~NLME$Treatment,main="SMAGMK",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-4,4),xlab="Treatment")

plot(NLME$TOPHT,NLME$SAP.Resid,main="SAP",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-8,4),xlab="TOPHT (m)")
plot(NLME$TOPHT,NLME$SCL.Resid,main="SCL",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-8,4),xlab="TOPHT (m)")
plot(NLME$TOPHT,NLME$SMAG.Resid,main="SMAG",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-8,4),xlab="TOPHT (m)")
plot(NLME$TOPHT,NLME$MAG.Resid,main="MAG",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-8,4),xlab="TOPHT (m)")
plot(NLME$TOPHT,NLME$VAL.Resid,main="VAL",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-8,4),xlab="TOPHT (m)")
plot(NLME$TOPHT,NLME$DCL.Resid,main="DCL",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-8,4),xlab="TOPHT (m)")
plot(NLME$TOPHT,NLME$SAPMK.Resid,main="SAP.MK",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-8,4),xlab="TOPHT (m)")
plot(NLME$TOPHT,NLME$SMAGMK.Resid,main="SMAG.MK",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-8,4),xlab="TOPHT (m)")

###########
# Subset WNLS to 
head(WNLS)
WNLS2008=subset(WNLS,Year=="2008")
# 2008 has no litterfall data 
WNLS2008

WNLS1=subset(WNLS,SCL.WNLS.LAI!="NA")
WNLS2=subset(WNLS1,LAI.spline!="NA")
WNLS2=subset(WNLS2,,select=c(Year,Plot,PlotID,LAI.spline,VAL.WNLS.LAI))
WNLS2
write.csv(WNLS2, file = ("C:/Users/Frank/Desktop/TreePLALAI.csv"), row.names = FALSE)
str(WNLS2)
par(mar=c(5,5,4,2))
plot(WNLS2$VAL.WNLS.LAI~WNLS2$LAI.spline,pch=3,col="red",ylab=expression(paste("Predicted LAI (", m^2,m^-2 ,")"),sep=""),xlab=expression(paste("Litterfall LAI (", m^2,m^-2 ,")"),sep=""),xlim=c(0,8),ylim=c(0,9))
points(WNLS2$VAL.WNLS.LAI~WNLS2$LAI.spline,pch=3,col="red")
# Adds a nonlinear trendline
lines(lowess(WNLS2$LAI.spline,WNLS2$VAL.WNLS.LAI),lty=3);abline(0,1)


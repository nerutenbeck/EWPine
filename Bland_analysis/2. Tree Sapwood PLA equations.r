library(RODBC);library(nlme);library(car)
temp=odbcConnectAccess("C:/Users/Frank/Desktop/Bland WhitePine Leaf Area Database.mdb")
TreeBLA=sqlFetch(temp,"Bland Branch PLA model_Tree LAs");close(temp)
'Location'=as.factor('Location');'Thinned'=as.factor('Thinned');'Tree'=as.factor('Tree')
names(TreeBLA)[names(TreeBLA)=="PLA"] = "TLA";names(TreeBLA)[names(TreeBLA)=="HT"] = "Ht";names(TreeBLA)[names(TreeBLA)=="Tree ID"] = "Tree"
# Tree GEA
head(TreeBLA)
# Calculates DBH/Ht ratioo
TreeBLA$DBH.Ht=TreeBLA$DBH/TreeBLA$Ht
# Calculates CL/SBA ratio
TreeBLA$CL.SBA=TreeBLA$CL/TreeBLA$SBA
# Calculates BA*mLCR
TreeBLA$BA.mLCR=TreeBLA$BAcm/TreeBLA$mLCR
##############  DIagnostics  #####################
# Are Barker's trees outliers due to inaccurate SLA measurements?
# Sapwood area R2=0.9292, R2adj=0.928 
plot(TLA~SBA,data=TreeBLA,pch=16,ylab=expression(paste("PLA(",m^2,")",sep="")),xlab=expression(paste("SBA(",cm^2,")",sep="")))
SBA.lm=lm(TLA~SBA,data=TreeBLA)
summary(SBA.lm)
# Plots PLA over CL
plot(TLA~CL, data=TreeBLA,ylab=expression(paste("PLA(",m^2,")",sep="")),xlab="CL(m)",pch=16)
plot(TLA~mLCR, data=TreeBLA)

# Basal area - R2=0.9439, R2adj=0.943
plot(TLA~BAcm, data=TreeBLA)
BA.lm=lm(TLA~BAcm,data=TreeBLA)
summary(BA.lm)
# DBH
plot(TLA~DBH,data=TreeBLA,ylab=expression(paste("PLA(",m^2,")",sep="")),xlab="DBH(cm)",pch=20)
DBH.lm=lm(TLA~DBH,data=TreeBLA)
summary(DBH.lm)

# DBH:Ht ratio
plot(TLA~DBH.Ht,data=TreeBLA,ylab=expression(paste("PLA(",m^2,")",sep="")),xlab="DBH:Ht",pch=16)
head(TreeBLA)

# Plots BA*mLCR
plot(TLA~BA.mLCR,data=TreeBLA,ylab=expression(paste("PLA(",m^2,")",sep="")),xlab="BA*mLCR",pch=16)
BAmLCR.lm=lm(TLA~BA.mLCR,data=TreeBLA)
summary(BAmLCR.lm)
# Plots CL/SBA
plot(TLA~CL.SBA,data=TreeBLA,ylab=expression(paste("PLA(",m^2,")",sep="")),xlab="CL/SBA",pch=16)

######################################################################################################
# Summary Statistics for TreeBLA
head(TreeBLA)
mean(TreeBLA$TLA)
sd((TreeBLA$TLA)/sqrt(61))
min(TreeBLA$TLA)
max(TreeBLA$TLA)
#######################################################################################################################
# AIC= 577.487
SAPNLS=gnls((TLA)~(a1*(SBA^a2)),data=TreeBLA,na.action=na.omit,start=c(a1=0.8,a2=1.26),control=nlmeControl(minScale=1e-10,returnObject=T))
summary(SAPNLS)
plot(SAPNLS)
# AIC = 491.3681
SAPWNLS=gnls((TLA)~(a1*(SBA^a2)),data=TreeBLA,na.action=na.omit,weights=varPower(0.5,form="SBA"),start=c(a1=0.07,a2=1.3),control=nlmeControl(minScale=1e-10,returnObject=T))
summary(SAPWNLS)  
plot(SAPWNLS)
#AIC= 469.8445 3 warning messages Error step halving factor reduced below minimum IN pnls step
SAPNLMER=nlme((TLA)~(a1*(SBA^a2)),data=TreeBLA,na.action=na.omit,weights=varPower(0.5,form="SBA"),verbose=T,fixed=a1+a2~1, random=a2~1|Location/Thinned,start=c(a1=9.7,a2=0.4),control=nlmeControl(minScale=1e-10,returnObject=T))
summary(SAPNLMER)
plot(SAPNLMER)
ranef(SAPNLMER)
fixef(SAPNLMER)

##########################################################################################################################
# Maguire and Bennett (1996)- AIC= 539.9221
MAGNLS=gnls((TLA)~(a1*CL^a2)*exp(a3*(DBH/Ht)),data=TreeBLA, na.action=na.omit,start=c(a1=4.197871,a2=0.232025,a3=1.043624),control=nlmeControl(minScale=1e-10,returnObject=T)) 
summary(MAGNLS)
# AIC=469.7052
MAGWNLS=gnls((TLA)~(a1*CL^a2)*exp(a3*(DBH/Ht)),data=TreeBLA,na.action=na.omit,weights=varPower(0.5,form="CL"),start=c(a1=4.197871,a2=0.232025,a3=1.043624),control=nlmeControl(minScale=1e-10,returnObject=T)) 
summary(MAGWNLS)
# AIC= 453.5994 a1~1
MAGNLMER=nlme((TLA)~(a1*CL^a2)*exp(a3*(DBH/Ht)),data=TreeBLA,na.action=na.omit,weights=varPower(0.5,form="CL"),fixed=a1+a2+a3~1,random=a1~1|Location/Thinned,start=c(a1=0.1,a2=2.9,a3=0.5),control=nlmeControl(minScale=1e-10,returnObject=T)) 
plot(MAGNLMER)
plot(TreeBLA$Ht,resid(MAGNLMER))
summary(MAGNLMER)
fixef(MAGNLMER)
ranef(MAGNLMER)
###########################################################################################################################
#AIC 572.8211
SMAGNLS=gnls((TLA)~(a1*SBA^a2)*exp(a3*(DBH/Ht)),data=TreeBLA, na.action=na.omit,start=c(a1=0.1,a2=1.3,a3=-0.1),control=nlmeControl(minScale=1e-10,returnObject=T)) 
summary(SMAGNLS)  
#AIC= 445.6874
SMAGWNLS=gnls((TLA)~(a1*SBA^a2)*exp(a3*(DBH/Ht)),data=TreeBLA, weights=varPower(0.5,form="SBA"),start=c(a1=0.6,a2=1,a3=0.6),na.action=na.omit,control=nlmeControl(minScale=1e-10,returnObject=T)) 
summary(SMAGWNLS) 
#SMAG NLME-R- AIC= 435.9749 a3~1, AIC=407.2883-a3~1, AIC= 415.8007 a1~1 
SMAGNLMER=nlme((TLA)~(a1*SBA^a2)*exp(a3*(DBH/Ht)),data=TreeBLA,na.action=na.omit,weights=varPower(0.5,form="SBA"),fixed=a1+a2+a3~1, random=a2~1|Location/Thinned,start=c(a1=0.19,a2=0.8,a3=1.2),control=nlmeControl(minScale=1e-10,returnObject=T)) 
summary(SMAGNLMER) 
ranef(SMAGNLMER)
fixef(SMAGNLMER)
############################################################################################################################
#AIC= 531.2726
BACLNLS=gnls((TLA)~(a1*BAcm^a2*(CL^a3)),data=TreeBLA, na.action=na.omit,start=c(a1=0.3970,a2=0.6149,a3=1.6450),control=nlmeControl(minScale=1e-10,returnObject=T)) 
summary(BACLNLS)
#AIC=404.3337
BACLWNLS=gnls((TLA)~(a1*BAcm^a2*(CL^a3)),data=TreeBLA,na.action=na.omit,weights=varPower(0.5,form="BAcm"),start=c(a1=0.3970,a2=0.6149,a3=1.6450),control=nlmeControl(minScale=1e-10,returnObject=T)) 
summary(BACLWNLS) 
# AIC=408.3337
BACLNLMER=nlme((TLA)~(a1*(BAcm^a2)*(CL^a3)),data=TreeBLA, na.action=na.omit,weights=varPower(0.5,form="BAcm"),fixed=a1+a2+a3~1,random=a3~1|Location/Thinned,start=c(a1=0.06,a2=0.88,a3=1.97),control=nlmeControl(minScale=1e-10,returnObject=T))
summary(BACLNLMER) 
ranef(BACLNLMER)
fixef(BACLNLMER)
################################################################################################################################
#AIC= 531.2726
DCLNLS=gnls((TLA)~(a1*DBH^a2*(CL^a3)),data=TreeBLA, na.action=na.omit,start=c(a1=0.3970,a2=0.6149,a3=1.6450),control=nlmeControl(minScale=1e-10,returnObject=T)) 
summary(DCLNLS)
plot(DCLNLS)
#AIC=  404.3339
DCLWNLS=gnls((TLA)~(a1*DBH^a2*(CL^a3)),data=TreeBLA,na.action=na.omit,weights=varPower(0.5,form="DBH"),start=c(a1=0.3970,a2=0.6149,a3=1.6450),control=nlmeControl(minScale=1e-10,returnObject=T)) 
summary(DCLWNLS) 
plot(DCLWNLS)
# AIC=407.5144 a2~1 OR a3~1, AIC=387.476 when a1~1 
DCLNLMER=nlme((TLA)~(a1*(DBH^a2)*(CL^a3)),data=TreeBLA, na.action=na.omit,weights=varPower(0.5,form="DBH"),fixed=a1+a2+a3~1,random=a3~1|Location/Thinned,start=c(a1=0.06,a2=0.88,a3=1.97),control=nlmeControl(minScale=1e-10,returnObject=T))
summary(DCLNLMER) 
ranef(DCLNLMER)
fixef(DCLNLMER)
plot(DCLNLMER)
###########################################################################################################################
# AIC= 531.2884
SCLNLS=gnls((TLA)~(a1*SBA^a2*CL^a3),data=TreeBLA,na.action=na.omit,start=c(a1=0.15,a2=0.5,a3=1.6),control=nlmeControl(minScale=1e-10,returnObject=T)) 
summary(SCLNLS)  
plot(SCLNLS)
# WNLS AIC= 421.437- underpredicts leaf area
SCLWNLS=gnls((TLA)~(a1*SBA^a2*CL^a3),data=TreeBLA,na.action=na.omit,weights=varPower(0.5,form="SBA"),start=c(a1=0.1383,a2=1.1947,a3=1),control=nlmeControl(minScale=1e-10,returnObject=T)) 
summary(SCLWNLS)
plot(SCLWNLS)
# AIC= 419.5364 when a3~1, AIC= 390.396 when a2~1, AIC=391.4981 when a1~1, AIC= 
SCLNLMER=nlme((TLA)~(a1*(SBA^a2)*(CL^a3)),data=TreeBLA,na.action=na.omit,weights=varPower(0.5,form="SBA"),        
fixed=a1+a2+a3~1,random=a3~1|Location/Thinned,start=c(a1=0.1,a2=0.8,a3=1.2),control=nlmeControl(minScale=1e-10,returnObject=T)) 
summary(SCLNLMER)  
ranef(SCLNLMER) 
fixef(SCLNLMER) 
plot(SCLNLMER)

###########################################################################################################################
# AIC= 570.5339- overpredicts small trees, underpredicts for mid-sized trees
VALNLS=gnls((TLA)~(a1)*((BAcm*mLCR)^a2),data=TreeBLA, na.action=na.omit,start=c(a1=0.9,a2=0.84),control=nlmeControl(minScale=1e-10,returnObject=T)) 
summary(VALNLS)
plot(VALNLS)
# AIC= 371.2438
VALWNLS=gnls((TLA)~(a1)*((BAcm*mLCR)^a2),data=TreeBLA,na.action=na.omit,weights=varPower(0.5,form="BAcm"),start=c(a1=0.9054,a2=0.8374),control=nlmeControl(minScale=1e-10,returnObject=T)) 
summary(VALWNLS)
plot(VALWNLS)
plot(resid$VALWNLS
#AIC=368.1343 when a2~1,AIC=344.0343 when a1~1, 
VALNLMER=nlme((TLA)~(a1*(BAcm*mLCR)^a2),data=TreeBLA,na.action=na.omit,weights=varPower(0.5,form="BAcm"),      
fixed=a1+a2~1,random=a2~1|Location/Thinned,start=c(a1=0.01,a2=1.6),control=nlmeControl(minScale=1e-10,returnObject=T)) 
summary(VALNLMER)
ranef(VALNLMER)
fixef(VALNLMER)
#############  SAP.MK  ###############################3
# AIC= 531.2884
SAPNLS.MK=gnls(((TLA)~(a1*(SBA^a2)*(CL/SBA)^a3)),data=TreeBLA,na.action=na.omit,start=c(a1=0.3,a2=.1,a3=2),control=nlmeControl(minScale=1e-10,returnObject=T))
summary(SAPNLS.MK)
#AIC=421.437
SAPWNLS.MK=gnls(((TLA)~(a1*(SBA^a2)*(CL/SBA)^a3)),data=TreeBLA,na.action=na.omit,weights=varPower(0.5,form="SBA"),start=c(a1=0.3,a2=.1,a3=2),control=nlmeControl(minScale=1e-10,returnObject=T))
summary(SAPWNLS.MK)
# AIC 420.7144- 
SAPNLME.MK=nlme(((TLA)~(a1*(SBA^a2)*(CL/SBA)^a3)),data=TreeBLA,na.action=na.omit,weights=varPower(0.5,form="SBA"),fixed=a1+a2+a3~1,random=a2~1|Location/Thinned,start=c(a1=0.3,a2=.1,a3=2),control=nlmeControl(minScale=1e-10,returnObject=T))
summary(SAPNLME.MK)
fixef(SAPNLME.MK)
ranef(SAPNLME.MK)
##########################  SMAG.MK  #############################################################
# AIC=420.7144
SMAGNLS.MK=gnls(((TLA)~(a1*(SBA^a2)*(CL/SBA)^a3)+((DBH/Ht)^a4)),data=TreeBLA,na.action=na.omit,
start=c(a1=0.14,a2=2.13,a3=1.6,a4=.13),control=nlmeControl(minScale=1e-10,returnObject=T))
summary(SMAGNLS.MK)
# AIC=  404.5328
SMAGWNLS.MK=gnls(((TLA)~(a1*(SBA^a2)*(CL/SBA)^a3)+((DBH/Ht)^a4)),data=TreeBLA,na.action=na.omit,weights=varPower(0.5,form="SBA"),
start=c(a1=0.04,a2=2.2,a3=1.4,a4=3),control=nlmeControl(minScale=1e-10,returnObject=T))
summary(SMAGWNLS.MK)
# AIC= 408.5327- 
SMAGNLME.MK=nlme(((TLA)~(a1*(SBA^a2)*(CL/SBA)^a3)+((DBH/Ht)^a4)),data=TreeBLA,na.action=na.omit,weights=varPower(0.5,form="SBA"),
fixed=a1+a2+a3+a4~1,random=a2~1|Location/Thinned,start=c(a1=0.3,a2=.1,a3=2,a4=3),control=nlmeControl(minScale=1e-10,returnObject=T))
summary(SMAGNLME.MK)
ranef(SMAGNLME.MK)
fixef(SMAGNLME.MK)
########################################################################################################################
# Predict Tree leaf area using parameter estimates generated above
#################################################################################
# Reads in parameter estimates for NLME models for fixed and random effects
SAP<-read.csv("C:/Users/Frank/Desktop/CSV Files/SAP parameter estimates.csv");MAG<-read.csv("C:/Users/Frank/Desktop/CSV Files/MAG parameter estimates.csv")
SMAG<-read.csv("C:/Users/Frank/Desktop/CSV Files/SMAG parameter estimates.csv");DCL<-read.csv("C:/Users/Frank/Desktop/CSV Files/DCL parameter estimates.csv")
SCL<-read.csv("C:/Users/Frank/Desktop/CSV Files/SCL parameter estimates.csv");VAL<-read.csv("C:/Users/Frank/Desktop/CSV Files/VAL parameter estimates.csv")
SAP.MK<-read.csv("C:/Users/Frank/Desktop/CSV Files/SAPNLME.MK parameter estimates.csv");SMAG.MK<-read.csv("C:/Users/Frank/Desktop/CSV Files/SMAGNLME.MK parameter estimates.csv")

# Merges model parameter estimates with Tree attribute data
SAP.=merge(SAP,TreeBLA,by='Tree',all=T);MAG.=merge(MAG,TreeBLA,by='Tree',all=T);SMAG.=merge(SMAG,TreeBLA,by='Tree',all=T);DCL.=merge(DCL,TreeBLA,by='Tree',all=T)
SCL.=merge(SCL,TreeBLA,by='Tree',all=T);VAL.=merge(VAL,TreeBLA,by='Tree',all=T);SAP.MK.=merge(SAP.MK,TreeBLA,by='Tree',all=T);SMAG.MK.=merge(SMAG.MK,TreeBLA,by='Tree',all=T)

# Predicts Tree leaf area for NLME-F and NLME-R
SAP.$TLASAPNLMER<-SAP.$a1*SAP.$SBA^SAP.$a2;SAP.$TLASAPNLMEF<-SAP.$a1*SAP.$SBA^SAP.$a2_fixed
MAG.$TLAMAGNLMER<-MAG.$a1*(MAG.$CL^MAG.$a2)*exp(MAG.$a3*(MAG.$DBH/MAG.$Ht));MAG.$TLAMAGNLMEF<-MAG.$a1_fixed*(MAG.$CL^MAG.$a2)*exp(MAG.$a3*(MAG.$DBH/MAG.$Ht))
SMAG.$TLASMAGNLMER<-SMAG.$a1*(SMAG.$SBA^SMAG.$a2)*exp(SMAG.$a3*(SMAG.$DBH/SMAG.$Ht));SMAG.$TLASMAGNLMEF<-SMAG.$a1*(SMAG.$SBA^SMAG.$a2_fixed)*exp(SMAG.$a3*(SMAG.$DBH/SMAG.$Ht))
DCL.$TLADCLNLMER<-DCL.$a1*(DCL.$DBH^DCL.$a2)*(DCL.$CL^DCL.$a3);DCL.$TLADCLNLMEF<-DCL.$a1*(DCL.$DBH^DCL.$a2)*(DCL.$CL^DCL.$a3_fixed)
SCL.$TLASCLNLMER<-SCL.$a1*(SCL.$SBA^SCL.$a2)*(SCL.$CL^SCL.$a3);SCL.$TLASCLNLMEF<-SCL.$a1*(SCL.$SBA^SCL.$a2)*(SCL.$CL^SCL.$a3_fixed)
VAL.$TLAVALNLMER<-VAL.$a1*((VAL.$BAcm*VAL.$mLCR)^VAL.$a2);VAL.$TLAVALNLMEF<-VAL.$a1*((VAL.$BAcm*VAL.$mLCR)^VAL.$a2_fixed)
SAP.MK.$TLASAPMKNLMER<-SAP.MK.$a1*(SAP.MK.$SBA^SAP.MK.$a2)*((SAP.MK.$CL/SAP.MK.$SBA)^SAP.MK.$a3);SAP.MK.$TLASAPMKNLMEF<-SAP.MK.$a1*(SAP.MK.$SBA^SAP.MK.$a2_fixed)*((SAP.MK.$CL/SAP.MK.$SBA)^SAP.MK.$a3)
SMAG.MK.$TLASMAGMKNLMER<-(SMAG.MK.$a1*(SMAG.MK.$SBA^SMAG.MK.$a2))*((SMAG.MK.$CL/SMAG.MK.$SBA)^SMAG.MK.$a3)+((SMAG.MK.$DBH/SMAG.MK.$Ht)^SMAG.MK.$a4);SMAG.MK.$TLASMAGMKNLMEF<-(SMAG.MK.$a1*(SMAG.MK.$SBA^SMAG.MK.$a2_fixed))*((SMAG.MK.$CL/SMAG.MK.$SBA)^SMAG.MK.$a3)+((SMAG.MK.$DBH/SMAG.MK.$Ht)^SMAG.MK.$a4)
# Removes NAs caused by climbed trees
SAP.=subset(SAP.,!is.na(TLASAPNLMEF),select=,)
MAG.=subset(MAG.,!is.na(TLAMAGNLMEF),select=,)
SMAG.=subset(SMAG.,!is.na(TLASMAGNLMEF),select=,)
DCL.=subset(DCL.,!is.na(TLADCLNLMEF),select=,)
SCL.=subset(SCL.,!is.na(TLASCLNLMEF),select=,)
VAL.=subset(VAL.,!is.na(TLAVALNLMEF),select=,)
SAP.MK.=subset(SAP.MK.,!is.na(TLASAPMKNLMEF),select=,)
SMAG.MK.=subset(SMAG.MK.,!is.na(TLASMAGMKNLMEF),select=,)

# Predict Tree leaf area for NLS and WNLS model    
TreeBLA$TLASAPNLS=(predict(SAPNLS,newdata=TreeBLA,na.action=na.omit,levels=0));TreeBLA$TLASAPWNLS=(predict(SAPWNLS,newdata=TreeBLA,na.action=na.omit,levels=0))
TreeBLA$TLASCLNLS=(predict(SCLNLS,newdata=TreeBLA,na.action=na.omit,levels=0));TreeBLA$TLASCLWNLS=(predict(SCLWNLS,newdata=TreeBLA,na.action=na.omit,levels=0))
TreeBLA$TLAMAGNLS=(predict(MAGNLS,newdata=TreeBLA,na.action=na.omit,levels=0));TreeBLA$TLAMAGWNLS=(predict(MAGWNLS,newdata=TreeBLA,na.action=na.omit,levels=0))
TreeBLA$TLASMAGNLS=(predict(SMAGNLS,newdata=TreeBLA,na.action=na.omit,levels=0));TreeBLA$TLASMAGWNLS=(predict(SMAGWNLS,newdata=TreeBLA,na.action=na.omit,levels=0))
TreeBLA$TLADCLNLS=(predict(DCLNLS,newdata=TreeBLA,na.action=na.omit,levels=0));TreeBLA$TLADCLWNLS=(predict(DCLWNLS,newdata=TreeBLA,na.action=na.omit,levels=0))
TreeBLA$TLAVALNLS=(predict(VALNLS,newdata=TreeBLA,na.action=na.omit,levels=0));TreeBLA$TLAVALWNLS=(predict(VALWNLS,newdata=TreeBLA,na.action=na.omit,levels=0))
TreeBLA$TLAVALNLS=(predict(VALNLS,newdata=TreeBLA,na.action=na.omit,levels=0));TreeBLA$TLAVALWNLS=(predict(VALWNLS,newdata=TreeBLA,na.action=na.omit,levels=0))

# New equations
TreeBLA$TLASAPNLS.MK=(predict(SAPNLS.MK,newdata=TreeBLA,levels=0,na.action=na.omit));TreeBLA$TLASAPWNLS.MK=(predict(SAPWNLS.MK,newdata=TreeBLA,levels=0,na.action=na.omit))
TreeBLA$TLASMAGNLS.MK=(predict(SMAGNLS.MK,newdata=TreeBLA,levels=0,na.action=na.omit));TreeBLA$TLASMAGWNLS.MK=(predict(SMAGWNLS.MK,newdata=TreeBLA,levels=0,na.action=na.omit))

############## Plot Residuals over Ht #################################################
par(mfrow=c(2,4))
     
SCLWNLS.res = resid(SCLWNLS) 
plot(TreeBLA$Ht,SCLWNLS.res) # overpredicts Tree 120SH
SCLWNLS.res    

SCLNLMER.res=resid(SAPNLMER)
plot(TreeBLA$Ht,SCLNLMER.res)
         
DCLWNLS.res=resid(DCLWNLS) # Overpredicts Tree 120SH
plot(TreeBLA$Ht,DCLWNLS.res) 

DCLNLMER.res=resid(DCLNLMER)
plot(TreeBLA$Ht,DCLNLMER.res)
     
DCLWNLS.res=resid(DCLWNLS) 
plot(TreeBLA$Ht,DCLWNLS.res)
     
MAGNLMEF.res=resid(MAGNLMER)
plot(TreeBLA$Ht,MAGNLMEF.res)
################################################################################
##################     RMSE and MAB Calculations     ###########################
################################################################################
rmse<-function(obs,pred) {
sqrt(sum((obs-pred)^2,na.rm=T)/length(obs))
    }
mab<-function(obs,pred) {
    (sum(abs(obs-pred),na.rm=T))/length(obs)
    }
# SAP RMSE and MAB
rmse.SAPNLS<-rmse(obs=TreeBLA$TLA,pred=TreeBLA$TLASAPNLS);rmse.SAPWNLS<-rmse(obs=TreeBLA$TLA,pred=TreeBLA$TLASAPWNLS)
rmse.SAPNLMER<-rmse(obs=TreeBLA$TLA,pred=SAP.$TLASAPNLMER);rmse.SAPNLMEF<-rmse(obs=TreeBLA$TLA,pred=SAP.$TLASAPNLMEF)
MAB.SAPNLS<-mab(obs=TreeBLA$TLA,pred=TreeBLA$TLASAPNLS);MAB.SAPWNLS<-mab(obs=TreeBLA$TLA,pred=TreeBLA$TLASAPWNLS)
MAB.SAPNLMER<-mab(obs=SAP.$TLA,pred=SAP.$TLASAPNLMER);MAB.SAPNLMEF<-mab(obs=SAP.$TLA,pred=SAP.$TLASAPNLMEF)
rmse.SAPNLS;rmse.SAPWNLS;rmse.SAPNLMER;rmse.SAPNLMEF
MAB.SAPNLS;MAB.SAPWNLS;MAB.SAPNLMER;MAB.SAPNLMEF
# SCL RMSE and MAB
rmse.SCLNLS<-rmse(obs=TreeBLA$TLA,pred=TreeBLA$TLASCLNLS);rmse.SCLWNLS<-rmse(obs=TreeBLA$TLA,pred=TreeBLA$TLASCLWNLS)
rmse.SCLNLMER<-rmse(obs=SCL.$TLA,pred=SCL.$TLASCLNLMER);rmse.SCLNLMEF<-rmse(obs=SCL.$TLA,pred=SCL.$TLASCLNLMEF)
MAB.SCLNLS<-mab(obs=TreeBLA$TLA,pred=TreeBLA$TLASCLNLS);MAB.SCLWNLS<-mab(obs=TreeBLA$TLA,pred=TreeBLA$TLASCLWNLS)
MAB.SCLNLMER<-mab(obs=SCL.$TLA,pred=SCL.$TLASCLNLMER);MAB.SCLNLMEF<-mab(obs=SCL.$TLA,pred=SCL.$TLASCLNLMEF)
rmse.SCLNLS;rmse.SCLWNLS;rmse.SCLNLMER;rmse.SCLNLMEF
MAB.SCLNLS;MAB.SCLWNLS;MAB.SCLNLMER;MAB.SCLNLMEF
# MAG RMSE and MAB
rmse.MAGNLS<-rmse(obs=TreeBLA$TLA,pred=TreeBLA$TLAMAGNLS);rmse.MAGWNLS<-rmse(obs=TreeBLA$TLA,pred=TreeBLA$TLAMAGWNLS)
rmse.MAGNLMER<-rmse(obs=MAG.$TLA,pred=MAG.$TLAMAGNLMER);rmse.MAGNLMEF<-rmse(obs=MAG.$TLA,pred=MAG.$TLAMAGNLMEF)
MAB.MAGNLS<-mab(obs=TreeBLA$TLA,pred=TreeBLA$TLAMAGNLS);MAB.MAGWNLS<-mab(obs=TreeBLA$TLA,pred=TreeBLA$TLAMAGWNLS)
MAB.MAGNLMER<-mab(obs=MAG.$TLA,pred=MAG.$TLAMAGNLMER);MAB.MAGNLMEF<-mab(obs=MAG.$TLA,pred=MAG.$TLAMAGNLMEF)
rmse.MAGNLS;rmse.MAGWNLS;rmse.MAGNLMER;rmse.MAGNLMEF
MAB.MAGNLS;MAB.MAGWNLS;MAB.MAGNLMER;MAB.MAGNLMEF
# SMAG RMSE and MAB
rmse.SMAGNLS<-rmse(obs=TreeBLA$TLA,pred=TreeBLA$TLASMAGNLS);rmse.SMAGWNLS<-rmse(obs=TreeBLA$TLA,pred=TreeBLA$TLASMAGWNLS)
rmse.SMAGNLMER<-rmse(obs=SMAG.$TLA,pred=SMAG.$TLASMAGNLMER);rmse.SMAGNLMEF<-rmse(obs=SMAG.$TLA,pred=SMAG.$TLASMAGNLMEF)
MAB.SMAGNLS<-mab(obs=TreeBLA$TLA,pred=TreeBLA$TLASMAGNLS);MAB.SMAGWNLS<-mab(obs=TreeBLA$TLA,pred=TreeBLA$TLASMAGWNLS)
MAB.SMAGNLMER<-mab(obs=SMAG.$TLA,pred=SMAG.$TLASMAGNLMER);MAB.SMAGNLMEF<-mab(obs=SMAG.$TLA,pred=SMAG.$TLASMAGNLMEF)
rmse.SMAGNLS;rmse.SMAGWNLS;rmse.SMAGNLMER;rmse.SMAGNLMEF
MAB.SMAGNLS;MAB.SMAGWNLS;MAB.SMAGNLMER;MAB.SMAGNLMEF
# DCL RMSE and MAB
rmse.DCLNLS<-rmse(obs=TreeBLA$TLA,pred=TreeBLA$TLADCLNLS);rmse.DCLWNLS<-rmse(obs=TreeBLA$TLA,pred=TreeBLA$TLADCLWNLS)
rmse.DCLNLMER<-rmse(obs=DCL.$TLA,pred=DCL.$TLADCLNLMER);rmse.DCLNLMEF<-rmse(obs=DCL.$TLA,pred=DCL.$TLADCLNLMEF)
MAB.DCLNLS<-mab(obs=TreeBLA$TLA,pred=TreeBLA$TLADCLNLS);MAB.DCLWNLS<-mab(obs=TreeBLA$TLA,pred=TreeBLA$TLADCLWNLS)
MAB.DCLNLMER<-mab(obs=DCL.$TLA,pred=DCL.$TLADCLNLMER);MAB.DCLNLMEF<-mab(obs=DCL.$TLA,pred=DCL.$TLADCLNLMEF)
rmse.DCLNLS;rmse.DCLWNLS;rmse.DCLNLMER;rmse.DCLNLMEF
MAB.DCLNLS;MAB.DCLWNLS;MAB.DCLNLMER;MAB.DCLNLMEF
# VAL RMSE and MAB
rmse.VALNLS<-rmse(obs=TreeBLA$TLA,pred=TreeBLA$TLAVALNLS);rmse.VALWNLS<-rmse(obs=TreeBLA$TLA,pred=TreeBLA$TLAVALWNLS)
rmse.VALNLMER<-rmse(obs=VAL.$TLA,pred=VAL.$TLAVALNLMER);rmse.VALNLMEF<-rmse(obs=VAL.$TLA,pred=VAL.$TLAVALNLMEF)
MAB.VALNLS<-mab(obs=TreeBLA$TLA,pred=TreeBLA$TLAVALNLS);MAB.VALWNLS<-mab(obs=TreeBLA$TLA,pred=TreeBLA$TLAVALWNLS)
MAB.VALNLMER<-mab(obs=VAL.$TLA,pred=VAL.$TLAVALNLMER);MAB.VALNLMEF<-mab(obs=VAL.$TLA,pred=VAL.$TLAVALNLMEF)
rmse.VALNLS;rmse.VALWNLS;rmse.VALNLMER;rmse.VALNLMEF
MAB.VALNLS;MAB.VALWNLS;MAB.VALNLMER;MAB.VALNLMEF
# SAP.MK RMSE and MAB
rmse.SAPNLS.MK<-rmse(obs=TreeBLA$TLA,pred=TreeBLA$TLASAPNLS.MK);rmse.SAPWNLS.MK<-rmse(obs=TreeBLA$TLA,pred=TreeBLA$TLASAPWNLS.MK)
rmse.SAPNLMER.MK<-rmse(obs=TreeBLA$TLA,pred=SAP.MK.$TLASAPMKNLMER);rmse.SAPNLMEF.MK<-rmse(obs=TreeBLA$TLA,pred=SAP.MK.$TLASAPMKNLMEF)
MAB.SAPNLS.MK<-mab(obs=TreeBLA$TLA,pred=TreeBLA$TLASAPNLS.MK);MAB.SAPWNLS.MK<-mab(obs=TreeBLA$TLA,pred=TreeBLA$TLASAPWNLS.MK)
MAB.SAPNLMER.MK<-mab(obs=TreeBLA$TLA,pred=SAP.MK.$TLASAPMKNLMER);MAB.SAPNLMEF.MK<-mab(obs=TreeBLA$TLA,pred=SAP.MK.$TLASAPMKNLMEF)
rmse.SAPNLS.MK;rmse.SAPWNLS.MK;rmse.SAPNLMER.MK;rmse.SAPNLMEF.MK
MAB.SAPNLS.MK;MAB.SAPWNLS.MK;MAB.SAPNLMER.MK;MAB.SAPNLMEF.MK    
# SMAG.MK RMSE and MAB
rmse.SMAGNLS.MK<-rmse(obs=TreeBLA$TLA,pred=TreeBLA$TLASMAGNLS.MK);rmse.SMAGWNLS.MK<-rmse(obs=TreeBLA$TLA,pred=TreeBLA$TLASMAGWNLS.MK)
rmse.SMAGNLMER.MK<-rmse(obs=TreeBLA$TLA,pred=SMAG.MK.$TLASMAGMKNLMER);rmse.SMAGNLMEF.MK<-rmse(obs=TreeBLA$TLA,pred=SMAG.MK.$TLASMAGMKNLMEF)
MAB.SMAGNLS.MK<-mab(obs=TreeBLA$TLA,pred=TreeBLA$TLASMAGNLS.MK);MAB.SMAGWNLS.MK<-mab(obs=TreeBLA$TLA,pred=TreeBLA$TLASMAGWNLS.MK)
MAB.SMAGNLMER.MK<-mab(obs=TreeBLA$TLA,pred=SMAG.MK.$TLASMAGMKNLMER);MAB.SMAGNLMEF.MK<-mab(obs=TreeBLA$TLA,pred=SMAG.MK.$TLASMAGMKNLMEF)
rmse.SMAGNLS.MK;rmse.SMAGWNLS.MK;rmse.SMAGNLMER.MK;rmse.SMAGNLMEF.MK
MAB.SMAGNLS.MK;MAB.SMAGWNLS.MK;MAB.SMAGNLMER.MK;MAB.SMAGNLMEF.MK
###########################################################################################################################
################### Graphical Residual Analysis  #################################
# Calculates residuals for NLME-F and NLME-R
SAPNLMER.resid=SAP.$TLA-SAP.$TLASAPNLMER;SAPNLMEF.resid=SAP.$TLA-SAP.$TLASAPNLMEF;na.omit(SAPNLMER.resid)
SCLNLMER.resid=SCL.$TLA-SCL.$TLASCLNLMER;SCLNLMEF.resid=SCL.$TLA-SCL.$TLASCLNLMEF;na.omit(SCLNLMER.resid)
MAGNLMER.resid=MAG.$TLA-MAG.$TLAMAGNLMER;MAGNLMEF.resid=MAG.$TLA-MAG.$TLAMAGNLMEF;na.omit(MAGNLMER.resid)
SMAGNLMER.resid=SMAG.$TLA-SMAG.$TLASMAGNLMER;SMAGNLMEF.resid=SMAG.$TLA-SMAG.$TLASMAGNLMEF;na.omit(SMAGNLMER.resid)
VALNLMER.resid=VAL.$TLA-VAL.$TLAVALNLMER;VALNLMEF.resid=VAL.$TLA-VAL.$TLAVALNLMEF;na.omit(VALNLMER.resid)
DCLNLMER.resid=DCL.$TLA-DCL.$TLADCLNLMER;DCLNLMEF.resid=DCL.$TLA-DCL.$TLADCLNLMEF;na.omit(DCLNLMER.resid)  
SAP.MKNLMER.resid=SAP.MK.$TLA-SAP.MK.$TLASAPMKNLMER;SAP.MKNLMEF.resid=SAP.$TLA-SAP.MK.$TLASAPMKNLMEF;na.omit(SAP.MKNLMER.resid)
SMAG.MKNLMER.resid=SMAG.MK.$TLA-SMAG.MK.$TLASMAGMKNLMER;SMAG.MKNLMEF.resid=SMAG.MK.$TLA-SMAG.MK.$TLASMAGMKNLMEF;na.omit(SMAG.MKNLMER.resid)  

# "" in axis labels will remove axis titles, cex and bold font in labels 
par(mar=c(5,5,4,2))
par(mfrow=c(2,4))
plot(SAPNLMER.resid~SAP.$TLA,main="SAP",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-100,100),xlim=c(0,400),xlab=expression(paste("Projected Leaf Area(",m^2,")",sep="")))
lines(lowess(SAP.$TLA,SAPNLMER.resid),lty=3);abline(h=0)     
points(SAPNLMEF.resid~SAP.$TLA,pch=16);abline(h=0);lines(lowess(SAP.$TLA,SAPNLMEF.resid))
   
plot(SCLNLMER.resid~SCL.$TLA,main="SCL",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-100,100),xlim=c(0,400),xlab=expression(paste("Projected Leaf Area(",m^2,")",sep="")))
lines(lowess(SCL.$TLA,SCLNLMER.resid),lty=3);abline(h=0)
points(SCLNLMEF.resid~SCL.$TLA,pch=16);lines(lowess(SCL.$TLA,SCLNLMEF.resid))
#SCLNLMER.resid;SCLNLMEF.resid         

plot(MAGNLMER.resid~MAG.$TLA,main="MAG",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-100,100),xlim=c(0,400),xlab=expression(paste("Projected Leaf Area(",m^2,")",sep="")))
lines(lowess(MAG.$TLA,MAGNLMER.resid),lty=3);abline(h=0)
points(MAGNLMEF.resid~MAG.$TLA,pch=16);lines(lowess(MAG.$TLA,MAGNLMEF.resid))

#
plot(SMAGNLMER.resid~SMAG.$TLA,main="SMAG",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-100,100),xlim=c(0,400),xlab=expression(paste("Projected Leaf Area(",m^2,")",sep="")))
lines(lowess(SMAG.$TLA,SMAGNLMER.resid),lty=3);abline(h=0)
points(SMAGNLMEF.resid~SMAG.$TLA,pch=16);lines(lowess(MAG.$TLA,MAGNLMEF.resid))

plot(DCLNLMER.resid~DCL.$TLA,main="DCL",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-100,100),xlim=c(0,400),xlab=expression(paste("Projected Leaf Area(",m^2,")",sep="")))
lines(lowess(DCL.$TLA,DCLNLMER.resid),lty=3);abline(h=0)
points(DCLNLMEF.resid~DCL.$TLA,pch=16);lines(lowess(DCL.$TLA,DCLNLMEF.resid))

plot(VALNLMER.resid~VAL.$TLA,main="VAL",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-100,100),xlim=c(0,400),xlab=expression(paste("Projected Leaf Area(",m^2,")",sep="")))
lines(lowess(VAL.$TLA,VALNLMER.resid),lty=3);abline(h=0)
points(VALNLMEF.resid~VAL.$TLA,pch=16);lines(lowess(VAL.$TLA,VALNLMEF.resid))
# SAP.MK
plot(SAP.MKNLMER.resid~SAP.MK.$TLA,main="SAP.MK",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-100,100),xlim=c(0,400),xlab=expression(paste("Projected Leaf Area(",m^2,")",sep="")))
lines(lowess(SAP.MK.$TLA,SAP.MKNLMER.resid),lty=3);abline(h=0)
points(SAP.MKNLMEF.resid~SAP.$TLA,pch=16);lines(lowess(SAP.MK.$TLA,SAP.MKNLMEF.resid))
     
plot(SMAG.MKNLMER.resid~SMAG.MK.$TLA,main="SMAG.MK",ylab=expression(paste("Residuals(",m^2,")",sep="")),ylim=c(-100,100),xlim=c(0,400),xlab=expression(paste("Projected Leaf Area(",m^2,")",sep="")))
lines(lowess(SMAG.MK.$TLA,SMAG.MKNLMER.resid),lty=3);abline(h=0)
points(SMAG.MKNLMEF.resid~SMAG.MK.$TLA,pch=16);lines(lowess(SMAG.MK.$TLA,SMAG.MKNLMEF.resid))
     
########################################################################################################     
# Resuduals over Ht
######################################################################################
par(mfrow=c(2,4))
plot(SAPNLMER.resid~SAP.$Ht,main="SAP",ylab=expression(paste("Residuals(",m^2,")",sep="")),xlim=c(0,40),xlab="TOPHT (m)")
points(SAPNLMEF.resid~SAP.$Ht,pch=16)
SAPNLMER.resid
SAPNLMEF.resid
plot(SCLNLMER.resid~SCL.$Ht,main="SCL",ylab=expression(paste("Residuals(",m^2,")",sep="")),xlim=c(0,40),xlab="TOPHT (m)")
points(SCLNLMEF.resid~SCL.$Ht,pch=16)
SCLNLMER.resid
SCLNLMEF.resid         
# Overestimates tree w/ Tree=43 PLA=85? resid=-140.8, Tree# 2 resid =-351, Tree 9 resid = 116
plot(MAGNLMER.resid~MAG.$Ht,main="MAG",ylab=expression(paste("Residuals(",m^2,")",sep="")),xlim=c(0,40),xlab="TOPHT (m)")
points(MAGNLMEF.resid~MAG.$Ht,pch=16)
MAGNLMEF.resid
#
plot(SMAGNLMER.resid~SMAG.$Ht,main="SMAG",ylab=expression(paste("Residuals(",m^2,")",sep="")),xlim=c(0,40),xlab="TOPHT (m)")
points(SMAGNLMEF.resid~SMAG.$Ht,pch=16)
plot(DCLNLMER.resid~DCL.$Ht,main="DCL",ylab=expression(paste("Residuals(",m^2,")",sep="")),xlim=c(0,40),xlab="TOPHT (m)")
points(DCLNLMEF.resid~DCL.$Ht,pch=16)

plot(VALNLMER.resid~VAL.$Ht,main="VAL",ylab=expression(paste("Residuals(",m^2,")",sep="")),xlim=c(0,40),xlab="TOPHT (m)")
points(VALNLMEF.resid~VAL.$Ht,pch=16)
# SAP.MK
plot(SAP.MKNLMER.resid~SAP.MK.$Ht,main="SAP.MK",ylab=expression(paste("Residuals(",m^2,")",sep="")),xlim=c(0,40),xlab="TOPHT (m)")
points(SAP.MKNLMEF.resid~SAP.$Ht,pch=16)

plot(SMAG.MKNLMER.resid~SMAG.MK.$Ht,main="SMAG.MK",ylab=expression(paste("Residuals(",m^2,")",sep="")),xlim=c(0,40),xlab="TOPHT (m)")
points(SMAG.MKNLMEF.resid~SMAG.MK.$Ht,pch=16)
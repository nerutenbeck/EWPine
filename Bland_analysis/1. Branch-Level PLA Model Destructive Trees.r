library(RODBC)
temp=odbcConnectAccess("C:/Users/Frank/Desktop/Bland WhitePine Leaf Area Database.mdb")
branch.data=sqlFetch(temp,"Branch Data")
sample.branches=sqlFetch(temp,"Sample Branches for BLA model_to R")
Tree.at=sqlFetch(temp,"Tree-Level Attributes")
BLARMSE=sqlFetch(temp,"BLARMSE")
close(temp)    
'Year'=as.factor('Year')
'Location'=as.factor('Location')
'Tree'=as.factor('Tree')
'Season'=as.factor('Season')
library(nlme)


###########################################################################
Tree.at=subset(Tree.at,Study!="Guiterman")
head(Tree.at)
mean(Tree.at$TreeAge)
min(Tree.at$TreeAge)
max(Tree.at$TreeAge)
SE<-sd(Tree.at$TreeAge)/sqrt(length(Tree.at$TreeAge))
SE


################################  Preliminary Analysis#################################
head(sample.branches)
plot(BLAm2~RDINC,data=sample.branches)
plot(BLAm2~BD,data=sample.branches)
plot(BLAm2~DINC,data=sample.branches)
plot(BLAm2~Season,data=sample.branches)
plot(BLAm2~Location,data=sample.branches)
plot(BLAm2~Year,data=sample.branches)
plot(BLAm2~Tree,data=sample.branches)
plot(BLAm2~Ht,data=sample.branches)
plot(BLAm2~DBH,data=sample.branches)
###########################   Regressions for branch leaf area (BLA)  ##############################################
#AIC= -264.646 No a3 parameter
brlafixednoa3=gnls((BLAm2)~(a0*BD^a1)*(RDINC^(a2-1)),data=sample.branches,weights=varPower(0.5,form="BD"),  
start=c(a0=29.336702,a1=1.053348,a2=0.976824),control=nlmeControl(minScale=1e-10,returnObject=T))
summary(brlafixednoa3)
# Fixed effects only AIC = -319.3024
brlafixed=gnls((BLAm2)~(a0*BD^a1)*(RDINC^(a2-1))*exp(-((a3)*RDINC^a2)),data=sample.branches,weights=varPower(0.5,form="BD"),
start=c(a0=0.8295703 ,a1=1.0262252,a2=1.5297188 ,a3=1.3671490),control=nlmeControl(minScale=1e-10,returnObject=T))
summary(brlafixed)
####################################################################################
#AIC= -347.1729   Guiterman's branch model  
brlaranG=nlme((BLAm2)~(a0*BD^a1)*(RDINC^(a2-1))*exp(-((a3)*RDINC^a2)),data=sample.branches,
start=c(a0=.8,a1=1.8798685,a2=1.8722691,a3=1.7450574),control=nlmeControl(minScale=1e-10,returnObject=T),
fixed=a0+a1+a2+a3~1, weights=varPower(0.5,form="BD"),random=a1+a3~1|thinType/Tree)
summary(brlaranG)
# One outlying value, possible FES-2
plot(brlaranG)
pairs(brlaranG)
ranef(brlaranG)
fixef(brlaranG)
source("C:/Users/Frank/Desktop/Tinn R Files/furnival.r")
furnival(brlaranG)
######################################################################################
# Random Location AIC= -325.6742
brlaranL=nlme((BLAm2)~(a0*BD^a1)*(RDINC^(a2-1))*exp(-((a3)*RDINC^a2)),data=sample.branches,
start=c(a0=.8,a1=1.8798685,a2=1.8722691,a3=1.7450574),control=nlmeControl(minScale=1e-10,returnObject=T),
fixed=a0+a1+a2+a3~1, weights=varPower(0.5,form="BD"),random=a3~1|Location)
summary(brlaranL)
# Best model w/ random effect on Location and Tree AIC= -340.2372
brlaranLT=nlme((BLAm2)~(a0*BD^a1)*(RDINC^(a2-1))*exp(-((a3)*RDINC^a2)),data=sample.branches,
start=c(a0=1,a1=1.8798685,a2=1.8722691,a3=1.7450574),control=nlmeControl(minScale=1e-10,returnObject=T),
fixed=a0+a1+a2+a3~1, weights=varPower(0.5,form="BD"),random=a1+a3~1|Location/Tree)
summary(brlaranLT)
plot(brlaranLT)
qqnorm(brlaranLT)
ranef(brlaranLT)
fixef(brlaranLT)
# Random =a3~1|Location/Year/Tree) AIC = -275.6531
brlaranLYT=nlme((BLAm2)~(a0*BD^a1)*(RDINC^(a2-1))*exp(-((a3)*RDINC^a2)),data=sample.branches,
start=c(a0=80,a1=1.8798685,a2=1.8722691,a3=1.7450574),control=nlmeControl(minScale=1e-10,returnObject=T),
fixed=a0+a1+a2+a3~1, weights=varPower(0.5,form="BD"),random=a1+a3~1|Location/Year/Tree)
summary(brlaranLYT)
########################################Branch Foliar Mass#####################################
##Branch foliar mass equation   
brfmG=nlme((Fmass)~(a0*BD^a1)*(RDINC^(a2-1))*exp(((a3)*RDINC^a2)),data=sample.branches, ,fixed=a0+a1+a2+a3~1,random=a1+a3~1|thinType/Tree, weights=varPower(0.5,form="BD"),
start=c(a0=80.66,a1=1.91,a2=1.87,a3=1.903),control=nlmeControl(minScale=1e-10,returnObject=T))
summary(brfmG)
plot(brfmG)
qqnorm(brfmG)
ranef(brfmG)
fixef(brlaran)
# Predict pBLA and pBFM for sample.branches
sample.branches$pBLAnoa3=(predict(brlafixednoa3,newdata=sample.branches,na.action=na.omit,levels=0))^2
sample.branches$pBLAfixed=(predict(brlafixed,newdata=sample.branches,na.action=na.omit,levels=2))^2
# Predict
sample.branches$pBLAranG=(predict(brlaranG,newdata=sample.branches,na.action=na.omit,levels=2))^2
sample.branches$pBLALT=(predict(brlaranLT,newdata=sample.branches,na.action=na.omit,levels=2))^2
sample.branches$pBLALYT=(predict(brlaranLYT,newYdata=sample.branches,na.action=na.omit,levels=3))^2
sample.branches$pBFM=(predict(brfm,newdata=sample.branches,na.action=na.omit,levels=2))^2
head(sample.branches)
  
sample.branches$pTFM=(predict(brfm,newdata=sample.branches,na.action=na.omit,levels=2))^2
############# Residual analysis by Season, Height, DBH and BD  ###################
# Non-linear relationship
plot(sample.branches$BLA,sample.branches$pBLA,pch=16,col="black",xlab='BLA',ylab='Predicted BLA')
abline(0,1)
sample.branches$resid = sample.branches$BLA - sample.branches$pBLA   #obs-pred
# Plot residuals by Season   
sample.branches$residuals<-resid(brlafixednoa3)                                    
bla.winter<-sample.branches[sample.branches$Season=="Feb",]
bla.summer<-sample.branches[sample.branches$Season=="Aug",]
plot(bla.winter$BLA,bla.winter$resid,pch=16,col="black",ylab="Residuals", xlab="Fitted",ylim=c(-5000,5000))
points(bla.summer$BLA,bla.summer$resid,pch=16,col="red",) 
legend("topleft",col=c("black","red"), pch=c(16,16), legend=c("Feb", "Aug"))
# Plot residuals by Height
Ht.res<- resid(brlaranLT)
scatter.smooth(sample.branches$Ht_Tot, Ht.res, ylab="Residuals", xlab="Total Height", ylim=c(-50,50)) 
# Plot residuals by DBH
DBH.res<- resid(brlaranLT)
scatter.smooth(sample.branches$DBH, DBH.res, ylab="Residuals", xlab="DBH", ylim=c(-50,50))
#Residuals by BD
###
# Residuals by RDINC- DINC
# sequence of 
bd=seq(0.1,200,0.5)
rdinc=seq(0.01,0.99,0.01)
aaa=merge(bd,dinc)
predict
# Residuals over Guiterman's 
#abline(h=0)
BD.res<-resid(brlaranLT)
plot(sample.branches$BD, BD.res, ylab="Residuals", xlab="BD", ylim=c(-50,50)) 

############################## Tree-level Equations ####### ####################
#BLA
treeran<-read.csv("C:/Users/Frank/Desktop/CSV Files/tree noa3 fixed only.csv")   
tree=merge(branch.data,treeran,by='Tree',all=T)
tree$pBLA=((tree$a0*tree$BD^tree$a1)*(tree$RDINC^(tree$a2-1)))
head(tree)
#########################summation of branch LA by tree using predict and aggregate
tree.LA=aggregate(x=tree$pBLA,by=list(tree$Tree,tree$DBH,tree$Ht_Tot,tree$CL),FUN=sum)
head(tree.LA)
#the 'aggregate' function deletes names (don't know why), must rename columns accordingly
names(tree.LA)[1:5] <- c("Tree","DBH","Ht","CL","predLA")
head(tree.LA)

#BFM
brfm.ran<-read.csv("C:/Users/Frank/Desktop/CSV Files/brfm.ran.csv")
treeTFM=merge(branch.pred,brfm.ran,by='Tree',all=T) 
head(treeTFM)
treeTFM$pBFM=((treeTFM$a0*treeTFM$BD^treeTFM$a1)*(treeTFM$RDINC^(treeTFM$a2-1))*exp(-(treeTFM$a3)*treeTFM$RDINC^treeTFM$a2))
head(treeTFM)
tree.TFM=aggregate(x=treeTFM$pBFM,by=list(treeTFM$Tree,treeTFM$DBH,treeTFM$Ht_Tot,treeTFM$CL),FUN=sum)
names(tree.TFM)[1:5] <- c("Tree","DBH","Ht","CL","TFM")
head(tree.TFM)                                                                                           
tree.totals=merge(tree.LA,tree.TFM,by=c('Tree','DBH','Ht'),all=T)
head(tree.totals)
#converts g values to kg values
tree.totals$TFM=tree.totals$TFM*0.001
head(tree.totals)
write.csv(tree.LA, file = ("C:/Users/Frank/Desktop/TLAL-T.csv"), row.names = FALSE) 

###############################  Model Validation  #####################################
# Predicts BLA using various model forms
sample.branches$pBLAfixed=(predict(brlafixed,newdata=sample.branches,na.action=na.omit,levels=0))^2
sample.branches$pBLAL=(predict(brlaranL,newdata=sample.branches,na.action=na.omit,levels=1))^2
sample.branches$pBLALT=(predict(brlaranLT,newdata=sample.branches,na.action=na.omit,levels=2))^2
sample.branches$pBLALYT=(predict(brlaranLYT,newdata=sample.branches,na.action=na.omit,levels=3))^2


rmse<-function(obs,pred) {
    sqrt(sum((obs-pred)^2,na.rm=T)/length(obs))
    }
mab<-function(obs,pred) {
    (sum(abs(obs-pred),na.rm=T))/length(obs)
    }
bias<-function(obs,pred) {
    (sum((obs-pred),na.rm=T))/length(obs)
    }
mpb<-function(obs,pred) {
    100*((sum(abs(obs-pred),na.rm=T))/(sum(obs,na.rm=T)))
    }

rmse.G<-rmse(obs=BLARMSE$BLAmeas,pred=BLARMSE$BLApred)
rmse.G
rmse.fixed<-rmse(obs=sample.branches$BLA,pred=sample.branches$pBLAfixed)
rmse.L<-rmse(obs=sample.branches$BLA,pred=sample.branches$pBLAL)
rmse.LT<-rmse(obs=sample.branches$BLA,pred=sample.branches$pBLALT)
rmse.LYT<-rmse(obs=sample.branches$BLA,pred=sample.branches$pBLALYT)
rmse.fixed
rmse.L
rmse.LT
rmse.LYT
MAB.fixed<-mab(obs=sample.branches$BLA,pred=sample.branches$pBLAfixed)
MAB.L<-mab(obs=sample.branches$BLA,pred=sample.branches$pBLAL)
MAB.LT<-mab(obs=sample.branches$BLA,pred=sample.branches$pBLALT)
MAB.LYT<-mab(obs=sample.branches$BLA,pred=sample.branches$pBLALYT)
MAB.fixed
MAB.L
MAB.LT
MAB.LYT

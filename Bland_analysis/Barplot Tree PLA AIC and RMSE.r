AIC=read.csv("C:/Users/Frank/Desktop/CSV Files/Tree Model AIC.csv")
RMSE=read.csv("C:/Users/Frank/Desktop/CSV Files/Tree Model RMSE.csv")
AAD=read.csv("C:/Users/Frank/Desktop/CSV Files/AAD Litterfall.csv")
AIC1=read.csv("C:/Users/Frank/Desktop/CSV Files/Stand Model AIC.csv")
RMSE1=read.csv("C:/Users/Frank/Desktop/CSV Files/Stand Model RMSE.csv")
StandAIC=read.csv("C:/Users/Frank/Desktop//CSV Files/Stand LAI AIC.csv")
############################################################################################
par(mfrow=c(2,1))

# This plots AIC from stand level LAI models from chapter 2
head(AIC1)
str(AIC1)
AIC1.BP=AIC1[,-1]
barplot(as.matrix(t(AIC1.BP)),names.arg=c("8","7","6","5","4","3","2","1"),las=1,cex.names=0.9,col=c("black"),horiz=T,xlab="AIC",beside=TRUE,xlim=c(0,250))
# mean AIC value =191.4
abline(v=191.4,lty=2)

# This plots RMSE and MAB from stand level LAI models in chapter 2
head(RMSE1)
str(RMSE1)
RMSE1.BP=RMSE1[,-1]
barplot(as.matrix(t(RMSE1.BP)),names.arg=c("8","7","6","5","4","3","2","1"),las=1,cex.names=0.8,col=c("grey","black"),horiz=T,xlab=expression(paste("RMSE(", m^2,m^-2 ,")"),sep=""),beside=TRUE,xlim=c(0,1))
legend("bottomright",c("RMSE","MAB"),fill=c("black","gray"),cex=0.9, bty="n")
# mean AIC value =191.4
abline(v=0.705,lty=2)
abline(v=0.570375,lty=2,col="gray")

head(AIC)
str(AIC)
AIC.BP=AIC[,-1]
barplot(as.matrix(t(AIC.BP)),names.arg=c("SMAG.MK","SAP.MK","VAL","MAG","DCL","SMAG","SCL","SAP"),las=1,cex.names=0.75,col=c("black","gray","white"),horiz=T,xlab="AIC",beside=TRUE,xlim=c(0,750))
legend("bottomright",c("NLS","WNLS","NLME-R"),fill=c("white","gray","black"),cex=0.9, bty="n")
# Adds average AIC lines
abline(v=534.4,lty=2)#NLS
abline(v=428.7,lty=1,col="gray")#WNLS
abline(v=422.9,lty=1)#NLME

head(RMSE)                                                                                             
str(RMSE)
RMSE.BP=RMSE[,-1]

barplot(as.matrix(t(RMSE.BP)),names.arg=c("SMAG.MK","SAP.MK","DCL","VAL","MAG","SMAG","SCL","SAP"),las=1,cex.names=0.75,beside=TRUE,xlim=c(0,75), 
col=c("black","gray25","gray75","white"),horiz=T,xlab=expression(paste("RMSE(",m^2,")",sep="")))
legend("bottomright", c("NLS","WNLS","NLME-R","NLME-F"),cex=1,bty="n", fill=c("white","gray75","gray25","black"))

plot(AIC$NLS,AIC$WNLS)

head(AAD)
str(AAD)
AAD.BP=AAD[,-1]
barplot(as.matrix(t(AAD.BP)),names.arg=c("SMAG.MK","SAP.MK","DCL","VAL","SMAG","MAG","SCL","SAP"),las=1,cex.names=0.7,beside=TRUE,xlim=c(0,3),
col=c("black","gray","white"),horiz=T,xlab=expression(paste("Average Absolute Deviation(", m^2,m^-2 ,")"),sep=""))
legend("bottomright", c("NLS","WNLS","NLME-F"),cex=1,bty="n", fill=c("white","gray","black"))

# This works, but not with legend
barplot(as.matrix(t(AAD.BP)),names.arg=c("SAP","SCL","MAG","SMAG","VAL","DCL"),beside=TRUE,xlim=c(0,9),
col=c("black","gray","white"),horiz=T,main="2011(16 plots)",xlab=expression(paste("Average Absolute Deviation(", m^2,m^-2 ,")"),sep=""))

# This plots AIC for stand-level models including those with TOPHT
head(StandAIC)
str(StandAIC)
SAIC.BP=StandAIC[,-1]
barplot(as.matrix(t(SAIC.BP)),names.arg=c("SAP.MK","SMAG","SCL","SAP"),las=1,cex.names=0.75,col=c("black","gray","white"),horiz=T,xlab="AIC",beside=TRUE,xlim=c(0,95))
legend("bottomright",c("Base","TOPHT","TOPHT-R"),fill=c("white","gray","black"),cex=0.9, bty="n")
# Adds average AIC lines
abline(v=534.4,lty=2)#NLS
abline(v=428.7,lty=1,col="gray")#WNLS
abline(v=422.9,lty=1)#NLME

# This plots RMSE for stand-level models including those with TOPHT
StandRMSE=read.csv("C:/Users/Frank/Desktop/CSV Files/Stand LAI RMSE.csv")
head(StandRMSE)
str(StandRMSE)
SRMSE.BP=StandRMSE[,-1]
barplot(as.matrix(t(SRMSE.BP)),names.arg=c("SAP.MK","SMAG","SCL","SAP"),las=1,cex.names=0.75,col=c("black","gray","white"),horiz=T,xlab=expression(paste("RMSE(", m^2,m^-2 ,")"),sep=""),beside=TRUE,xlim=c(0,0.42))
legend("bottomright",c("Base","TOPHT","TOPHT-R"),fill=c("white","gray","black"),cex=0.9, bty="n")

# This plots MAB for stand-level models including those with TOPHT
StandMAB=read.csv("C:/Users/Frank/Desktop/CSV Files/Stand LAI MAB.csv")
head(StandMAB)
str(StandMAB)
SMAB.BP=StandMAB[,-1]
barplot(as.matrix(t(SMAB.BP)),names.arg=c("SAP.MK","SMAG","SCL","SAP"),las=1,cex.names=0.75,col=c("black","gray","white"),horiz=T,xlab=expression(paste("Mean Absolute Bias(", m^2,m^-2 ,")"),sep=""),beside=TRUE,xlim=c(0,0.2))
legend("bottomright",c("Base","TOPHT","TOPHT-R"),fill=c("white","gray","black"),cex=0.9, bty="n")

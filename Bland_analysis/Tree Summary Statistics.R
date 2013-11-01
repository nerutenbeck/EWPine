library(RODBC)
temp=odbcConnectAccess("C:/Users/Frank/Desktop/Demeritt WP- Bland.mdb")
trees=sqlFetch(temp,"Trees")
plots=sqlFetch(temp,"Plots")
close(temp)

trees=subset(trees,Status=='1')
# Merges tree data with plot information
TSS=merge(trees,plots,by=c('Plot'))
# Renames duplicate outputs
names(TSS)[names(TSS)=="Stand.x"] = "Stand"
names(TSS)[names(TSS)=="PlotID.x"] = "PlotID"
#This extracts important information from resulting dataframe
TSS=subset(TSS,,select=c(ID,Year,Stand,Stand_ID,Plot,PlotID,thinType,Subplot,TreeNo,TreeNoOld,DBH,LLB,LLW,TotHt,Status,AreaExpansion,AreaExpansionEng,Planted,SIeng,SI))
# Changes Area expansion to TPH and TPA based on plot size
names(TSS)[names(TSS)=="AreaExpansion"] = "TPH"
names(TSS)[names(TSS)=="AreaExpansionEng"] = "TPA"
TSS$CL_LLB=TSS$TotHt-TSS$LLB
TSS$CL_LLW=TSS$TotHt-TSS$LLW
TSS$CR_LLB=TSS$CL_LLB/TSS$TotHt
TSS$CR_LLW=TSS$CL_LLW/TSS$TotHt
TSS$mLCR=TSS$CL_LLB/(TSS$TotHt-1.37)

# This grabs Tree volume calculations
source("C:/Users/Frank/Desktop/Tinn R Files/Volume Calculations.r")

Volume=subset(Volume,,select=c(ID,Year,Plot,TreeNo,KozakVOB,LWKHonerVOB,HonerVIB,KozakVIB,LWKHonerVIB,KozakVOBeng,LWKHonerVOBeng,HonerVIBeng,KozakVIBeng,LWKHonerVIBeng))
TSS=merge(TSS,Volume,by=c('ID','Year','Plot','TreeNo'))

TSS$TreeBAm=(0.0000785398*(TSS$DBH^2))
TSS$TreeBAft=(0.005454*(TSS$DBH/2.54)^2)

TSS$TreeBAperHa=(0.0000785398*(TSS$DBH^2))*TSS$TPH
TSS$TreeBAperAc=(0.005454*(TSS$DBH/2.54)^2)*TSS$TPA
###############################################################################################

write.csv(TSS, file = ("C:/Users/Frank/Desktop/R output/TSS.csv"), row.names = FALSE)
################################################################################################
###########################    Density Metrics    ##############################################
library(plyr);library(doBy)
# See allometric LAI.r for details on these plyr and doBy functions

# Sums TPA, TPH TreeBAperHa and TreeBAperAc by Year and Plot from TSS dataframe
Density=ddply(TSS,.(Year,Stand,PlotID,Plot,thinType),colwise(sum,c('TPA','TPH','TreeBAperHa','TreeBAperAc')),.progress="text")
# Calculates density from TSS table
Volume=ddply(TSS,.(Year,Stand,PlotID,Plot,thinType),colwise(mean,c('KozakVOB','LWKHonerVOB','HonerVIB','KozakVIB','LWKHonerVIB','KozakVOBeng','LWKHonerVOBeng','HonerVIBeng','KozakVIBeng','LWKHonerVIBeng')),.progress="text")
# Merges Density and Volume dataframes
DM=merge(Density,Volume,by=c('Year','Stand','PlotID','Plot','thinType'))

# Rename tree volumes and basal area to plot means or totals 
names(DM)[names(DM)=="KozakVOB"] ="MeanKozakVOB";names(DM)[names(DM)=="LWKHonerVOB"]="MeanLWKHonerVOB";
names(DM)[names(DM)=="HonerVIB"]="MeanHonerVIB";names(DM)[names(DM)=="KozakVIB"]= "MeanKozakVIB";
names(DM)[names(DM)=="LWKHonerVIB"]= "MeanLWKHonerVIB";names(DM)[names(DM)=="KozakVOBeng"]="MeanKozakVOBeng"
names(DM)[names(DM)=="LWKHonerVOBeng"]= "MeanLWKHonerVOBeng";names(DM)[names(DM)=="HonerVIBeng"]= "MeanHonerVIBeng"
names(DM)[names(DM)=="KozakVIBeng"]= "MeanKozakVIBeng";names(DM)[names(DM)=="LWKHonerVIBeng"]= "MeanLWKHonerVIBeng"
names(DM)[names(DM)=="TreeBAperHa"]= "BAperHa";names(DM)[names(DM)=="TreeBAperAc"]= "BAperAc"
# Calculates QMD
DM$QMD=((DM$BAperHa/DM$TPH)/0.0000785398)^0.5
DM$QMDeng=DM$QMD*0.3937008
# Add SI, Planted from plot table
SI=ddply(TSS,.(Year,Stand,PlotID,Plot,thinType),colwise(mean,c('SI','SIeng','Planted')),.progress="text")
DM=merge(SI,DM,by=c('Year','Stand','PlotID','Plot','thinType'))
DM$LogSI=log(DM$SI)/log(10)
DM$LogSIeng=log(DM$SIeng)/log(10)

# Calculate TOPHT 
# Calculate SI from TOPHT trees
top1a=subset(TSS,TPH==250,select=,);top1a$Year=factor(top1a$Year)
top1b=subset(TSS,TPH==131.578947,select=,);top1b$Year=factor(top1b$Year)
top1=subset(TSS,TPH==100,select=,);top1$Year=factor(top1$Year)
top4=subset(TSS,TPH==25,select=,);top4$Year=factor(top4$Year)
top3=subset(TSS,TPH==33.333,select=,);top3$Year=factor(top3$Year)

TOPHT1a=ddply(top1a,.(Year,PlotID),subset,rank(desc(TotHt),ties.method= "random")<=1,.progress="text")
TOPHT1b=ddply(top1b,.(Year,PlotID),subset,rank(desc(TotHt),ties.method= "random")<=1,.progress="text")
TOPHT1=ddply(top1,.(Year,PlotID),subset,rank(desc(TotHt),ties.method= "random")<=1,.progress="text")
TOPHT3=ddply(top3,.(Year,PlotID),subset,rank(desc(TotHt),ties.method= "random")<=3,.progress="text")
TOPHT4=ddply(top4,.(Year,PlotID),subset,rank(desc(TotHt),ties.method= "random")<=4,.progress="text")  

MTOPHT1a=ddply(TOPHT1a,.(Year,PlotID),summarise,TOPHT=mean(TotHt,na.rm=F),.progress="text")
MTOPHT1b=ddply(TOPHT1b,.(Year,PlotID),summarise,TOPHT=mean(TotHt,na.rm=F),.progress="text")
MTOPHT1=ddply(TOPHT1,.(Year,PlotID),summarise,TOPHT=mean(TotHt,na.rm=F),.progress="text")
MTOPHT3=ddply(TOPHT3,.(Year,PlotID),summarise,TOPHT=mean(TotHt,na.rm=F),.progress="text")
MTOPHT4=ddply(TOPHT4,.(Year,PlotID),summarise,TOPHT=mean(TotHt,na.rm=F),.progress="text")

TOPHT=rbind(MTOPHT1a,MTOPHT1b)
TOPHT=rbind(TOPHT,MTOPHT1)
TOPHT=rbind(TOPHT,MTOPHT3)
TOPHT=rbind(TOPHT,MTOPHT4)
# Merges TOPHT with rest of DM
DM=merge(DM,TOPHT,by=c('Year','PlotID'),all=T)
###########################################################
DM$TOPHTeng=(DM$TOPHT)*3.2808399
DM$LogTPA=log(DM$TPA)/log(10)
DM$LnTPA=log(DM$TPA)
DM$LnTPH=log(DM$TPH)
DM$LogTPH=log(DM$TPH)/log(10)
DM$LogMeanKozakVOB=log(DM$MeanKozakVOB)/log(10)
DM$LogMeanLWKHonerVOB=log(DM$MeanLWKHonerVOB)/log(10)
DM$LogMeanHonerVIB=log(DM$MeanHonerVIB)/log(10)
DM$LogMeanKozakVIB=log(DM$MeanKozakVIB)/log(10)
DM$LogMeanLWKHonerVIB=log(DM$MeanLWKHonerVIB)/log(10)
DM$LnMeanLWKHonerVIB=log(DM$MeanLWKHonerVIB)
DM$LogMeanKozakVOBeng=log(DM$MeanKozakVOBeng)/log(10)
DM$LogMeanLWKHonerVOBeng=log(DM$MeanLWKHonerVOBeng)/log(10)
DM$LogMeanHonerVIBeng=log(DM$MeanHonerVIBeng)/log(10)
DM$LnMeanLWKHonerVIBeng=log(DM$MeanLWKHonerVIBeng)
DM$LogMeanKozakVIBeng=log(DM$MeanKozakVIBeng)/log(10)
DM$LogMeanLWKHonerVIBeng=log(DM$MeanLWKHonerVIBeng)/log(10)
########### Relative Density Calculations  ##############################
# Calculate parameters from Density Management Diagram-Equations.R
# MaxLogTPA = (MeanVol-Intercept)/Slope
# Max TPA = 10^MaxLogTPA
# RD = TPA/MaxTPA

# Kozak SFA VOB 
DM$SFAKozakVOBMaxLogTPA=(DM$LogMeanKozakVOBeng-4.92180034)/-1.31657
DM$SFAKozakVOBMaxTPA=10^DM$SFAKozakVOBMaxLogTPA
DM$SFAKozakVOBRD=DM$TPA/DM$SFAKozakVOBMaxTPA
# Kozak SFA VIB
DM$SFAKozakVIBMaxLogTPA=(DM$LogMeanKozakVIBeng-4.841510523)/-1.31277
DM$SFAKozakVIBMaxTPA=10^DM$SFAKozakVIBMaxLogTPA
DM$SFAKozakVIBRD=DM$TPA/DM$SFAKozakVIBMaxTPA
# LWK Honer SFA VIB
DM$SFALWKHonerVIBMaxLogTPA=(DM$LogMeanLWKHonerVIBeng-4.771724345)/-1.30562
DM$SFALWKHonerVIBMaxTPA=10^DM$SFALWKHonerVIBMaxLogTPA
DM$SFALWKHonerVIBRD=DM$TPA/DM$SFALWKHonerVIBMaxTPA
# Honer SFA VOB 
DM$SFALWKHonerVOBMaxLogTPA=(DM$LogMeanLWKHonerVOBeng-4.703160681)/-1.27356
DM$SFALWKHonerVOBMaxTPA=10^DM$SFALWKHonerVOBMaxLogTPA
DM$SFALWKHonerVOBRD=DM$TPA/DM$SFALWKHonerVOBMaxTPA
# Honer SFA VIB
DM$SFAHonerVIBMaxLogTPA=(DM$LogMeanHonerVIBeng-4.566475712)/-1.2381
DM$SFAHonerVIBMaxTPA=10^DM$SFAHonerVIBMaxLogTPA
DM$SFAHonerVIBRD=DM$TPA/DM$SFAHonerVIBMaxTPA
####################################################################################
# Kozak VOB 
DM$KozakVOBMaxLogTPA=(DM$LogMeanKozakVOBeng-4.97718)/-1.26855
DM$KozakVOBMaxTPA=10^DM$KozakVOBMaxLogTPA
DM$KozakVOBRD=DM$TPA/DM$KozakVOBMaxTPA
# Kozak VIB
DM$KozakVIBMaxLogTPA=(DM$LogMeanKozakVIBeng-4.90034)/-1.26532
DM$KozakVIBMaxTPA=10^DM$KozakVIBMaxLogTPA
DM$KozakVIBRD=DM$TPA/DM$KozakVIBMaxTPA
# LWK Honer VIB
DM$LWKHonerVIBMaxLogTPA=(DM$LogMeanLWKHonerVIBeng-4.83978)/-1.31731
DM$LWKHonerVIBMaxTPA=10^DM$LWKHonerVIBMaxLogTPA
DM$LWKHonerVIBRD=DM$TPA/DM$LWKHonerVIBMaxTPA
# LWK Honer VOB 
DM$LWKHonerVOBMaxLogTPA=(DM$LogMeanLWKHonerVOBeng-4.76304)/-1.28224
DM$LWKHonerVOBMaxTPA=10^DM$LWKHonerVOBMaxLogTPA
DM$LWKHonerVOBRD=DM$TPA/DM$LWKHonerVOBMaxTPA
# Honer VIB
DM$HonerVIBMaxLogTPA=(DM$LogMeanHonerVIBeng-4.80447)/-1.29832
DM$HonerVIBMaxTPA=10^DM$HonerVIBMaxLogTPA
DM$HonerVIBRD=DM$TPA/DM$HonerVIBMaxTPA
######################################################################
DM$InnesMaxLogTPAeng=(DM$LnMeanLWKHonerVIBeng-10.2)/-1.51
DM$InnesMaxTPAeng=2.71828^DM$InnesMaxLogTPAeng
DM$InnesRD=DM$TPA/DM$InnesMaxTPAeng

# Curtis RD (1982)- 0.00007854*TPH*QMD^1.5)
DM$CurtisRD=0.00007854*DM$TPH*(DM$QMD^1.5)
# SDI- TPH*(QMD/25)^1.6 
DM$SDI=DM$TPH*(DM$QMD/25)^1.6
# Wilson's (1946) spacing coefficient
DM$WilsonSpacing=(10000*DM$TPH)^0.5/DM$TOPHT

##########################################################################################
write.csv(DM, file = ("C:/Users/Frank/Desktop/R output/DM.csv"), row.names = FALSE)

# Summary statistics for Table 1 in chapter 3
# Merge litterfall LAI data with DM
library(RODBC)
temp=odbcConnectAccess("C:/Users/Frank/Desktop/Demeritt WP- Bland.mdb")
LAI=sqlFetch(temp,"Plot_litterfall LAI_all years")
DM=merge(DM,LAI,by=c("Year","Stand","PlotID","Plot"))


###########
S1990=subset(DM,Stand=="1990")
min(S1990$TOPHTeng)
max(S1990$TOPHTeng)
mean(S1990$TOPHTeng)
mean(S1990$Litterfall_LAI)

S22=subset(DM,Plot=="2x2")
min(S22$TOPHTeng)
max(S22$TOPHTeng)
mean(S22$TOPHTeng)
mean(S22$Litterfall_LAI)

S22n=subset(DM,Plot=="2x2 New")
min(S22n$TOPHTeng)
max(S22n$TOPHTeng)
mean(S22n$TOPHTeng)
mean(S22n$Litterfall_LAI)

S33=subset(DM,Plot=="3x3")
min(S33$TOPHTeng)
max(S33$TOPHTeng)
mean(S33$TOPHTeng)
mean(S33$Litterfall_LAI)

S70=subset(DM,Plot=="1970")
min(S70$TOPHTeng)
max(S70$TOPHTeng)
mean(S70$TOPHTeng)
mean(S70$Litterfall_LAI)
     
S70L=subset(DM,Plot=="LD")
min(S70L$TOPHTeng)
max(S70L$TOPHTeng)
mean(S70L$TOPHTeng)
mean(S70L$Litterfall_LAI)

S70L=subset(DM,thinType=="LD")
S70L=subset(S70L,Stand=="Thinning")
min(S70L$TOPHTeng)
max(S70L$TOPHTeng)
mean(S70L$TOPHTeng)
mean(S70L$Litterfall_LAI)

Nutting=subset(DM,Plot=="Nutting")
min(Nutting$TOPHTeng)
max(Nutting$TOPHTeng)
mean(Nutting$TOPHTeng)
mean(Nutting$Litterfall_LAI)

# Thinning study plots
T=subset(DM,Stand=="Thinning")
head(T)

write.csv(C1, file = ("C:/Users/Frank/Desktop/R output/Control.TOPHT.csv"), row.names = FALSE)

TLD=subset(T,thinType=="LD")
mean(TLD$TOPHTeng)
mean(TLD$Litterfall_LAI)

TB=subset(T,thinType=="B")
mean(TB$TOPHTeng)
mean(TB$Litterfall_LAI)

MM=subset(DM,Plot=="MM")
mean(MM$TOPHTeng)
mean(MM$Litterfall_LAI)

# Old-Field
OF1=subset(TOPHT,PlotID=="1")
OF2=subset(TOPHT,PlotID=="2")
OF=merge(OF1,OF2,by="Year")
mean(OF$TOPHT)

# HD
HD1=subset(DM,Plot=="HD1")
HD2=subset(DM,Plot=="HD2")
HD3=subset(DM,Plot=="HD3")
HDone=merge(HD1,HD2,by="Year")
HD=merge(HD3,HDone,by="Year")


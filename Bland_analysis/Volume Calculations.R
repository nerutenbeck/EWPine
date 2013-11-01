library(RODBC);library(plyr);temp=odbcConnectAccess("C:/Users/Frank/Desktop/Demeritt WP- Bland.mdb")
trees=sqlFetch(temp,"Trees");plots=sqlFetch(temp,"Plots");close(temp)

# Loads Kozak VOB and eventually VIB calculations from ## KozakVOB and KOZAKVIB.R files ##
KozakVOB=read.csv("C:/Users/Frank/Desktop/R Output/KozakVOB.csv")
KozakVOB=subset(KozakVOB,,select=c(ID,Year,Stand,PlotID,Plot,Subplot,TreeNo,TreeNoOld,DBH,LLB,LLW,TotHt,Status,KozakVOB))
KozakVIB=read.csv("C:/Users/Frank/Desktop/R Output/KozakVIB.csv")
KozakVIB=subset(KozakVIB,,select=c(ID,Year,Stand,PlotID,Plot,Subplot,TreeNo,TreeNoOld,DBH,LLB,LLW,TotHt,Status,KozakVIB))

# Takes out extraneous information
trees=subset(trees,,select=c(ID,Year,Stand,PlotID,Plot,Subplot,TreeNo,TreeNoOld,DBH,LLB,LLW,TotHt,Status))
# Removes dead trees
trees=subset(trees,Status=="1")

#  Metric versions- 2 VOB and 3 VIB equations total Only 2 originally in metric
# HonerVIB (1967)
trees$HonerVIB=((trees$DBH/2.54)^2)/(0.691+(363.676/(trees$TotHt*3.2808399)))/35.31466672
# Honer VOB in m^3 from Li and Weiskittel
trees$LWKHonerVOB=((trees$DBH/2.54)^2)/(0.971+(346.08/(trees$TotHt*3.2808399)))/35.31466672
# Honer VIB in m^3  from Li and Weiskitell 2012
trees$LWKHonerVIB=((trees$DBH/2.54)^2)/(0.397+(410.15/(trees$TotHt*3.2808399)))/35.31466672
# Honer VIB in m^3 from Honer (1967)- I shouldn't use this because BH is 1.37 meters
#trees$HonerVIB=(0.004319*(trees$DBH^2)/(0.691+(110.848/(trees$TotHt))))
                    
# HonerVIB from Honer (1963)
trees$HonerVIBeng=((trees$DBH/2.54)^2)/(0.691+(363.676/(trees$TotHt*3.2808399)))
# Honer VOB in ft^3 from Honer (1963)
trees$LWKHonerVOBeng=((trees$DBH/2.54)^2)/(0.971+(346.08/(trees$TotHt*3.2808399)))
# Honer VIB in ft^3 from Honer (1967)
trees$LWKHonerVIBeng=((trees$DBH/2.54)^2)/(0.397+(410.15/(trees$TotHt*3.2808399)))
OF=subset(trees,Stand=="Old-Field")
head(OF)
OF
########################################################################################
# Adds Kozak VOB and VIB volumes from Li and Weiskittell 2012
# This loses 1 tree
Volume=merge(trees,KozakVOB,by=c('ID','Year','Stand','PlotID','Plot','Subplot','TreeNo','TreeNoOld','DBH','Status'))
Volume=merge(Volume,KozakVIB,by=c('ID','Year','Stand','PlotID','Plot','Subplot','TreeNo','TreeNoOld','DBH','Status'))

names(Volume)[names(Volume)=="TotHt.x"] = "TotHt"
Volume$KozakVOBeng=Volume$KozakVOB*35.31466672;Volume$KozakVIBeng=Volume$KozakVIB*35.31466672
Volume=subset(Volume,,select=c(ID,Year,Stand,PlotID,Plot,Subplot,TreeNo,TreeNoOld,DBH,TotHt,KozakVOB,LWKHonerVOB,KozakVOBeng,LWKHonerVOBeng,KozakVIB,LWKHonerVIB,HonerVIB,KozakVIBeng,LWKHonerVIBeng,HonerVIBeng))         
#write.csv(Volume, file = ("C:/Users/Frank/Desktop/R output/Volume.csv"), row.names = FALSE)

####################################################################################################
# Merge with plot data
head(plots)
head(Volume)
# Something is wrong with the Old-Field plots
Volume=merge(Volume,plots,by=c('Stand','Plot','PlotID'),all=T)
head(Volume)
# Sum tree volume to plot-level
PlotVol=ddply(Volume,.(Year,Stand,PlotID,Plot),colwise(sum,c('KozakVIB','KozakVIBeng','LWKHonerVIB','LWKHonerVIBeng')),.progress="text")






# Merge plot volume with plot attribute info
Vol=merge(PlotVol,plots,by=c('Plot','Stand','PlotID'),all=T) 
# subset relevant info
Vol=subset(Vol,,select=c('Year','Stand','Plot','PlotID','KozakVIB','KozakVIBeng','LWKHonerVIB','LWKHonerVIBeng','AreaExpansionEng'))
head(Vol)

Vol$KozakVIB=Vol$KozakVIB*Vol$AreaExpansionEng
Vol$KozakVIBeng=Vol$KozakVIBeng*Vol$AreaExpansionEng
Vol$LWKHonerVIB=Vol$LWKHonerVIB*Vol$AreaExpansionEng
Vol$LWKHonerVIBeng=Vol$LWKHonerVIBeng*Vol$AreaExpansionEng

Vol












###########################################################
###
### White Pine Growth Efficiecy Over Time
### Nathan E. Rutenbeck
### University of Maine School of Forest Resources
### 2013
###
###########################################################

rm(list=ls())
options(show.signif.stars=F)
windows()

library(ggplot2)
library(reshape2)
library(doBy)
library(plyr)

################## Data preparation

# Import plots
plots.load<-read.table('Plots.txt',head=T,sep='\t') # Load plot data
plots.cl<-plots.load[,c(1:8,15)] # Drop unused columns
names(plots.cl)<-c('plotid','stand','plot','standid','trt','plot.ha','plot.m2','exp.fac','SI') # rename variables
plots<-plots.cl # Clean plots object

summary(plots)

# Import litterfall leaf areas from plots, clean variables, change names, reorder
litter.load<-read.table('Plot_litterfall.txt',head=T,sep='\t') # litterfall data
litter.cl<-litter.load[,c(1:7)] # Clean litter object
names(litter.cl)<-c('year','stand','trt','plot','plotid','LAI','LAI.SE') # rename variables
litter<-orderBy(~plotid+year,data=litter.cl) # clean up order

summary(litter)

lai.plot<-ggplot(litter,aes(x=year,y=LAI))+geom_point()+
  geom_errorbar(aes(ymax=LAI+LAI.SE,ymin=LAI-LAI.SE))
lai.plot+facet_wrap(~plot) # LAI trends over time at the plot level

# Import tree list and measurements from database
DBtrees.load<-read.table('Trees.txt',head=T,sep='\t') # tree measurements
DBtrees.cl<-DBtrees.load[,c(3:8,2,9:12)] # tree.old matches tree number in windendro data....
names(DBtrees.cl)<-c('stand','plotid','plot','subplot','tree','tree.old','year','DBH','LLB','LLW','ht') # rename variables
DBtrees<-orderBy(~subplot+tree+year,data=DBtrees.cl) #clean up order

summary(DBtrees)

DBtrees<-subset(DBtrees,DBH>0) # Select trees with diameters >0

# Read in ring data, get unique index, add calculated values

win<-read.table('windendro.txt',head=T,sep='\t') # ring data
win.cl<-win[,c(1:6,13,38:ncol(win))] # Remove unused columns
winmelt<-melt(win.cl,id=1:8,na.rm=T) # Reshape data to desired structure
winmelt$Age<-unclass(winmelt$variable) # translate variables into ages
winmelt.cl<-winmelt[,c(4:6,1,3,7:8,10:11)] # Drop unused columns
names(winmelt.cl)<-c('stand','plot','subplot','tree', # Rename columns
                     'pathid','ringcount','pathlength','seglength','age')
summary(winmelt.cl)
winmelt.cl[winmelt.cl$seglength>20,] #Look for outliers
winmelt.cl<-winmelt.cl[winmelt.cl$seglength<20,] # Remove ridiculous outliers (only two)
hist(winmelt.cl$seglength) # Look at histogram of segment lengths

# Make sure ring data identifiers match those of the Database trees
names(winmelt.cl)[4]<-'tree.old' # rename tree number to match field in database

levels(winmelt.cl$stand)%in%levels(DBtrees$stand) # Confirm that stand names match
levels(winmelt.cl$stand)
levels(DBtrees$stand)
levels(DBtrees$stand)[7]<-'Nutting' # Change stand name in DBtrees object
levels(winmelt.cl$stand)%in%levels(DBtrees$stand) # Re-confirm that stand names match

levels(winmelt.cl$plot)%in%levels(DBtrees$plot) # Confirm that plot names match
levels(winmelt.cl$subplot)%in%levels(DBtrees$subplot) # Confirm that subplot names match

treeids<-data.frame(treeid=1:nrow(unique(winmelt.cl[1:4])), # Give each tree a unique ID number
                    orderBy(~stand+plot+subplot+tree.old,
                            data=unique(winmelt.cl[1:4])))
head(treeids)
nrow(treeids)

rings.temp1<-merge(winmelt.cl,treeids,by=c('stand','plot','subplot','tree.old')) # Merge ring data with treeids
rings.temp2<-ddply(rings.temp1,.(treeid),mutate,year=2011-max(age)+age) # Calculate year for each ring

rings.temp3<-orderBy(~stand+plot+subplot+tree.old+year,data=rings.temp2) # Order correctly by tree ID and year

# Calculate mean segment length and area increment for each year.
rings.mean<-ddply(rings.temp3,.(stand,plot,subplot,tree.old,treeid,year,age),summarize,
                  mean.seg=mean(seglength),
                  ainc.mm=pi*I(mean(seglength)^2)) 

# Predicted DIB in a given year will be the cumulative sum of average diameter increments
rings<-ddply(rings.mean,.(treeid),mutate,
             RIB=cumsum(mean.seg)/10,
             DIB=2*cumsum(mean.seg)/10,
             ainc.cm=ainc.mm/100) # Calculate predicted radius and diameter
summary(rings)
treeids$treeid%in%rings$treeid # nothing dropped...


# Merge to get measured diameters and ring data all in the same place

trees.both<-merge(DBtrees,rings,by=c('stand','plot','subplot','tree.old','year'))

summary(trees.both)


levels(rings$stand)%in%levels(trees.both$stand) # Confirm stands match
levels(rings$plot)%in%levels(trees.both$plot) # confirm plots match
levels(rings$subplot)%in%levels(trees.both$subplot) # Confirm subplots match
miss<-treeids$treeid%in%trees.both$treeid # Lost some trees!
lost<-treeids[!miss,]
lost # Only two trees in rings but not in database
# I have to live with it. I can't find plausible matches

###########################################################
###
### White Pine Growth Efficiecy Over Time
### Nathan E. Rutenbeck
### University of Maine School of Forest Resources
### 2013
###
###########################################################

source('lme4_models.R')
windows()

######################### Growth Efficiency

#Calculate annual volume increments

all.trees<-merge(DBtrees[,c(1,3:4,6:8,11)],rings[-c(8:10)],
                 by=c('stand','plot','subplot','tree.old','year'),all.y=T) # Merge database trees with rings information
names(all.trees)
head(all.trees)
mod.coefs<-merge(ht.mod.coefs,DBH.mod.coefs,by=('treeid'),all.y=T) # Merge coefficients from DBH and height models
summary(mod.coefs)
treeids$treeid%in%mod.coefs$treeid # Just missing two trees...

tree.master<-merge(mod.coefs,all.trees,by='treeid') # Merge model coefficients with trees
tree.master<-orderBy(~treeid+age,data=tree.master) # reorder
head(tree.master,100)
summary(tree.master)
treeids$treeid%in%tree.master$treeid # Still just missing two trees

names(tree.master)

# diameter

tree.master$ht.pred<-with(tree.master,ht.b0+I(ht.b1*DIB)+I(ht.b2*age)) # predict heights from model coefficients
tree.master$dia.pred<-with(tree.master,d.b0+d.b1*DIB+d.b2*ht.pred+d.b3*age) # predict diameters from model coefficients

tree.master$vol.pred<-with(tree.master,
                           ((dia.pred/2.54)^2)/(0.971+(346.08/(ht.pred*3.2808399)))/35.31466672) # predict volumes from model coefficients


# Calculate plot-level volumes

plot.vols<-ddply(tree.master,.(stand,plot,year),summarize,
                 ba=sum(ainc.cm),
                 vol=sum(vol.pred)) # Predict volumes and basal areas at the plot level

head(plot.vols,100)
levels(plot.vols$plot) # These are PLOT VOLUMES in each given year, not increment!!!!

plot.vols$vinc<-c(plot.vols$vol[1],rep(NA,nrow(plot.vols)-1)) # Make empty vector to store increment info

for(i in 2:nrow(plot.vols)){ 
  ifelse(plot.vols$plot[i]==plot.vols$plot[i-1], # Check if ID is same as previous
         plot.vols$vinc[i]<-plot.vols$vol[i]-plot.vols$vol[i-1], # If same as previous, VINC= present-previous
         plot.vols$vinc[i]<-plot.vols$vol[i]) # If it's not the same as previous, it equals itself
} # 
head(plot.vols)
plot.vols$sge<-with(plot.vols,vinc/vol) # Calculate Size Growth Efficiency (increment/present volume)


summary(plot.vols)
volsmelt<-na.omit(melt(plot.vols[,-c(4:5)],id=1:3)) # Melt plot volumes for plotting purposes

long.sge<-ggplot(volsmelt,aes(y=value,x=year,lty=variable))+geom_line(cex=0.7)+
  facet_wrap(~plot)+ylim(0,2)+xlim(1950,2012)+
  scale_linetype_discrete(labels=c(expression('VINC ('*m^3*')'),'Size GE'))
long.sge # not sure how much this means looking way back, but it's interesting to see.

# Merge with litter data

levels(litter$plot)%in%levels(all.trees$plot) # Chapman plot is missing

ge<-merge(litter,plot.vols,by=c('plot','year')) # Merge litter data with plot volumes
levels(ge$plot)

head(ge,100)
ge$F.GE<-with(ge,vinc/LAI) # Calculate foliar growth efficiency
ge$S.GE<-with(ge,vinc/vol) # Calculate size growth efficiency
names(ge)[3]<-'stand' # Rename a column
names(ge)

ge.rshp<-melt(ge[,c(1:5,11,12)],id=1:5) # Reshape the data for plotting
head(ge.rshp)

fge.plot<-ggplot(ge.rshp[ge.rshp$plot!='Blue-Spring',],aes(y=value,x=year,color=trt,lty=variable))+
  geom_line(cex=0.7)+facet_wrap(~plot)+
  scale_linetype_discrete(labels=c(expression('VINC ('*m^3*')'),
                                   'Foliar GE'))+
  scale_color_discrete(labels=c('low density','light thin','no thin','PCT'),
                       name='treatment')

ge.rshp2<-melt(ge[,c(1:5,11,13)],id=1:5)
head(ge.rshp2)

sge.plot<-ggplot(ge.rshp2[ge.rshp2$plot!='Blue-Spring',],aes(y=value,x=year,color=trt,lty=variable))+
  geom_line(cex=0.7)+facet_wrap(~plot)+
  scale_linetype_discrete(labels=c(expression('VINC ('*m^3*')'),
                                   'Size GE'))+
  scale_color_discrete(labels=c('low density','light thin','no thin','PCT'),
                       name='treatment')


pdf('GE_plots.pdf',height=8,width=10)

lai.plot+facet_wrap(~plot) # LAI trends over time at the plot level
fge.plot # Foliar growth efficiency
sge.plot # Size growth efficiency (modern)
# long.sge # Size growth efficiency (historic and probably uninformative)

# Make individual plots
plts<-levels(ge$plot)[c(1:9,26:34)]

par(ask=F)
for(i in 1:length(plts)){
  tmp.i=subset(ge.rshp,plot==plts[i])
  tmp2.i=subset(ge.rshp2,plot==plts[i])
  print(ggplot(tmp.i,aes(y=value,x=year,lty=variable))+
          geom_line(cex=0.8)+ggtitle(plts[i])+labs(linetype='')+ylab('')+
          scale_linetype_discrete(labels=c(expression('VINC ('*m^3*')'),'Foliar GE')))
  print(ggplot(tmp2.i,aes(y=value,x=year,lty=variable))+
          geom_line(cex=0.8)+ggtitle(plts[i])+labs(linetype='')+ylab('')+
          scale_linetype_discrete(labels=c(expression('VINC ('*m^3*')'),'Size GE')))
   
}


dev.off()
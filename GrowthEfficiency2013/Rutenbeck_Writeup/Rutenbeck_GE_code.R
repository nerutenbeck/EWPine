###########################################################
###
### Influence of White Pine Growth Efficiecy Over Time
### Nathan E. Rutenbeck
### University of Maine School of Forest Resources
### Spring 2013
###
###########################################################

rm(list=ls())
update.packages()
options(show.signif.stars=F)
windows()

library(ggplot2)
library(reshape)
library(doBy)
library(arm)
library(R2jags)
library(MCMCpack)



################## Data preparation

# Import plots
plots.load<-read.table('Plots.txt',head=T,sep='\t')
plots.cl<-plots.load[,c(1:8,15)]
names(plots.cl)<-c('plotid','stand','plot','standid','trt','plot.ha','plot.m2','exp.fac','SI')
plots<-plots.cl

summary(plots)

# Import litterfall leaf areas from plots, clean variables, change names, reorder
litter.load<-read.table('Plot_litterfall.txt',head=T,sep='\t') # litterfall data
litter.cl<-litter.load[,c(1:7)]
names(litter.cl)<-c('year','stand','trt','plot','plotid','LAI','LAI.SE')
litter<-orderBy(~plotid+year,data=litter.cl)

summary(litter)

lai.plot<-ggplot(litter,aes(x=year,y=LAI))+geom_point()+
  geom_errorbar(aes(ymax=LAI+LAI.SE,ymin=LAI-LAI.SE))
lai.plot+facet_wrap(~plot) # LAI trends over time at the plot level

# Import tree list and measurements from database
DBtrees.load<-read.table('Trees.txt',head=T,sep='\t') # tree measurements
DBtrees.cl<-DBtrees.load[,c(3:8,2,9:12)] # tree.old matches tree number in windendro data....
names(DBtrees.cl)<-c('stand','plotid','plot','subplot','tree','tree.old','year','DBH','LLB','LLW','ht')
DBtrees<-orderBy(~subplot+tree+year,data=DBtrees.cl) #clean up order

summary(DBtrees)

DBtrees<-subset(DBtrees,DBH>0)

# Read in ring data, get unique index, add calculated values

win<-read.table('windendro.txt',head=T,sep='\t') # ring data
win.cl<-win[,c(1:6,13,38:ncol(win))] # Remove unused columns
winmelt<-melt(win.cl,id=1:8,na.rm=T) # Reshape data to desired structure
winmelt$Age<-unclass(winmelt$variable) # translate variables into ages
winmelt.cl<-winmelt[,c(4:6,1,3,7:8,10:11)]
names(winmelt.cl)<-c('stand','plot','subplot','tree',
                     'pathid','ringcount','pathlength','seglength','age')
summary(winmelt.cl)
winmelt.cl[winmelt.cl$seglength>20,]
winmelt.cl<-winmelt.cl[winmelt.cl$seglength<20,] # Remove ridiculous outliers (only two)
hist(winmelt.cl$seglength)

# Make sure ring data identifiers match those of the Database trees
names(winmelt.cl)[4]<-'tree.old' # rename tree number to match field in database

levels(winmelt.cl$stand)%in%levels(DBtrees$stand)
levels(winmelt.cl$stand)
levels(DBtrees$stand)
levels(DBtrees$stand)[7]<-'Nutting'
levels(winmelt.cl$stand)%in%levels(DBtrees$stand)

levels(winmelt.cl$plot)%in%levels(DBtrees$plot)
levels(winmelt.cl$subplot)%in%levels(DBtrees$subplot)

treeids<-data.frame(treeid=1:nrow(unique(winmelt.cl[1:4])),
                                  orderBy(~stand+plot+subplot+tree.old,
                                          data=unique(winmelt.cl[1:4])))
head(treeids)
nrow(treeids)

rings.temp1<-merge(winmelt.cl,treeids,by=c('stand','plot','subplot','tree.old'))
rings.temp2<-ddply(rings.temp1,.(treeid),mutate,year=2011-max(age)+age) #Calculate year ring formed

rings.temp3<-orderBy(~stand+plot+subplot+tree.old+year,data=rings.temp2) # Order correctly by individual and year

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


levels(rings$stand)%in%levels(trees.both$stand)
levels(rings$plot)%in%levels(trees.both$plot)
levels(rings$subplot)%in%levels(trees.both$subplot)
miss<-treeids$treeid%in%trees.both$treeid # Lost some!
lost<-treeids[!miss,]
lost # Only two trees in rings but not in database
# I have to live with it. I can't find plausible matches


########################### Model Fitting

# my own residual plot function
rplot<-function(mod,...){
  plot(resid(mod)~fitted(mod),pch=20,...)
  abline(0,0,lty=3)
}

# Mean Absolute Bias function
MAB<-function(model){
  mab<-sum(abs(resid(model)))/length(resid(model))
  return(mab)
}

# Root Mean Square Error
RMSE<-function(model){
  rmse<-sqrt(I(sum(I(resid(model)^2)))/length(resid(model)))
  return(rmse)
}

# Diameter model

DBH.mod1<-lmer(DBH~DIB+(1+DIB|treeid),data=trees.both)
DBH.mod2<-lmer(DBH~DIB+ht+(1+DIB+ht|treeid),data=trees.both)
DBH.mod3<-lmer(DBH~DIB+ht+age+(1+DIB+ht+age|treeid),data=trees.both)
anova(DBH.mod1,DBH.mod2,DBH.mod3)

pdf('./DBH_fit.pdf',height=8,width=10)
par(mfrow=c(3,1))
rplot(DBH.mod1,ylim=c(-3,3), main='DBH~DIB')
rplot(DBH.mod2,ylim=c(-3,3), main='DBH~DIB+ht')
rplot(DBH.mod3,ylim=c(-3,3), main='DBH~DIB+ht+age') # Definitely the best.
par(mfrow=c(1,1))

dbh.resids<-melt(data.frame('DIB'=resid(DBH.mod1),'DIB+ht'=resid(DBH.mod2),'DIB+ht+age'=resid(DBH.mod3)))

dbh.resid.hist<-ggplot(dbh.resids,aes(x=value))+geom_histogram(binwidth=0.1)+
  xlab('DBH bias (cm)')+ggtitle('Residuals of Diameter Models')+facet_wrap(~variable)

dbh.fitted<-melt(data.frame('DIB'=fitted(DBH.mod1),
                           'DIB.ht'=fitted(DBH.mod2),'DIB.ht.age'=fitted(DBH.mod3)))
dbh.fitted$obs<-rep(trees.both$DBH,3)

dbh.fit.plot<-ggplot(dbh.fitted,aes(x=obs,y=value))+geom_point()+
  geom_abline(slope=1,intercept=0)+facet_wrap(~variable)+xlab('Observed DBH')+ylab('Fitted DBH')

dbh.fit.num<-melt(data.frame(model=c('DIB','DIB+ht','DIB+ht+age'),
                            'AIC'=AIC(DBH.mod1,DBH.mod2,DBH.mod3)[[2]],
                            'BIC'=BIC(DBH.mod1,DBH.mod2,DBH.mod3)[[2]]),id='model')

dbh.fit.plot2<-ggplot(dbh.fit.num,aes(x=variable,y=value))+
  geom_bar(stat='identity',aes(fill=variable))+facet_wrap(~model)

dbh.fit.num2<-melt(data.frame(model=c('DIB','DIB+ht','DIB+ht+age'),
                              'RMSE'=c(RMSE(DBH.mod1),RMSE(DBH.mod2),RMSE(DBH.mod3)),
                              'MAB'=c(MAB(DBH.mod1),MAB(DBH.mod2),MAB(DBH.mod3))))

dbh.fit.plot3<-ggplot(dbh.fit.num2,aes(x=variable,y=value))+
  geom_bar(stat='identity',aes(fill=variable))+facet_wrap(~model)

dbh.resid.hist
dbh.fit.plot
dbh.fit.plot2
dbh.fit.plot3 
dev.off()

print(dbh.fit.num)
print(dbh.fit.num2) # compare numerically

display(DBH.mod1)
display(DBH.mod1)
display(DBH.mod3) # Will use this model for volume prediction


DBH.mod.coefs<-data.frame(treeid=rownames(coef(DBH.mod3)$treeid),
                          coef(DBH.mod3)$treeid[1],coef(DBH.mod3)$treeid[2],
                          coef(DBH.mod3)$treeid[3],coef(DBH.mod3)$treeid[4])
names(DBH.mod.coefs)<-c('treeid','d.b0','d.b1','d.b2','d.b3')

head(DBH.mod.coefs) # these are the varying coefficients for each of the J trees

# I would like to display the model with the data for each tree, but that is hard to do
# with so many predictors

# DBH.plot<-ggplot(trees.both,aes(x=DIB,y=DBH))+
#   geom_point()+geom_abline(data=DBH.mod.coefs,aes(intercept=DBH.ahat,slope=DBH.bhat))+
#   facet_wrap(~treeid)

# DBH.plot # shows tree level slope, intercept differences, individual prediction bias

# Height models

ht.mod1<-lmer(ht~age+(1+age|treeid),trees.both)
ht.mod2<-lmer(ht~DIB+(1+DIB|treeid),trees.both)
ht.mod3<-lmer(ht~DIB+age+(1+DIB+age|treeid),trees.both)
anova(ht.mod1,ht.mod2,ht.mod3) # Here it looks like model 3 performs the best.

pdf('./Ht_fit.pdf', height=8,width=10)

par(mfrow=c(3,1))
rplot(ht.mod1,ylim=c(-3,3),main='ht~age')
rplot(ht.mod2,ylim=c(-3,3),main='ht~DIB')
rplot(ht.mod3,ylim=c(-3,3),main='ht~DIB+age')
par(mfrow=c(1,1))

# Residuals above could be better...
# I also wonder if I could improve height predictions by incorporating stand and plot effects.
ht.mod4<-update(ht.mod3,.~.+stand)
ht.mod5<-update(ht.mod3,.~.+plot)
ht.mod6<-update(ht.mod3,.~.+subplot)
# Can't add multiple layers of plot/subplot because of fitting problems from correlation matrix.

anova(ht.mod3,ht.mod4,ht.mod5,ht.mod6) # Here it looks like model 5 would improve things
# Theoretically I should add plot factors as random effects to reflect experimental structure? 
# I don't have time now to explore

#rplot(ht.mod5,ylim=c(-3,3),main='ht~plot+DIB+age')

ht.resids<-melt(data.frame('age'=resid(ht.mod1),'DIB'=resid(ht.mod2),
                           'DIB+age'=resid(ht.mod3),'DIB+age+plot'=resid(ht.mod5)))

ht.resid.hist<-ggplot(ht.resids,aes(x=value))+geom_histogram(binwidth=0.1)+
  xlab('ht bias (meters)')+ggtitle('Residuals of Height Models')+facet_wrap(~variable,nrow=1)

ht.fitted<-melt(data.frame('age'=fitted(ht.mod1),
                           'DIB'=fitted(ht.mod2),
                           'DIB.age'=fitted(ht.mod3),
                           'DIB.age.plot'=fitted(ht.mod5)))
ht.fitted$obs<-rep(trees.both$ht,4)

ht.fit.plot<-ggplot(ht.fitted,aes(x=obs,y=value))+geom_point()+
  geom_abline(slope=1,intercept=0)+facet_wrap(~variable)+
  xlab('Observed Height')+ylab('Fitted Height')

ht.fit.num<-melt(data.frame(model=c('age','DIB','DIB.age','DIB.age.plot'),
                            'AIC'=AIC(ht.mod1,ht.mod2,ht.mod3,ht.mod5)[[2]],
                            'BIC'=BIC(ht.mod1,ht.mod2,ht.mod3,ht.mod5)[[2]]),id='model')

ht.fit.plot2<-ggplot(ht.fit.num,aes(x=variable,y=value))+
  geom_bar(stat='identity',aes(fill=variable))+facet_wrap(~model,nrow=1)

ht.fit.num2<-melt(data.frame(model=c('age','DIB','DIB.age','DIB.age.plot'),
                             'RMSE'=c(RMSE(ht.mod1),RMSE(ht.mod2),RMSE(ht.mod3),RMSE(ht.mod5)),
                             'MAB'=c(MAB(ht.mod1),MAB(ht.mod2),MAB(ht.mod3),MAB(ht.mod5))))

ht.fit.plot3<-ggplot(ht.fit.num2,aes(x=variable,y=value))+
  geom_bar(stat='identity',aes(fill=variable))+facet_wrap(~model,nrow=1)

ht.resid.hist
ht.fit.plot
ht.fit.plot2
ht.fit.plot3

dev.off()

print(ht.fit.num)
print(ht.fit.num2) # compare numerically

display(ht.mod1)
display(ht.mod2)
display(ht.mod3)
display(ht.mod5) # Definitely select this one eventually, 
# but I don't have time to build the model matrix right now...

ht.mod.coefs<-data.frame(treeid=rownames(coef(ht.mod3)$treeid),coef(ht.mod3)$treeid)
                         coef(ht.mod5)$treeid[1],coef(ht.mod5)$treeid[2]
names(ht.mod.coefs)<-c('treeid','ht.b0','ht.b1','ht.b2') #coefficients for heights

head(ht.mod.coefs)

########## Bayesian fit for the simple diameter model. Eventually I'd like to try it with the full model
# I went through a couple iterations of this model fitting, and finally
# had to scale the variables and model the correlation with the inverse Wishart distribution to get it
# to converge properly.

DBH.jm<-function(){
  for (i in 1:n){
    y[i]~dnorm(yhat[i],tau.y)
    yhat[i]<-a[treeid[i]]+b[treeid[i]]*x[i]    
  }
  tau.y<-pow(sigma.y,-2)
  sigma.y~dunif(0,100)
  
  for(j in 1:J){
    a[j]<-xi.a*B.raw[j,1] 
    b[j]<-xi.b*B.raw[j,2]
    B.raw[j,1:2]~dmnorm(B.raw.hat[j,],Tau.B.raw[,])
    B.raw.hat[j,1]<-mu.a.raw
    B.raw.hat[j,2]<-mu.b.raw
  }
  mu.a<-xi.a*mu.a.raw # scaled mean slope
  mu.b<-xi.b*mu.b.raw # scaled mean slope
  mu.a.raw~dnorm(0,0.0001) # raw mean for intercepts
  mu.b.raw~dnorm(0,0.0001) # raw mean for slopes
  
  xi.a~dunif(0,100) # Scaling parameters
  xi.b~dunif(0,100)
  
  Tau.B.raw[1:2,1:2]~dwish(W[,],df) # modeling Tau directly with the Wishart distribution
  df<-3 # degrees of freedom for the distribution of Tau.B.raw
  
  Sigma.B.raw[1:2,1:2]<-inverse(Tau.B.raw[,]) # Calculating raw sigmas
  sigma.a<-xi.a*sqrt(Sigma.B.raw[1,1]) # caluculating scaled standard deviations for the intercepts
  sigma.b<-xi.b*sqrt(Sigma.B.raw[2,2]) # calculating scaled standard deviations for the slopes
  rho<-Sigma.B.raw[1,2]/sqrt(Sigma.B.raw[1,1]*Sigma.B.raw[2,2]) # Correlation
}

dbh.prep<-orderBy(~treeid,data=na.omit(trees.both[,c(12,8,17)]))
head(dbh.prep)

dbh.dat<-with(dbh.prep,list(n=nrow(dbh.prep),
                            J=length(unique(treeid)),
                            y=DBH,
                            x=DIB, # This makes the model not converge... Need to scale
                            #x=scale(DIB,scale=T)[,1], # No convergence by just subt. mean...
                            treeid=as.numeric(unclass(as.factor(treeid))),
                            W=diag(2)))
str(dbh.dat)
dbh.inits<-function(){
  list(B.raw=array(rnorm(2*540),c(540,2)), # This is a Jx2 matrix of rnorm samples
       mu.a.raw=rnorm(1),
       mu.b.raw=rnorm(1),
       sigma.y=runif(1),
       Tau.B.raw=rwish(3,diag(2)), # rwish() from the MCMCpack package
       xi.a=runif(1),
       xi.b=runif(1))
}

dbh.params<-c('a','b','mu.a','mu.b','sigma.a','sigma.b','sigma.y','rho')

dbh.fit<-jags(data=dbh.dat,inits=dbh.inits,dbh.params,model.file=DBH.jm,n.iter=500)
print(dbh.fit)

max(dbh.fit$BUGSoutput$summary[,8]) # Max Rhat - should be below 1.1 

dbh.upd<-update(dbh.fit,n.iter=1000) # If convergence is problematic

str(dbh.fit)

dbh.jagscoefs<-list(mean=dbh.fit$BUGSoutput$mean,sd=dbh.fit$BUGSoutput$sd,median=dbh.fit$BUGSoutput$median)
str(dbh.jagscoefs)
dbh.jagscoefs$mean

### Try the more complex model. As of the end of the term, this still doesn't work.

DBH.jm2<-function(){
  for (i in 1:n){
    y[i]~dnorm(yhat[i],tau.y)
    yhat[i]<-a[treeid[i]]+b[treeid[i]]*DIB[i]+c[treeid[i]]*ht[i]+d[treeid[i]]*age[i]    
  }
  tau.y<-pow(sigma.y,-2)
  sigma.y~dunif(0,100)
  
  for(j in 1:J){
    a[j]<-xi.a*B.raw[j,1] 
    b[j]<-xi.b*B.raw[j,2]
    c[j]<-xi.c*B.raw[j,3]
    d[j]<-xi.d*B.raw[j,4]
    B.raw[j,1:4]~dmnorm(B.raw.hat[j,],Tau.B.raw[,])
    B.raw.hat[j,1]<-mu.a.raw
    B.raw.hat[j,2]<-mu.b.raw
    B.raw.hat[j,3]<-mu.c.raw
    B.raw.hat[j,4]<-mu.d.raw
  }
  mu.a<-xi.a*mu.a.raw 
  mu.b<-xi.b*mu.b.raw 
  mu.c<-xi.c*mu.c.raw
  mu.d<-xi.d*mu.c.raw
  mu.a.raw~dnorm(0,0.0001) 
  mu.b.raw~dnorm(0,0.0001) 
  mu.c.raw~dnorm(0,0.0001)
  mu.d.raw~dnorm(0,0.0001)
  
  xi.a~dunif(0,100) 
  xi.b~dunif(0,100)
  xi.c~dunif(0,100)
  xi.d~dunif(0,100)
  
  Tau.B.raw[1:4,1:4]~dwish(W[,],df) 
  df<-5 
  
  Sigma.B.raw[1:4,1:4]<-inverse(Tau.B.raw[,]) 
  sigma.a<-xi.a*sqrt(Sigma.B.raw[1,1]) 
  sigma.b<-xi.b*sqrt(Sigma.B.raw[2,2]) 
  sigma.c<-xi.c*sqrt(Sigma.B.raw[3,3])
  sigma.d<-xi.d*sqrt(Sigma.B.raw[4,4])

  rho.ab<-Sigma.B.raw[1,2]/sqrt(Sigma.B.raw[1,1]*Sigma.B.raw[2,2]) 
  rho.ac<-Sigma.B.raw[1,3]/sqrt(Sigma.B.raw[1,1]*Sigma.B.raw[3,3])
  rho.ad<-Sigma.B.raw[1,4]/sqrt(Sigma.B.raw[1,1]*Sigma.B.raw[4,4])
  rho.bc<-Sigma.B.raw[2,3]/sqrt(Sigma.B.raw[2,2]*Sigma.B.raw[3,3])
  rho.bd<-Sigma.B.raw[2,4]/sqrt(Sigma.B.raw[2,2]*Sigma.B.raw[4,4])
  rho.cd<-Sigma.B.raw[3,4]/sqrt(Sigma.B.raw[3,3]*Sigma.B.raw[4,4])
}

dbh.prep2<-orderBy(~treeid,data=na.omit(trees.both[,c(12,8,11,13,17)]))
head(dbh.prep2)

dbh.dat2<-with(dbh.prep2,list(n=nrow(dbh.prep2),
                            J=length(unique(treeid)),
                            y=DBH,
                            DIB=scale(DIB)[,1],
                            ht=scale(ht)[,1],
                            age=scale(age)[,1],
                            treeid=as.numeric(unclass(as.factor(treeid))),
                            W=diag(4)))
str(dbh.dat2)
dbh.inits2<-function(){
  list(B.raw=array(rnorm(4*540),c(540,4)), # This is a Jx2 matrix of rnorm samples
       mu.a.raw=rnorm(1),
       mu.b.raw=rnorm(1),
       mu.c.raw~rnorm(1),
       mu.d.raw~rnorm(1),
       sigma.y=runif(1),
       Tau.B.raw=rwish(5,diag(4)), # rwish() from the MCMCpack package
       xi.a=runif(1),
       xi.b=runif(1),
       xi.c=runif(1),
       xi.d=runif(1))
}

dbh.params2<-c('mu.a','mu.b','mu.c','mu.d','sigma.a','sigma.b','sigma.c','sigma.d','sigma.y',
              'rho.ab','rho.ac','rho.ad','rho.bc','rho.bd','rho.cd')

# This still doesn't fit for some reason. Something wrong with the initialization?

dbh.fit2<-jags(data=dbh.dat2,inits=dbh.inits2,dbh.params2,model.file=DBH.jm2,n.iter=50)

# After all that, I have to say that you gotta love lmer()....


######################### Growth Efficiency, finally...

#Calculate annual volume increments

all.trees<-merge(DBtrees[,c(1,3:4,6:8,11)],rings[-c(8:10)],
                 by=c('stand','plot','subplot','tree.old','year'),all.y=T)
names(all.trees)
head(all.trees)
mod.coefs<-merge(ht.mod.coefs,DBH.mod.coefs,by=('treeid'),all.y=T)
summary(mod.coefs)
treeids$treeid%in%mod.coefs$treeid # Just missing two trees...

tree.master<-merge(mod.coefs,all.trees,by='treeid')
tree.master<-orderBy(~treeid+age,data=tree.master)
head(tree.master,100)
summary(tree.master)
treeids$treeid%in%tree.master$treeid

names(tree.master)

# diameter

tree.master$ht.pred<-with(tree.master,ht.b0+I(ht.b1*DIB)+I(ht.b2*age))
tree.master$dia.pred<-with(tree.master,d.b0+d.b1*DIB+d.b2*ht.pred+d.b3*age)

tree.master$vol.pred<-with(tree.master,
                           ((dia.pred/2.54)^2)/(0.971+(346.08/(ht.pred*3.2808399)))/35.31466672)


# Calculate plot-level volumes

plot.vols<-ddply(tree.master,.(stand,plot,year),summarize,
                 ba=sum(ainc.cm),
                 vol=sum(vol.pred))

head(plot.vols,100)
levels(plot.vols$plot) # So these are volumes in each given year, not increment!!!!

plot.vols$vinc<-c(plot.vols$vol[1],rep(NA,nrow(plot.vols)-1))

for(i in 2:nrow(plot.vols)){
  ifelse(plot.vols$plot[i]==plot.vols$plot[i-1],
         plot.vols$vinc[i]<-plot.vols$vol[i]-plot.vols$vol[i-1],
         plot.vols$vinc[i]<-plot.vols$vol[i])
} # 
head(plot.vols)
plot.vols$sge<-with(plot.vols,vinc/vol)


summary(plot.vols)
volsmelt<-na.omit(melt(plot.vols[,-c(4:5)],id=1:3))

long.sge<-ggplot(volsmelt,aes(y=value,x=year,lty=variable))+geom_line(cex=0.7)+
  facet_wrap(~plot)+ylim(0,2)+xlim(1950,2012)+
  scale_linetype_discrete(labels=c(expression('VINC ('*m^3*')'),'Size GE'))
long.sge # not sure how much this means looking way back, but it's interesting to see.

# Merge with litter data

levels(litter$plot)%in%levels(all.trees$plot) # I don't know what Chapman is, so....

ge<-merge(litter,plot.vols,by=c('plot','year'))
levels(ge$plot)

head(ge,100)
ge$F.GE<-with(ge,vinc/LAI)
ge$S.GE<-with(ge,vinc/vol)
names(ge)[3]<-'stand'
names(ge)
ge.rshp<-melt.data.frame(ge[,c(1:5,11,12)],id=1:5)
head(ge.rshp)


ge.rshp<-melt.data.frame(ge[,c(1:5,11,12)],id=1:5)
head(ge.rshp)

fge.plot<-ggplot(ge.rshp[ge.rshp$plot!='Blue-Spring',],aes(y=value,x=year,color=trt,lty=variable))+
  geom_line(cex=0.7)+facet_wrap(~plot)+
  scale_linetype_discrete(labels=c(expression('VINC ('*m^3*')'),
                                   'Foliar GE'))+
  scale_color_discrete(labels=c('low density','light thin','no thin','PCT'),
                       name='treatment')

ge.rshp2<-melt.data.frame(ge[,c(1:5,11,13)],id=1:5)
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
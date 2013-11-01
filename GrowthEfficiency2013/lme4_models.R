###########################################################
###
### White Pine Growth Efficiecy Over Time
### Nathan E. Rutenbeck
### University of Maine School of Forest Resources
### 2013
###
###########################################################

source('Data_prep.R')
library(arm)

########################### Model Fitting with lme4

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

# Fit a series of diameter models to predict diameter in a given year

DBH.mod1<-lmer(DBH~DIB+(1+DIB|treeid),data=trees.both)
DBH.mod2<-lmer(DBH~DIB+ht+(1+DIB+ht|treeid),data=trees.both)
DBH.mod3<-lmer(DBH~DIB+ht+age+(1+DIB+ht+age|treeid),data=trees.both)
anova(DBH.mod1,DBH.mod2,DBH.mod3)

# Plot model fits for visual comparison

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
#                          coef(ht.mod5)$treeid[1],coef(ht.mod5)$treeid[2]

names(ht.mod.coefs)<-c('treeid','ht.b0','ht.b1','ht.b2') #coefficients for heights

head(ht.mod.coefs)
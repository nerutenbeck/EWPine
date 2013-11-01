###########################################################
###
### White Pine Growth Efficiecy Over Time
### Nathan E. Rutenbeck
### University of Maine School of Forest Resources
### 2013
###
###########################################################


library(R2jags)
library(MCMCpack)

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

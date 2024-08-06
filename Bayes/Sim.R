library('rstan')
library('boot')
library('Metrics')
nreps=100 #Number of simulations
tmax<-30 

#Containers
btrue<-matrix(NA,nreps,tmax)
strue<-matrix(NA,nreps,tmax)
ntrue<-matrix(NA,nreps,tmax+1)
BdataOut<-matrix(NA,nreps,tmax)

best<-bdiff<-array(NA,c(nreps,3000,tmax))
sest<-sdiff<-array(NA,c(nreps,3000,tmax))
nest<-ndiff<-array(NA,c(nreps,3000,tmax+1))

#loop over simulations
for (i in 1:nreps){
b<-exp(rnorm(tmax,-2,1)) ###Generate establishment rates
s<-inv.logit(rnorm(tmax,3,1)) ###Generate survival rates
lam<-b+s #lambda
Nstart<-sample(20:60,1) ###Generate starting population size
N<-c(Nstart,Nstart*cumprod(lam)) ###Predict population size
Bstar<-N[-length(N)]*b ###Establishment each interval
B<-Bstar*c(rev(cumprod(rev(s)))[-1],1) ###Take into account survival to observation period


Bdata<-rpois(length(B),B) ###Generate data

#save truth
btrue[i,]<-b
strue[i,]<-s
ntrue[i,]<-N
BdataOut[i,]<-Bdata

###Fit model
init<-list(list("s"=rep(0.95,length(B)),"b"=rep(0.2,length(B)),'Nstart'=50),list("s"=rep(0.95,length(B)),"b"=rep(0.2,length(B)),'Nstart'=50),list("s"=rep(0.95,length(B)),"b"=rep(0.2,length(B)),'Nstart'=50))

modeldata<-list("tmax"=length(BdataOut[i,]), "B"=BdataOut[i,])
options(mc.cores = parallel::detectCores())
fit1<-stan(file='SimFit.stan',data=modeldata, chains=3,iter=2000, init=init)
#launch_shinystan(fit1)

###Save pars from model
pars<-extract(fit1,pars=c('b','N','s'))

best[i,,]<-pars$b
sest[i,,]<-pars$s
nest[i,,]<-pars$N
bdiff[i,,]<-t(t(pars$b)-btrue[i,])
sdiff[i,,]<-t(t(pars$s)-strue[i,])
ndiff[i,,]<-t(t(pars$N)-ntrue[i,])

}



library('rstan')

wideB<-round(as.data.frame.matrix(t(xtabs(B~Time+dataset,s1,na.action = na.pass))))

Tmin<-apply(wideB,1,function(x){min(which(x>-1))})
Tmax<-apply(wideB,1,function(x){max(which(x>-1))})
Tmin1600<-Tmin
Tmin1600[which(Tmin1600<18)]<-18
stime<-cbind(Tmin,Tmax,Tmin1600)

Nout<-array(NA,c(3000,38,NPops))
bout<-bout.use<-array(NA,c(3000,37,NPops))
sout<-sout.use<-array(NA,c(3000,37,NPops))
for (i in 1:NPops){
tmax<-length(stime[i,1]:stime[i,2])  
init<-list(list("s"=rep(0.95,tmax),"b"=rep(0.2,tmax)),list("s"=rep(0.95,tmax),"b"=rep(0.2,tmax)),list("s"=rep(0.95,tmax),"b"=rep(0.2,tmax)))
modeldata<-list("tmax"=tmax, "B"=c(unlist(wideB[i,stime[i,1]:stime[i,2]])))

options(mc.cores = parallel::detectCores())
fit1<-stan(file='bayes/SimFit.stan',data=modeldata, chains=3,iter=2000, init=init,control = list(adapt_delta = 0.9, stepsize = 0.1))

pars<-rstan::extract(fit1,pars=c('b','N','s'))
Nout[,stime[i,1]:(stime[i,2]+1),i]<-pars$N
bout[,stime[i,1]:stime[i,2],i]<-pars$b
sout[,stime[i,1]:stime[i,2],i]<-pars$s

bout.use[,(stime[i,1]+1):stime[i,2],i]<-pars$b[,-1]
sout.use[,(stime[i,1]+1):stime[i,2],i]<-pars$s[,-1]
}


years<-seq(1260,1980,20)
post1600<-which(years>=1600)









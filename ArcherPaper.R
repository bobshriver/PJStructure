ArcherData<-read.csv('Archer1989.csv')


ArcherData$N<-cumsum(ArcherData$B)
years<-rev(seq(1985,1770,by=-20))
N<-ArcherData$N[match(years,ArcherData$Year)]

B<-round(c(diff(N)))



tmax<-length(B)  
init<-list(list("s"=rep(0.95,tmax),"b"=rep(0.2,tmax)),list("s"=rep(0.95,tmax),"b"=rep(0.2,tmax)),list("s"=rep(0.95,tmax),"b"=rep(0.2,tmax)))
modeldata<-list("tmax"=tmax, "B"=B)

options(mc.cores = parallel::detectCores())
fit1<-stan(file='bayes/SimFit.stan',data=modeldata, chains=3,iter=2000, init=init,control = list(adapt_delta = 0.9, stepsize = 0.1))

pars<-rstan::extract(fit1,pars=c('b','N','s'))

plot(years[-11],apply(pars$b,2,median),type='l',ylim=c(0,1.5),lwd=2,xlab='Year (CE)',ylab=expression("Establishment rate (20yrs"^-1~")"))
lines(years[-11],apply(pars$b,2,quantile,prob=0.975),type='l',lty=2,lwd=2)
lines(years[-11],apply(pars$b,2,quantile,prob=0.025),type='l',lty=2,lwd=2)

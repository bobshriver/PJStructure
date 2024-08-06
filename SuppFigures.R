###Code for producing supplemental figures#####3
library(boot)
library(RColorBrewer)
library(scales)
library(ggplot2)
library(gridExtra)
library(plyr)
library(plotrix)
library(segmented)


source("Functions.R")
agedataall<-read.csv("AgeSt2.csv") ####Load Data
agedataall<-agedataall[-which(agedataall$Census_keep=='n'),]###Remove any data after 1960 for datasets that are not complete census. Also remove trailing intervals that were not complete (e.g. last census was in 2003 and included trees from 2000-2010)

#####Aggregate data to 20 year intervals#####



bins<-seq(1260,2000,by=20)
agedataall20<-matrix(NA,1,15) ###Matrix that data will be appended to
colnames(agedataall20)<-colnames(agedataall)[-1]
for (z in 1:length(unique(agedataall$Dataset))){ #iterate through datasets
  
  if(agedataall$Interval[agedataall$Dataset==z][1]==20){ ##If dataset interval is already 20 all that is needed is to append that data on the table
    
    agedataall20=rbind.data.frame(agedataall20,agedataall[agedataall$Dataset==z,-1])
    
  }
  
  if(agedataall$Interval[agedataall$Dataset==z][1]==10){ #aggragate data for 10 year intervals
    duse<-agedataall[agedataall$Dataset==z,] #extract dataset
    tmin<-min(duse$Talign) 
    tmax<-max(duse$Talign)
    if(length(which(tmin==bins))==0){tmin=tmin-10} #if the the dataset starts on an odd interval (1710), move start year back by 10 (e.g. 1700). E.g. Since no trees were observed in 1700-1709, establishment is 0
    if(length(which(tmax==bins))>0){tmax=tmax-10} ## if the last year was not a complete 20 year interval (e.g. 1960-1969 not 1960-1979) then we have to cut off the last 10 year interval, making the end 1940-1959. Since 1960-1960 does not include a complete 20 year period it is not an accurate estimate of establshemnt. 
    zbins<-seq(tmin,tmax+10,20)
    
    for(t in zbins[-1]){ ###Iterate over bins
      B<-sum(duse$B[which(duse$Talign>=(t-20) & duse$Talign<(t))]) ##sum data from previous 20 years
      new<-t(as.matrix(c(t-20,B,agedataall[agedataall$Dataset==z,][1,4:16]))) #Add back in other info from data
      colnames(new)<-colnames(agedataall)[-1]
      
      agedataall20=rbind.data.frame(agedataall20,new) #append
      
    }
    
  }
  
}
agedataall20<-as.data.frame((agedataall20[-1,])) #Remove NA row
agedataall20$Interval<-20 ##all intervals are now 20 years
agedataall20<-as.data.frame(lapply(agedataall20,unlist)) ##For some reason R makes all the dataframe columns a list, so this fixes that. 
agedataall20$CensusOffset<-agedataall20$Time_collected-(agedataall20$Talign+20) ##Calculate new offset (number of years between the end of the interval and when the stands were censused), this is needed for establishment with survival. 
######



######


####Fig. S1####
source("AgeSt_mapping/Age_St_mapping.R") ###Sources other script to create figure
ggsave('FigS1.tiff', plot=map, device="tiff",width = 4, height = 5, dpi = 500)

####Fig. S2####

tiff("FigS2.tiff",width = 9,height=6,units="in", res=300)    


bins20<-seq(1260,2000,by=20)
RelCount<-numeric(length(bins20)-1)
for(i in 2:length(bins20)){
  RelCount[i-1]<-mean(s1$RelB[which(s1$Time>=bins20[i-1] & s1$Time<bins20[i])],na.rm=T) ##mean amont of estblishment in each 20 year interval across all datasets
}


par(mfrow=c(2,1),mar = c(2.5,4.6,2.5,4.5))
barplot(RelCount,space = 0,ylab="Proportion of total establishment",xlab='Year (CE)', ylim=c(0,.18))
at_tick <- seq_len(length(bins20)) 
axis(side = 1, at = at_tick - 1, labels = FALSE)
axis(side = 1, at = seq_along(1:(length(bins20)-1))[c(1,6,11,16)] - 0.5, tick = FALSE, labels = bins20[-length(bins20)][c(1,6,11,16)])

par(new=T)
barplot(RelCount,space = 0,axes=F,ylim=c(0,.18))
box()
fig_label("A",region="plot", cex=2) 

par(new = TRUE)

species<-unique(agedataall$Species)
spcol<-alpha(brewer.pal(n = length(species), name = "Set1"),0.5)
spline<-c(1,2,1,1,2)
plot(-1,-1,xlim=c(1260,2000),ylim=c(0,1),axes = FALSE, bty = "n", xlab = "", ylab = "")
axis(side=4, ylim=c(0,1), col.axis="#984EA3" ,col="#984EA3" )
mtext("Cumulative establishment",side=4,col="#984EA3",line=3) 


abline(h=.5,lty=2,col="#984EA3")

for(i in 1:max(s1$dataset)){ #plot relative popualtion size in each dataset
  tt<-s1$Time[which(s1$dataset==i)]
  pop<-s1$Relpop[which(s1$dataset==i)]
  tt<-tt[which(is.na(pop)==F)]
  pop<-pop[which(is.na(pop)==F)]
  lines(tt+10,pop,lwd=2,col=spcol[which(species==unique(s1$Species[which(s1$dataset==i)])[2])],lty=spline[which(species==unique(s1$Species[which(s1$dataset==i)])[2])])
  
}
legend(1450, 1, legend=species,
       col=alpha(brewer.pal(n = length(species), name = "Set1"),0.8), lwd=2,cex=0.8,lty=spline) 

par(mar = c(5.5,4.6,.5,4.5))

plot(bayespoint$Time,bayespoint$brate,xlim=c(1260,2000),ylim=c(0,2), ylab=expression("Establishment rate (20yrs"^-1~")"), xlab='Year (CE)',col=alpha('black',0.6),pch=19, cex=0.5)
polygon(c(years,rev(years)),c((allmeanupper),rev((allmeanlower))),col = alpha("grey75",.75), border = FALSE)
lines(years,allmean,lwd=4,col='steelblue')

fig_label("B",region="plot", cex=2) 

dev.off()

################
#############
#source('Bayes/Sim.R') ###This will take a while to run 10-30minutes, need to run to create workspace file below
load('Bayes/Simulations.Rdata') ###Simulations already run

tiff("FigS3.tiff",width = 8,height=8,units="in", res=300)    


bmedian<-apply(best,c(1,3),median)
bupper<-apply(best,c(1,3),quantile,probs=.975)
blower<-apply(best,c(1,3),quantile,probs=.025)

brmse<-apply(bdiff,c(2,3),function(x){sqrt(sum(x^2)/length(x))})

coverage<-apply((btrue<bupper)*(btrue>blower),2,mean)




par(mfrow=c(2,2))

plot(c(btrue),apply(best,c(1,3),median), log='xy',xlab='True Value',ylab='Posteior Median Establishement Rate')
abline(0,1)
text(.005,4,"A", cex=2) 

plot(apply(brmse,2,median),type='l',xlab='Simulation Time Interval', ylim=c(0,3),ylab='RMSE')
lines(apply(brmse,2,quantile,probs=.975),type='l',lty=2)
lines(apply(brmse,2,quantile,probs=.025),type='l',lty=2)
fig_label("B",region="plot", cex=2) 

plot(coverage,type='l',xlab='Simulation Time Interval', ylim=c(0,1),ylab='Coverage')
abline(h=0.95,col='tomato',lty=2)
fig_label("C",region="plot", cex=2) 


dev.off()
######################
tiff("FigS4.tiff",width = 8,height=4,units="in", res=300)    


par(mfrow=c(1,2))
plot(c(ntrue),apply(nest,c(1,3),median), log='xy',xlab='True Value',ylab='Posterior Median Population Size')
abline(0,1)
text(15,35000,"A", cex=2) 

plot(c(strue),apply(sest,c(1,3),median), log='xy',xlab='True Value',ylab='Posterior Median Survival Rate', ylim=c(.4,1),xlim=c(.4,1))
abline(0,1)
text(.41,.97,"B", cex=2) 

dev.off()

##############Fig S5##########
tiff("FigS5.tiff",width = 8,height=9, units="in", res=300)    

s1avg<-decadeavg(s1)

par(mfrow=c(3,1),mar = c(4.5,6,2.5,.5))
plot(0,0,ylim=c(0,1),xlim=c(1600,1990), xlab="Year (CE)", ylab=expression("Survival rate (20yrs"^-1~")"))


set.seed(100)
for(i in 1:50){
  sv<-svariable(x=agedataall20,0.95,backward = F)
  sv_mean<-na.omit(decadeavg(sv,var='srate'))
  lines(sv_mean$years, sv_mean$mean,lwd=.2,col='steelblue')
  sv<-svariable(x=agedataall20,0.95,backward = T)
  sv_mean<-na.omit(decadeavg(sv,var='srate'))
  lines(sv_mean$years, sv_mean$mean,lwd=.2,col='steelblue')
}

fig_label("A",region="plot", cex=2) 



plot(0,0,ylim=c(0,.7),xlim=c(1600,1990), ylab=expression("Establishment rate (20yrs"^-1~")"), xlab='Year (CE)')


set.seed(100) # for reproducible results 
for(i in 1:50){
  sv<-svariable(x=agedataall20,0.95,backward = F)
  sv_mean<-na.omit(decadeavg(sv,var='brate'))
  lines(sv_mean$years, sv_mean$mean,lwd=.2, col='steelblue')
  sv<-svariable(x=agedataall20,0.95,backward = T)
  sv_mean<-na.omit(decadeavg(sv,var='brate'))
  lines(sv_mean$years, sv_mean$mean,lwd=.2,col='steelblue')
}
lines(s1avg$years,s1avg$mean, lwd=2,lty=2)
s.90<-sconstant(agedataall20,srate=.90)
s.90avg<-mvavg(s.90,smoothyears = 20)


lines(s.90avg$years,s.90avg$mean, lwd=2,lty=2, col='coral')

fig_label("B",region="plot", cex=2) 




bins20<-seq(1600,2000,by=20)
RelCount<-numeric(length(bins20)-1)
for(i in 2:length(bins20)){
  RelCount[i-1]<-mean(s.90$RelB[which(s.90$Time>=bins20[i-1] & s.90$Time<bins20[i])],na.rm=T) ##mean amont of estblishment in each 20 year interval across all datasets
}

par(mar = c(4.5,6,2.5,.5))


species<-unique(agedataall$Species)
spcol<-alpha(brewer.pal(n = length(species), name = "Set1"),0.5)
spline<-c(1,2,1,1,2)
plot(-1,-1,xlim=c(1600,2000),ylim=c(0,1), xlab = "Year (CE)", ylab = "Relative population density")
fig_label("C",region="plot", cex=2) 

abline(h=.5,lty=2,col="#984EA3")

for(i in 1:max(s.90$dataset)){
  tt<-s.90$Time[which(s.90$dataset==i)]
  pop<-s.90$Relpop[which(s.90$dataset==i)]
  tt<-tt[which(is.na(pop)==F)]
  pop<-pop[which(is.na(pop)==F)]
  lines(tt+10,pop,lwd=2,col=spcol[which(species==unique(s.90$Species[which(s.90$dataset==i)])[2])],lty=spline[which(species==unique(s.90$Species[which(s.90$dataset==i)])[2])])
  
}
legend(1950, 0.4, legend=species,
       col=alpha(brewer.pal(n = length(species), name = "Set1"),0.8), lwd=2,cex=0.8,lty=spline) 

dev.off()



###Fig S6####
###40 year lag
tiff("FigS6.tiff",width = 8,height=3.5,units="in", res=300)    

s1lag<-blag(agedataall20,lag=40)

s1lag.lm <- lm(brate ~ Time, data = s1lag[which(s1lag$Time>=1600),])

s1lag.seg <- segmented(s1lag.lm, 
                      seg.Z = ~ Time, npsi=2)

plot(s1lag$Time,s1lag$brate,xlim=c(1600,1990),ylim=c(0,1.2), ylab=expression("Establishment rate (20yrs"^-1~")"), xlab='Year (CE)',col=alpha('black',0.6),pch=19, cex=0.5)
s1avg<-decadeavg(s1lag)

polygon(c(s1avg$years,rev(s1avg$years)),c((s1avg$upper),rev((s1avg$lower))),col = alpha("grey75",0.75), border = FALSE)
lines(s1avg$years,s1avg$mean,lwd=4,col='steelblue')
plot(s1lag.seg,add=TRUE,link=FALSE,lwd=3, col='coral',lty=2262)
dev.off()


####Figure s7
####Archer Figure ####
tiff("FigS7.tiff",width = 8,height=3.5,units="in", res=300)    
source('ArcherPaper.R') ##Sources separate script to create figure
dev.off()


####Fig S8####
###Digitization errors###
digdata<-read.csv('DigCompData.csv')
tiff("FigS8.tiff",width = 8,height=3.5,units="in", res=300)    

par(mfrow=c(1,2))
plot(digdata$Time,digdata$B1,pch=1, col='Steelblue', xlab="Year (CE)", ylab='Establishment')
points(digdata$Time,digdata$B2,pch=3, col='Coral')
fig_label("A",region="plot", cex=2) 
plot(digdata$B1,digdata$B2,pch=19, col='Steelblue', xlab="Measurement 1", ylab="Measurement 2")
abline(0,1)
fig_label("B",region="plot", cex=2) 
dev.off()



##Figure S9#########
############Bayes vs Determinisitc##########
tiff("FigS9.tiff",width = 4,height=4,units="in", res=300)    

bins20<-seq(1260,2000,by=20)#Calculate mean temp for each 
Detb<-matrix(NA,length(bins20)-1,NPops)
for(i in 1:NPops){
  Site_b<-s1[which(s1$dataset==i),c('Time','brate')]
  Detb[,i]<-Site_b$brate
}
plot(c(sitemeanB),Detb,ylab="Deterministic Estimate",xlab='Bayesian posterior median')
abline(0,1)
text(0.5,6,paste("r=",round(cor(c(sitemeanB),c(Detb),use='pairwise.complete.obs'),2)))
dev.off()


###Code for producing supplemental figures#####
library(boot)
library(RColorBrewer)
library(scales)
library(ggplot2)
library(gridExtra)
library(plyr)
library(plotrix)
library(segmented)


source("Functions.R")
###Read in fit models and data from the Analysis Script###
#load('Fitmodels_Rev.Rdata') ###Modelfits already run


s1<-sconstant(agedataall20,srate=1) 
spcol<-alpha(brewer.pal(n = length(species), name = "Set1"),0.7)

####Fig S1### Data Plot by site
tiff("FigS1.tiff",width = 8,height=12,units="in", res=300)    
par(mfrow=c(8,4),  mar = c(4,4,1,1))
for(p in 1:NPops){
  
  
  lim<-max(s1$B[which(s1$dataset==p)],na.rm=T)
  plot(seq(1260,1980,20),s1$B[which(s1$dataset==p)],ylim=c(0,lim),type='h',col=spcol[match(specieslist,species)][p],lwd=2,lend=1, ylab="Establishment (# of trees)",xlab='Year',cex.lab=0.8)
  fig_label(paste(p,specieslist[p]),region="plot", cex=1.5) 
  
  
  
}

dev.off()


####Fig. S2#### Full Time serie
s1<-s1[which(s1$dataset<30),]

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
spcol<-alpha(brewer.pal(n = length(species), name = "Set1"),0.35)

matplot(years,sitemeanB[,1:29],xlim=c(1260,2000),ylim=c(0,2),pch=19,col=spcol[match(specieslist,species)],lty=spline[match(specieslist,species)], ylab=expression("Establishment rate (20yrs"^-1~")"), xlab='Year (CE)', lwd=2)
polygon(c(years,rev(years)),c((allmeanupper),rev((allmeanlower))),col = alpha("grey75",.75), border = FALSE)
lines(years,allmean,lwd=4,col='steelblue')

fig_label("B",region="plot", cex=2) 

dev.off()


####Fig S3### Establishment rates by site

spcol<-alpha(brewer.pal(n = length(species), name = "Set1"),0.7)

load("Fitmodels_Rev.RData")
tiff("FigS3.tiff",width = 8,height=12,units="in", res=300)    
par(mfrow=c(8,4),  mar = c(4,4,1,1))

for (i in 1:32){
  med<-apply(bout.use[,,i],2,median)
  upper<-apply(bout.use[,,i],2,quantile,probs=.975,na.rm=T)
  lower<-apply(bout.use[,,i],2,quantile,probs=.025,na.rm=T)
  plot(years,med,type='l',lwd=2, ylim=c(0,1.3), col=spcol[match(specieslist,species)][i], ylab=expression("Establishment rate (20yrs"^-1~")"), xlab='Year (CE)')
  lines(years,upper, lty=2, col=spcol[match(specieslist,species)][i])
  lines(years,lower, lty=2,, col=spcol[match(specieslist,species)][i])
  fig_label(paste(i,specieslist[i]),region="plot", cex=1.5) 
  
}
dev.off()









####### Fig S4##### Compare establishment with constant rates
Bval<-array(NA,c(3000,20,29))
Bpredval<-array(NA,c(3000,20,29))
tiff("FigS4.tiff",width = 8,height=12,units="in", res=300)    
par(mfrow=c(8,4),  mar = c(4,4,1,1))
for(p in 1:29){
  
  for (i in 1:3000){
    
    
    if((stime[p,1]+1)<18){miny=18} #min if before 1600
    if ((stime[p,1]+1)>17){miny=(stime[p,1]+1)}
    
    lam<-sout.use[i,post1600,p]+bout.use[i,post1600,p] ###Calculate lambda estimates survival + establishment 
    lampred<-sout.use[i,post1600,p]+mean(bout.use[i,post1600,p][1:13],na.rm=T)###Calculate lambda estimates survival + mean establishment 
    N<-c(1,cumprod(na.omit(lam)))*Nout[i,miny-1,p]
    Npred<-c(1,cumprod(na.omit(lampred)))*Nout[i,miny-1,p]
    
    B=(head(N,-1)*na.omit(bout.use[i,post1600,p]))*rev(cumprod(rev(na.omit(sout.use[i,post1600,p])))) ###Observed Establishment is population size*estab. rate and then the cumulative survival of cohort is calculated to time of observation. 
    Bpred=(head(Npred,-1)*mean(bout.use[i,post1600,p][1:13],na.rm=T))*rev(cumprod(rev(na.omit(sout.use[i,post1600,p])))) ###Same as above, but with mean estab. rate
    Bval[i,(miny:stime[p,2])-17,p]<-B #Relativize to be proportion of total
    Bpredval[i,(miny:stime[p,2])-17,p]<-Bpred
    
    
  } 
  lim<-ifelse(max(apply(Bval[,,p],2,median,na.rm=T),na.rm = T)>max(apply(Bpredval[,,p],2,median,na.rm=T),na.rm = T),max(apply(Bval[,,p],2,median,na.rm=T),na.rm = T),max(apply(Bpredval[,,p],2,median,na.rm=T),na.rm = T))
  plot(seq(1600,1980,20)-4,apply(Bval[,,p],2,median,na.rm=T),ylim=c(0,lim),type='h',col='coral',lwd=2,lend=1, ylab="Establishment (# of trees)",xlab='Year',cex.lab=0.8)
  arrows(x0=seq(1600,1980,20)+4,y0=0,y1=apply(Bpredval[,,p],c(2),median,na.rm=T),length=0,lend=1,col='steelblue',lwd=2)
  fig_label(paste(p,specieslist[p]),region="plot", cex=1.5) 
  
  
}
dev.off()



#############
#source('Bayes/Sim.R') ###This will take a while to run 10-30minutes, need to run to create workspace file below
load('Bayes/Simulations.Rdata') ###Simulations already run

tiff("FigS5.tiff",width = 8,height=8,units="in", res=300)    


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
tiff("FigS6.tiff",width = 8,height=4,units="in", res=300)    


par(mfrow=c(1,2))
plot(c(ntrue),apply(nest,c(1,3),median), log='xy',xlab='True Value',ylab='Posterior Median Population Size')
abline(0,1)
text(15,35000,"A", cex=2) 

plot(c(strue),apply(sest,c(1,3),median), log='xy',xlab='True Value',ylab='Posterior Median Survival Rate', ylim=c(.4,1),xlim=c(.4,1))
abline(0,1)
text(.41,.97,"B", cex=2) 

dev.off()

######Fig. S7####


tiff("FigS7.tiff",width = 8,height=3.5,units="in", res=300)    
source('MastPaper.R') ##Sources separate script to create figure
dev.off()


####

################

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

######Fig S9: Climate mapping####


FIAplots<-read.csv("FIA_Clim_Normals.csv")
STRUCplots<-read.csv("PJstructure_Clim_Normals_all_months.csv")

plot(FIAplots$tmean,FIAplots$map, cex=0.01)
points(STRUCplots$tmean,STRUCplots$map, col='red', pch=19)

STRUCplots$Species<-factor(STRUCplots$Species,levels=c("PiMo-JuOs", "JuOc","PiMo","PiEd","JuOs"))

Cmap <- ggplot(FIAplots, mapping = aes(x = tmean, y = map))   

ClimateMap<-Cmap + 
  geom_density_2d(bins=10,show.legend=F)+geom_point(alpha=0.1) +geom_jitter(data=STRUCplots,width=.2, height=3,aes(x = tmean, y = map,col=Species,size=1.5))+guides(size = "none")+ylim(100,750)+xlim(4,17)+scale_colour_manual(values=spcol)+  
  geom_point(aes(x=mean(FIAplots$tmean), y=mean(FIAplots$map)), colour="blue", shape=4, size=3)+
  xlab("MAT (\u00B0C)")+ylab("MAP (mm)")+theme_bw()

ggsave('FigS9.tiff', plot=ClimateMap, device="tiff",width = 6, height = 5, dpi = 500)

########

####Fig S10###
##Made by rerunning Fig 2 code with models fit with Uniform(1,1000) prior


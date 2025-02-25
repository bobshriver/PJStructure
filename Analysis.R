###Set Working directory to file location###
library(boot)
library(RColorBrewer)
library(scales)
library(ggplot2)
library(ggtext)
library(gridExtra)
library(plyr)
library(plotrix)
library(segmented)


source("Functions.R") ###pull in script with functions
agedataall<-read.csv("AgeSt2.csv") ####Load Data
agedataall<-agedataall[-which(agedataall$Census_keep=='n'),]###Remove data after 1940-1960 that are not complete census. Also remove trailing intervals that were not complete (e.g. last census was in 2003 and included trees from 2000-2010)
agedataall$B=agedataall$B*agedataall$DensityCon ###Convert back to raw counts from density
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




s1<-sconstant(agedataall20,srate=1) #### Calculate establishment rates assuming 100% survival

NPops<-32
Nyears<-400
bins20<-seq(1600,2000,by=20)#Calculate mean temp for each 

source('bayes/SiteDataPrep.R')###To fit Bayesian model for establishment rate. This takes about 10 minutes to run on my mac
#Save model fits
#save.image(file='Fitmodels_Rev.RData')

specieslist<-numeric(NPops)
for(i in 1:NPops){
specieslist[i]<-na.exclude(s1[which(s1$dataset==i),c('Species')])[1]
}

preSet<-18:30



sitemeanB<-apply(bout.use,c(2,3),median,na.rm=T)
sitemeanN<-apply(Nout,c(2,3),median,na.rm=T)
sitemeanS<-apply(sout.use,c(2,3),median,na.rm=T)
###Averages across sites for pre-1800 populations
bayespoint<-cbind.data.frame(rep(years,29), c(sitemeanB[,1:29]))
colnames(bayespoint)<-c("Time",'brate')
allmean<-apply(apply(bout.use[,,1:29],c(1,2),mean,na.rm=T),2,median,na.rm=T)
allmeanupper<-apply(apply(bout.use[,,1:29],c(1,2),mean.Date,na.rm=T),2,quantile,prob=c(0.975),na.rm=T)
allmeanlower<-apply(apply(bout.use[,,1:29],c(1,2),mean,na.rm=T),2,quantile,prob=c(0.025),na.rm=T)

#Segmented Regression
s1.lm <- lm(brate ~ Time, data = bayespoint[which(bayespoint$Time>=1600),])
#one break point
s1.seg <- segmented(s1.lm, 
                    seg.Z = ~ Time, npsi=1)
slope(s1.seg)
s1.seg$psi
AIC(s1.seg)

#two break point
s2.seg <- segmented(s1.lm, 
                    seg.Z = ~ Time, npsi=2)
slope(s2.seg)
s2.seg$psi
AIC(s2.seg)


s3.seg <- segmented(s1.lm, 
                    seg.Z = ~ Time, npsi=3)
slope(s3.seg)
s3.seg$psi
AIC(s3.seg)


s4.seg <- segmented(s1.lm, 
                    seg.Z = ~ Time, npsi=4)
slope(s4.seg)
s4.seg$psi
AIC(s4.seg)





###Deterministic Average
s1avg<-decadeavg(s1)
s1<-s1[which(s1$dataset<30),]#only pre-1800 populations
  

####Figure 2##########
tiff("Fig2.tiff",width = 8,height=6.5,units="in", res=300)    

bins20<-seq(1600,2000,by=20)
RelCount<-numeric(length(bins20)-1)
for(i in 2:length(bins20)){
RelCount[i-1]<-mean(s1$RelB[which(s1$Time>=bins20[i-1] & s1$Time<bins20[i])],na.rm=T) ##mean amount of establishment in each 20 year interval across all presettlment data sets
}


par(mfrow=c(2,1),mar = c(2.5,4.6,2.5,4.5))
barplot(RelCount,space = 0,ylab="Proportion of total establishment",xlab='Year', ylim=c(0,.18))
at_tick <- seq_len(length(bins20)) 
axis(side = 1, at = at_tick - 1, labels = FALSE)
axis(side = 1, at = seq_along(1:(length(bins20)-1))[c(1,6,11,16)] - 0.5, tick = FALSE, labels = bins20[-length(bins20)][c(1,6,11,16)])


par(new=T)
barplot(RelCount,space = 0,axes=F,ylim=c(0,.18))
box()
fig_label("B",region="plot", cex=2) 

par(new = TRUE)

species<-unique(agedataall$Species)
spcol<-alpha(brewer.pal(n = length(species), name = "Set1"),0.5)
spline<-c(1,2,1,1,2)
plot(-1,-1,xlim=c(1600,2000),ylim=c(0,1),axes = FALSE, bty = "n", xlab = "", ylab = "")
axis(side=4, ylim=c(0,1), col.axis="#984EA3" ,col="#984EA3" )
mtext("Cumulative establishment",side=4,col="#984EA3",line=3) 


#abline(h=.5,lty=2,col="#984EA3")
splist<-numeric(29)
for(i in 1:max(s1$dataset)){ #plot cumulative est. in each dataset
  tt<-s1$Time[which(s1$dataset==i)]
  pop<-s1$Relpop[which(s1$dataset==i)]
  tt<-tt[which(is.na(pop)==F)]
  pop<-pop[which(is.na(pop)==F)]
  lines(tt+10,pop,lwd=2,col=spcol[which(species==unique(s1$Species[which(s1$dataset==i)])[2])],lty=spline[which(species==unique(s1$Species[which(s1$dataset==i)])[2])])
}
legend(1650, 1, legend=species,
       col=alpha(brewer.pal(n = length(species), name = "Set1"),0.8), lwd=2,cex=0.8,lty=spline) 

par(mar = c(5.5,4.6,.5,4.5))
###Plot Establishment rates from bayesian estimation
spcol<-alpha(brewer.pal(n = length(species), name = "Set1"),0.35)

#plot pre-1800 populations
matplot(years,sitemeanB[,1:29],xlim=c(1600,2000),ylim=c(0,.7),pch=19,col=spcol[match(specieslist,species)],lty=spline[match(specieslist,species)], ylab=expression("Establishment rate (20yrs"^-1~")"), xlab='Year (CE)', lwd=2)

polygon(c(years,rev(years)),c((allmeanupper),rev((allmeanlower))),col = alpha("grey75",.8), border = FALSE)
plot.segmented(s2.seg,add=TRUE,lwd=4, col='coral',lty=2262,psi.lines=TRUE)
lines(years,allmean,lwd=4,col='steelblue')

fig_label("C",region="plot", cex=2) 

dev.off()

source("AgeSt_mapping/Age_St_mapping.R") ###Sources other script to create figure
ggsave('Fig2A.tiff', plot=map, device="tiff",width = 4, height = 5, dpi = 500)


####Figure 3#####
#Containers
bratio<-array(NA,c(3000,20,29))
Bval<-Bvaltot<-array(NA,c(3000,20,29))
Bpredval<-Bpredvaltot<-array(NA,c(3000,20,29))
###Loop over chains and populations
for (i in 1:3000){
  
  
  for(p in 1:29){
    if((stime[p,1]+1)<18){miny=18} #min if before 1600
    if ((stime[p,1]+1)>17){miny=(stime[p,1]+1)}
    
    lam<-sout.use[i,post1600,p]+bout.use[i,post1600,p] ###Calculate lambda estimates survival + establishment 
    lampred<-sout.use[i,post1600,p]+mean(bout.use[i,post1600,p][1:13],na.rm=T)###Calculate lambda estimates survival + mean establishment 
    N<-c(1,cumprod(na.omit(lam)))
    Npred<-c(1,cumprod(na.omit(lampred)))
    
    B=(head(N,-1)*na.omit(bout.use[i,post1600,p]))*rev(cumprod(rev(na.omit(sout.use[i,post1600,p])))) ###Observed Establishment is population size*estab. rate and then the cumulative survival of cohort is calculated to time of observation. 
    Bpred=(head(Npred,-1)*mean(bout.use[i,post1600,p][1:13],na.rm=T))*rev(cumprod(rev(na.omit(sout.use[i,post1600,p])))) ###Same as above, but with mean estab. rate
    Bval[i,(miny:stime[p,2])-17,p]<-B/sum(B) #Relativize to be proportion of total
    Bpredval[i,(miny:stime[p,2])-17,p]<-Bpred/sum(Bpred)
    Bvaltot[i,(miny:stime[p,2])-17,p]<-B 
    Bpredvaltot[i,(miny:stime[p,2])-17,p]<-Bpred
    
    bratio[i,(miny:stime[p,2])-17,p]<-Bpred/B ###Constant/Changing
    
    
  }  
  
}
ratiomedian<-(apply((apply((bratio),c(1,2),median,na.rm=T)),2,median,na.rm=T))
ratiomedianupper<-apply(((apply((bratio),c(1,2),median,na.rm=T))),2,quantile,prob=c(0.975),na.rm=T)
ratiomedianlower<-apply(((apply((bratio),c(1,2),median,na.rm=T))),2,quantile,prob=c(0.025),na.rm=T)


popcor<-numeric(29)
for(i in 1:29){popcor[i]<-cor(apply(Bpredvaltot,c(2,3),median,na.rm=T)[,i], apply(Bvaltot,c(2,3),median,na.rm=T)[,i], use='pairwise.complete.obs')

}
cortab<-table(specieslist[1:29],cut(popcor^2,breaks=c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1)))
colnames(cortab)[1:9]<-c('0.1-0.2','0.2-0.3','0.3-0.4','0.4-0.5','0.5-0.6','0.6-0.7','0.7-0.8','0.8-0.9','0.9-1')


tiff("Fig3.tiff",width = 8,height=6.5,units="in", res=300)    

par(mfrow=c(2,1),mar = c(1.5,5,2.5,2.5))


plot(seq(1600,1980,20)-4,apply(apply(Bval,c(2,3),median,na.rm=T),1,mean,na.rm=T),type='h',col='coral',lwd=11,lend=1, ylim=c(0,.25), ylab="Proportion of total establishment",xlab='Year')
arrows(x0=seq(1600,1980,20)+4,y0=0,y1=apply(apply(Bpredval,c(2,3),median,na.rm=T),1,mean,na.rm=T),length=0,lend=1,col='steelblue',lwd=11)
legend(1590, .1, legend=c("Changing rate","Constant rate"),
       col=c('coral','steel blue'), cex=0.75,pch=15) 

legend(1816, .258, legend=species,
       col=alpha(brewer.pal(n = length(species), name = "Set1"),0.8), cex=0.45,pt.cex=0.9,pch=15) 
fig_label("A",region="plot", cex=2) 

text(x=1780,y=.08,"Proportion of variance explained",cex=0.75)
text(x=1715,y=.19,"# of populations",cex=0.75,srt=90)

par(mar = c(4.5,5,1.5,2.5))
matplot(seq(1600,1980,20),apply(bratio,c(2,3),median,na.rm=T),log='y',ylim=c(.2,5),type='l',col=spcol[match(specieslist,species)],lty=spline[match(specieslist,species)],lwd=1.5, ylab='Ratio of establishement', xlab='Year (CE)', cex.lab=0.85)
polygon(c(seq(1600,1980,20),rev(seq(1600,1980,20))),c((ratiomedianupper),rev((ratiomedianlower))),col = alpha("grey75",.75), border = FALSE)
lines(seq(1600,1980,20),ratiomedian,lwd=3)

abline(h=1,lty=2,lwd=2)
legend(1750, .9, legend=species,
       col=alpha(brewer.pal(n = length(species), name = "Set1"),0.8), lwd=2,cex=0.8,lty=spline) 
text(1595,4.5,"C", cex=2) 


par(fig = c(0.45,.6,0.765,.915), new = TRUE,mar = c(0,0,0,0),cex.lab=0.75, cex.axis=0.75,mgp=c(3,1,0))

#plot(apply(Bval,c(2,3),median,na.rm=T),apply(Bpredval,c(2,3),median,na.rm=T),col="black",cex=.1,pch=19)
barplot(cortab,col=spcol[c(2,5,4,3,1)], las = 2, cex.lab=.75, cex.axis =.75, cex.names = .6)
axis(side = 1, at=seq(0.1,1.1*10,length.out=10), labels = FALSE)

fig_label("B",region="plot", cex=1.5) 
#abline(0,1)
#text(.4,.1,paste("r=",round(cor(c(apply(Bval,c(2,3),median,na.rm=T)),c(apply(Bpredval,c(2,3),median,na.rm=T)),use='pairwise.complete.obs'),2)),cex=0.75)

dev.off()

#Figure 4  ###################################################
tiff("Fig4.tiff",width = 8,height=6, units="in", res=300)    

par(mfrow=c(2,1),mar = c(1.5,5,2.5,2.5))
matplot(seq(1600,1980,by=20),t(sout.use[,18:37,15]),type='l',lty=1,lwd=.1,ylab="Survival Rate")
fig_label("A",region="plot", cex=2) 
par(mar = c(4.5,5,1.5,2.5))

matplot(seq(1600,1980,by=20),t(bout.use[,18:37,15]),type='l',lty=1,lwd=.1, ylim=c(0,1.5), xlab='Year (CE)',ylab="Establishment Rate")
fig_label("B",region="plot", cex=2) 
text(x=1915,y=.31,"Survival Rate",cex=0.75)
text(x=1845,y=1.1,"Estab. Rate",cex=0.75,srt=90)


par(fig = c(0.7,.85,0.29,.44), new = TRUE,mar = c(0,0,0,0),cex.lab=0.75, cex.axis=0.75,mgp=c(3,1,0))

plot(sout.use,bout.use,col="black",cex=.01,pch=19,ylim=c(0,8))
text(.6,6,paste("r=",round(cor(c(sout.use),c(bout.use),use='pairwise.complete.obs'),2)),cex=0.75)
fig_label("C",region="plot", cex=1) 

dev.off()



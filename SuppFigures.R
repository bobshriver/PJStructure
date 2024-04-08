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
agedataall20<-matrix(NA,1,14)
colnames(agedataall20)<-colnames(agedataall)[-1]
for (z in 1:length(unique(agedataall$Dataset))){
  
  if(agedataall$Interval[agedataall$Dataset==z][1]==20){
    
    agedataall20=rbind.data.frame(agedataall20,agedataall[agedataall$Dataset==z,-1])
    
  }
  
  if(agedataall$Interval[agedataall$Dataset==z][1]==10){
    duse<-agedataall[agedataall$Dataset==z,]
    tmin<-min(duse$Talign)
    tmax<-max(duse$Talign)
    if(length(which(tmin==bins))==0){tmin=tmin-10}
    if(length(which(tmax==bins))>0){tmax=tmax-10}
    zbins<-seq(tmin,tmax+10,20)
    
    for(t in zbins[-1]){
      B<-sum(duse$B[which(duse$Talign>=(t-20) & duse$Talign<(t))])
      new<-t(as.matrix(c(t-20,B,agedataall[agedataall$Dataset==z,][1,4:15])))
      colnames(new)<-colnames(agedataall)[-1]
      
      agedataall20=rbind.data.frame(agedataall20,new)
      
    }
    
  }
  
}
agedataall20<-as.data.frame((agedataall20[-1,]))
agedataall20$Interval<-20
agedataall20<-as.data.frame(lapply(agedataall20,unlist))
agedataall20$CensusOffset<-agedataall20$Time_collected-(agedataall20$Talign+20)
######


####Fig. S1####
source("AgeSt_mapping/Age_St_mapping.R") ###Sources other script to create figure
ggsave('FigS1.tiff', plot=map, device="tiff")

####Fig. S2####


tiff("FigS2.tiff",width = 8,height=6.5,units="in", res=300)    

s1<-sconstant(agedataall20,srate=1)


bins20<-seq(1260,2000,by=20)
RelCount<-numeric(length(bins20)-1)
for(i in 2:length(bins20)){
  RelCount[i-1]<-mean(s1$RelB[which(s1$Time>=bins20[i-1] & s1$Time<bins20[i])],na.rm=T)
}

par(mfrow=c(2,1),mar = c(2.5,4.6,2.5,4.5))

barplot(RelCount,space = 0,ylab="Proportion of total establishment",xlab='Year')
at_tick <- seq_len(length(bins20)) 
axis(side = 1, at = at_tick - 1, labels = FALSE)
axis(side = 1, at = seq_along(1:(length(bins20)-1)) - 0.5, tick = FALSE, labels = bins20[-length(bins20)])

rect(-100, -1000, 29.5, 1000, col=adjustcolor("grey",0.3 ), border='white',lwd=2)
rect(29.5, -1000, 2050, 1000, col=adjustcolor("coral", 0.2), border='white',lwd=2)
par(new=T)
barplot(RelCount,space = 0,axes=F)
box()

par(new = TRUE)


species<-unique(agedataall$Species)
spcol<-alpha(brewer.pal(n = length(species), name = "Set1"),0.5)
spline<-c(1,2,1,1,2)
plot(-1,-1,xlim=c(1260,2000),ylim=c(0,1),axes = FALSE, bty = "n", xlab = "", ylab = "")
axis(side=4, ylim=c(0,1), col.axis="#984EA3" ,col="#984EA3" )
mtext("Relative population density",side=4,col="#984EA3",line=3) 

for(i in 1:max(s1$dataset)){
  tt<-s1$Time[which(s1$dataset==i)]
  pop<-s1$Relpop[which(s1$dataset==i)]
  tt<-tt[which(is.na(pop)==F)]
  pop<-pop[which(is.na(pop)==F)]
  lines(tt+10,pop,lwd=2,col=spcol[which(species==unique(s1$Species[which(s1$dataset==i)])[2])],lty=spline[which(species==unique(s1$Species[which(s1$dataset==i)])[2])])
  
}
legend(1400, 1, legend=species,
       col=alpha(brewer.pal(n = length(species), name = "Set1"),0.8), cex=0.8,lty=spline)   
fig_label("A",region="plot", cex=2) 


plot(s1$Time,s1$brate,xlim=c(1260,1990),ylim=c(0,8.5), ylab=expression("Establishment rate (20yrs"^-1~")"), xlab='Year',col=alpha('black',0.6),pch=19, cex=0.5)
s1avg<-decadeavg(s1)
rect(0, -1000, 1850, 1000, col=adjustcolor("grey",0.3 ), border='white',lwd=2)
rect(1850, -1000, 2050, 1000, col=adjustcolor("coral", 0.2), border='white',lwd=2)
polygon(c(s1avg$years,rev(s1avg$years)),c((s1avg$upper),rev((s1avg$lower))),col = alpha("grey75",1), border = FALSE)
lines(s1avg$years,s1avg$mean,lwd=4,col='steelblue')
text(1750,6,substitute(paste(bold('Pre E-A Settlement'))),cex=0.75)
text(1930,6,substitute(paste(bold('Post E-A Settlement'))),cex=0.75)
fig_label("B",region="plot", cex=2) 
dev.off()

################


##########Figure S3#########


s.90<-sconstant(agedataall20,srate=.90)
s.90.lm <- lm(brate ~ Time, data = s.90[which(s.90$Time>1600),])

s.90.seg <- segmented(s.90.lm, 
                    seg.Z = ~ Time, npsi=2)


bins20<-seq(1600,2000,by=20)
RelCount<-numeric(length(bins20)-1)
for(i in 2:length(bins20)){
  RelCount[i-1]<-mean(s.90$RelB[which(s.90$Time>=bins20[i-1] & s.90$Time<bins20[i])],na.rm=T) ##mean amont of estblishment in each 20 year interval across all datasets
}

tiff("FigS3.tiff",width = 8,height=6.5,units="in", res=300)    

par(mfrow=c(2,1),mar = c(2.5,4.6,2.5,4.5))
barplot(RelCount,space = 0,ylab="Proportion of total establishment",xlab='Year')
at_tick <- seq_len(length(bins20)) 
axis(side = 1, at = at_tick - 1, labels = FALSE)
axis(side = 1, at = seq_along(1:(length(bins20)-1))[c(1,6,11,16)] - 0.5, tick = FALSE, labels = bins20[-length(bins20)][c(1,6,11,16)])
rect(-10, -1000, 12.5, 1000, col=adjustcolor("grey",0.3 ), border='white',lwd=2)
rect(12.5, -1000, 2050, 1000, col=adjustcolor("coral", 0.2), border='white',lwd=2)

par(new=T)
barplot(RelCount,space = 0,axes=F)
box()
fig_label("A",region="plot", cex=2) 

par(new = TRUE)

species<-unique(agedataall$Species)
spcol<-alpha(brewer.pal(n = length(species), name = "Set1"),0.5)
spline<-c(1,2,1,1,2)
plot(-1,-1,xlim=c(1600,2000),ylim=c(0,1),axes = FALSE, bty = "n", xlab = "", ylab = "")
axis(side=4, ylim=c(0,1), col.axis="#984EA3" ,col="#984EA3" )
mtext("Relative population density",side=4,col="#984EA3",line=3) 


abline(h=.5,lty=2,col="#984EA3")

for(i in 1:max(s.90$dataset)){
  tt<-s.90$Time[which(s.90$dataset==i)]
  pop<-s.90$Relpop[which(s.90$dataset==i)]
  tt<-tt[which(is.na(pop)==F)]
  pop<-pop[which(is.na(pop)==F)]
  lines(tt+10,pop,lwd=2,col=spcol[which(species==unique(s.90$Species[which(s.90$dataset==i)])[2])],lty=spline[which(species==unique(s.90$Species[which(s.90$dataset==i)])[2])])
  
}
legend(1650, 1, legend=species,
       col=alpha(brewer.pal(n = length(species), name = "Set1"),0.8), lwd=2,cex=0.8,lty=spline) 

par(mar = c(5.5,4.6,.5,4.5))

plot(s.90$Time,s.90$brate,xlim=c(1600,2000),ylim=c(0,0.8), ylab=expression("Establishment rate (20yrs"^-1~")"), xlab='Year (CE)',col=alpha('black',0.6),pch=19, cex=0.5)
s.90avg<-decadeavg(s.90)
rect(0, -1000, 1850, 1000, col=adjustcolor("grey",0.3 ), border='white',lwd=2)
rect(1850, -1000, 2050, 1000, col=adjustcolor("coral", 0.2), border='white',lwd=2)
polygon(c(s.90avg$years,rev(s.90avg$years)),c((s.90avg$upper),rev((s.90avg$lower))),col = alpha("grey75",.75), border = FALSE)
lines(s.90avg$years,s.90avg$mean,lwd=4,col='steelblue')
text(1750,0.75,substitute(paste(bold('Pre E-A Settlement'))))
text(1930,0.75,substitute(paste(bold('Post E-A Settlement'))))
plot(s.90.seg,add=TRUE,link=FALSE,lwd=3, col='coral',lty=2262)
fig_label("B",region="plot", cex=2) 

dev.off()




###Fig S4####
###40 year lag
tiff("FigS4.tiff",width = 8,height=3.5,units="in", res=300)    

s1lag<-blag(agedataall20,lag=40)

s1lag.lm <- lm(brate ~ Time, data = s1lag[which(s1lag$Time>=1600),])

s1lag.seg <- segmented(s1lag.lm, 
                      seg.Z = ~ Time, npsi=2)

plot(s1lag$Time,s1lag$brate,xlim=c(1600,1990),ylim=c(0,1.2), ylab=expression("Establishment rate (20yrs"^-1~")"), xlab='Year',col=alpha('black',0.6),pch=19, cex=0.5)
s1avg<-decadeavg(s1lag)
rect(0, -1000, 1850, 1000, col=adjustcolor("grey",0.3 ), border='white',lwd=2)
rect(1850, -1000, 2050, 1000, col=adjustcolor("coral", 0.2), border='white',lwd=2)
polygon(c(s1avg$years,rev(s1avg$years)),c((s1avg$upper),rev((s1avg$lower))),col = alpha("grey75",0.75), border = FALSE)
lines(s1avg$years,s1avg$mean,lwd=4,col='steelblue')
text(1750,3,substitute(paste(bold('Pre E-A Settlement'))))
text(1930,3,substitute(paste(bold('Post E-A Settlement'))))
plot(s1lag.seg,add=TRUE,link=FALSE,lwd=3, col='coral',lty=2262)
dev.off()
####Fig S5#######

yrs<-match(seq(1840,1940,by=20),s1avg$years)
brate1<-c(s1avg$mean[yrs][1],rep(s1avg$mean[which(s1avg$years==1860)],length(yrs)-1))+1
brateall<-s1avg$mean[yrs]+1
Btotal<-diff(c(0.1,0.1*cumprod(brateall)))[-1]
Bexpg<-diff(c(0.1,0.1*cumprod(brate1)))[-1]

yrs<-yrs[-1]
Bbrate<-Btotal-Bexpg
PBbrate<-Bbrate/Btotal
PBexpg<-Bexpg/Btotal
Percent<-c(PBexpg,PBbrate)
Pdata<-as.data.frame(Percent)
Pdata$Mechanism<-c(rep("Multiplicative growth", length(yrs)),rep("Changing estab. rate", length(yrs)))
Pdata$Year<-rep(s1avg$years[yrs],2)
Pdata$Establishment<-c(Bexpg,Bbrate)

Percent<-c(sum(Bexpg)/sum(Btotal),sum(Bbrate)/sum(Btotal))
Tdata<-as.data.frame(Percent)
Tdata$Mechanism<-c("Multiplicative growth","Changing estab. rate")
Tdata$Year<-rep(1,2)


p1<-ggplot(Pdata, aes(fill=Mechanism, y=Establishment, x=Year)) + 
  geom_bar(position="stack", stat="identity")+scale_fill_manual(values=c(alpha('coral',0.7),alpha('steelblue',0.7)))+scale_x_continuous(breaks=c(1860,1880, 1900, 1920, 1940))+ylab("Establishment")+xlab("Year (CE)")+theme_bw()+geom_text(x = 1855, y = .7
                                                                                                                                                                                                                                            , width = unit(.8, "inch")
                                                                                                                                                                                                                                            , label = "A"
                                                                                                                                                                                                                                            
                                                                                                                                                                                                                                            , color = "black"
                                                                                                                                                                                                                                            , size=6)

# Plot
p2<-ggplot(Pdata, aes(x=Year, y=Percent, fill=Mechanism)) + 
  geom_area(alpha=0.6 , size=1, colour="white")+scale_fill_manual(values=c('coral','steelblue'))+theme_bw()+theme(legend.position="none")+ylab("Proportional contribution")+xlab("Year (CE)")+geom_text(x = 1865, y = 1
                                                                                                                                                                                                        , width = unit(1, "inch")
                                                                                                                                                                                                        , label = "B"
                                                                                                                                                                                                        , color = "black"
                                                                                                                                                                                                        , size=6 )


p3<-ggplot(Tdata, aes(x=Year, y=Percent, fill=Mechanism)) + scale_fill_manual(values=c(alpha('coral',0.7),alpha('steelblue',0.7)))+geom_col()+ xlab("Total contribution")+ylab('')+theme_bw()+theme(axis.text.x=element_text(color='white'),axis.ticks.x=element_blank(),
                                                                                                                                                                                                    legend.position="none")+geom_text(x = 0.6, y = 1
                                                                                                                                                                                                                                      , width = unit(.8, "inch")
                                                                                                                                                                                                                                      , label = "C"
                                                                                                                                                                                                                                      , color = "black"
                                                                                                                                                                                                                                      , size=6 )
F3<-grid.arrange(p1,p2,p3,
                 layout_matrix = rbind(c(1, 1, 1, 1),
                                       c(1, 1 ,1, 1),
                                       c(2,2 ,2 ,3),
                                       c(2,2 ,2, 3))
)

ggsave('FigS5.tiff', plot=F3, device="tiff")

####Fig S6####
#####Precip Figure####
tiff("FigS6.tiff",width = 9,height=2.9,units="in", res=300)    

s1avg<-decadeavg(s1)

precip<-read.csv("PrecipRecords.csv",header = T)#import reconstructed precip data
bins20<-seq(1600,1980,by=20)
Precip<-numeric(length(bins20)-1)
for(i in 2:length(bins20)){
  Precip[i-1]<-mean(precip$Precip[which(precip$Year>=bins20[i-1] & precip$Year<bins20[i])],na.rm=T)
}

layout(mat=matrix(c(1,1,1,1,2,2),1,6))
par(mar = c(4.5,4.6,4,4.5))

plot(seq(1600,1960,by=20),s1avg$mean[18:36], type='l', ylim=c(0,0.8),lwd=3,col='steelblue', ylab=expression("Establishment rate (20yrs"^-1~")"), xlab='Year (CE)')
polygon(c(s1avg$years[18:37],rev(s1avg$years[18:37])),c((s1avg$upper[18:37]),rev((s1avg$lower[18:37]))),col = alpha("grey75",.75), border = FALSE)
lines(seq(1600,1960,by=20),s1avg$mean[18:36], type='l', ylim=c(0,0.8),lwd=3,col='steelblue')

par(new=T)

plot(seq(1600,1960,by=20),Precip, type='l',col='coral',lwd=3,axes = FALSE, bty = "n", xlab = "", ylab = "",ylim=c(22,33))
axis(side=4, col.axis="coral" ,col="coral" )
mtext(expression("Mean annual precipitation (cm 20yrs"^-1~")"),side=4,col="coral",line=2.5,cex=0.7) 
fig_label("A",region="plot", cex=2) 

cor(Precip,s1avg$mean[18:36])
Pmod<-lm(s1avg$mean[18:36]~Precip)
summary(Pmod)
newdata=as.data.frame(seq(14.5,16.5,length=20))
colnames(newdata)<-'Precip'
prednew<-predict.lm(Pmod, newdata=newdata,se.fit=T)

plot(Precip,s1avg$mean[18:36], pch=19,ylab=expression("Establishment rate (20yrs"^-1~")"),xlab=expression("Mean annual precipitation (cm 20yrs"^-1~")"))

text(26,0.47,labels=paste("Correlation=",round(cor(Precip,s1avg$mean[18:36]),2)))
fig_label("B",region="plot", cex=2) 

dev.off()
####Figure s7
####Archer Figure ####
tiff("FigS7.tiff",width = 8,height=3.5,units="in", res=300)    
source('ArcherPaper.R') ##Sources seperate script to create figure
dev.off()
####Fig s8
#####Mast et al figure#####
tiff("FigS8.tiff",width = 8,height=3.5,units="in", res=300)    
source('MastPaper.R') ##Sources seperate script to create figure
dev.off()


####Fig S9####
###Digitization errors###
digdata<-read.csv('DigCompData.csv')
tiff("FigS9.tiff",width = 8,height=3.5,units="in", res=300)    

par(mfrow=c(1,2))
plot(digdata$Time,digdata$B1,pch=1, col='Steelblue', xlab="Year", ylab='Establishment')
points(digdata$Time,digdata$B2,pch=3, col='Coral')
fig_label("A",region="plot", cex=2) 
plot(digdata$B1,digdata$B2,pch=19, col='Steelblue', xlab="Measurement 1", ylab="Measurement 2")
abline(0,1)
fig_label("B",region="plot", cex=2) 
dev.off()

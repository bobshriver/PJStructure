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

#####Aggregate data to 20 year intervals#####



bins<-seq(1260,2000,by=20)
agedataall20<-matrix(NA,1,14) ###Matrix that data will be appended to
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
      new<-t(as.matrix(c(t-20,B,agedataall[agedataall$Dataset==z,][1,4:15]))) #Add back in other info from data
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

#Segmented Regression
s1.lm <- lm(brate ~ Time, data = s1[which(s1$Time>=1600),])
#one break point
s1.seg <- segmented(s1.lm, 
                    seg.Z = ~ Time, npsi=1)
slope(s1.seg)
s1.seg$psi

#two break points
s1.seg.2 <- segmented(s1.lm, 
                    seg.Z = ~ Time, npsi=2)
slope(s1.seg.2)
s1.seg.2$psi

#three break points
s1.seg.3 <- segmented(s1.lm, 
                      seg.Z = ~ Time, npsi=3)
slope(s1.seg.3)
s1.seg.3$psi


#four break points# Beyond 4 does not converge
s1.seg.4 <- segmented(s1.lm, 
                      seg.Z = ~ Time, npsi=4)
slope(s1.seg.4)
s1.seg.4$psi

s1avg<-decadeavg(s1)


####Figure 2##########
tiff("Fig2.tiff",width = 8,height=6.5,units="in", res=300)    


bins20<-seq(1600,2000,by=20)
RelCount<-numeric(length(bins20)-1)
for(i in 2:length(bins20)){
RelCount[i-1]<-mean(s1$RelB[which(s1$Time>=bins20[i-1] & s1$Time<bins20[i])],na.rm=T) ##mean amont of estblishment in each 20 year interval across all datasets
}


par(mfrow=c(2,1),mar = c(2.5,4.6,2.5,4.5))
barplot(RelCount,space = 0,ylab="Proportion of total establishment",xlab='Year')
at_tick <- seq_len(length(bins20)) 
axis(side = 1, at = at_tick - 1, labels = FALSE)
axis(side = 1, at = seq_along(1:(length(bins20)-1))[c(1,6,11,16)] - 0.5, tick = FALSE, labels = bins20[-length(bins20)][c(1,6,11,16)])
rect(0, -1000, 12.5, 1000, col=adjustcolor("grey",0.3 ), border='white',lwd=2)
rect(12.5, -1000, 2050, 1000, col=adjustcolor("coral", 0.2), border='white',lwd=2)

par(new=T)
barplot(RelCount,space = 0,axes=F)
box()
fig_label("A",region="plot", cex=2) 

par(new = TRUE)

species<-unique(agedataall$Species)
spcol<-alpha(brewer.pal(n = length(species), name = "Set1"),0.5)
plot(-1,-1,xlim=c(1600,2000),ylim=c(0,1),axes = FALSE, bty = "n", xlab = "", ylab = "")
axis(side=4, ylim=c(0,1), col.axis="#984EA3" ,col="#984EA3" )
mtext("Cumulative establishment",side=4,col="#984EA3",line=3) 


abline(h=.5,lty=2,col="#984EA3")

for(i in 1:max(s1$dataset)){ #plot relative popualtion size in each dataset
  tt<-s1$Time[which(s1$dataset==i)]
  pop<-s1$Relpop[which(s1$dataset==i)]
  tt<-tt[which(is.na(pop)==F)]
  pop<-pop[which(is.na(pop)==F)]
  lines(tt+10,pop,lwd=2,col=spcol[which(species==unique(s1$Species[which(s1$dataset==i)])[2])])
  
}
legend(1650, 1, legend=species,
       col=alpha(brewer.pal(n = length(species), name = "Set1"),0.8), lwd=2,cex=0.8,lty=1) 

par(mar = c(5.5,4.6,.5,4.5))

plot(s1$Time,s1$brate,xlim=c(1600,2000),ylim=c(0,0.8), ylab=expression("Establishment rate (20yrs"^-1~")"), xlab='Year (CE)',col=alpha('black',0.6),pch=19, cex=0.5)
rect(0, -1000, 1850, 1000, col=adjustcolor("grey",0.3 ), border='white',lwd=2)
rect(1850, -1000, 2050, 1000, col=adjustcolor("coral", 0.2), border='white',lwd=2)
polygon(c(s1avg$years,rev(s1avg$years)),c((s1avg$upper),rev((s1avg$lower))),col = alpha("grey75",.75), border = FALSE)
lines(s1avg$years,s1avg$mean,lwd=4,col='steelblue')
text(1750,0.75,substitute(paste(bold('Pre E-A Settlement'))))
text(1930,0.75,substitute(paste(bold('Post E-A Settlement'))))
plot(s1.seg.2,add=TRUE,link=FALSE,lwd=3, col='coral',lty=2262)
fig_label("B",region="plot", cex=2) 

dev.off()

####Figure 3#####

yrs<-match(seq(1860,1940,by=20),s1avg$years)
brate1<-rep(s1avg$mean[which(s1avg$years==1860)],length(yrs))+1 ##Under 100%  survival the population growth rate is establshishment +1
brateall<-s1avg$mean[yrs]+1
Btotal<-diff(c(0.15,0.15*cumprod(brateall))) ##simulations of the populations under observed average (changing) growth rate. Estsblishment is the change in population size when survival is 100%
Bexpg<-diff(c(0.15,0.15*cumprod(brate1))) ##simulations of the populations under constant 1860 growth rate.


Bbrate<-Btotal-Bexpg ##amount of establishment due to changing estab rates is Total minus constant. This assumes any compounding population growth from changing establishment rates get propgated on towards the changing establishment rate estimates. 
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
  geom_bar(position="stack", stat="identity")+scale_fill_manual(values=c(alpha('coral',0.7),alpha('steelblue',0.7)))+scale_x_continuous(breaks=c(1860,1880, 1900, 1920, 1940))+ylab("Establishment")+xlab("Year (CE)")+theme_bw()+geom_text(x = 1855, y = .19
                                                                                                                                                                            , width = unit(.8, "inch")
                                                                                                                                                                            , label = "A"
                                                                                                                                                                           
                                                                                                                                                                            , color = "black"
                                                                                                                                                                            , size=6)

# Plot
p2<-ggplot(Pdata, aes(x=Year, y=Percent, fill=Mechanism)) + 
  geom_area(alpha=0.6 , size=1, colour="white")+scale_fill_manual(values=c('coral','steelblue'))+theme_bw()+theme(legend.position="none")+ylab("Proportional contribution")+xlab("Year (CE)")+geom_text(x = 1865, y = 1
                                                                                                                                                                                            , width = unit(.8, "inch")
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

ggsave('Fig3.tiff', plot=F3, device="tiff")

#####Figure 4
tiff("Fig4.tiff",width = 8,height=6, units="in", res=300)    

par(mfrow=c(2,1),mar = c(0.5,6,4.5,.5))
plot(0,0,ylim=c(0,1),xlim=c(1600,1990),xaxt="n", xlab='', ylab=expression("Survival rate (20yrs"^-1~")"))
rect(0, -1000, 1850, 1000, col=adjustcolor("grey",0.3 ), border='white',lwd=2)
rect(1850, -1000, 2050, 1000, col=adjustcolor("coral", 0.2), border='white',lwd=2)


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



par(mar = c(4.5,6,.5,.5))
plot(0,0,ylim=c(0,.7),xlim=c(1600,1990), ylab=expression("Establishment rate (20yrs"^-1~")"), xlab='Year (CE)')
rect(0, -1000, 1850, 1000, col=adjustcolor("grey",0.3 ), border='white',lwd=2)
rect(1850, -1000, 2050, 1000, col=adjustcolor("coral", 0.2), border='white',lwd=2)
text(1730,0.65,substitute(paste(bold('Pre E-A Settlement'))))
text(1940,0.65,substitute(paste(bold('Post E-A Settlement'))))


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

dev.off()
####################################################

#Figure 5###

temps<-read.csv("TempRecords.csv",header=T) ##import reconstructed temp data
bins20<-seq(1600,2000,by=20)
Temp<-numeric(length(bins20)-1)
for(i in 2:length(bins20)){
  Temp[i-1]<-mean(temps$MaxTemp[which(temps$Year>=bins20[i-1] & temps$Year<bins20[i])],na.rm=T) ####Average temp by 20 year bins
}



tiff("Fig5.tiff",width = 9,height=2.9,units="in", res=300)    

layout(mat=matrix(c(1,1,1,1,2,2),1,6))
par(mar = c(4.5,4.6,4,4.5))

plot(seq(1600,1980,by=20),s1avg$mean[18:37], type='l', ylim=c(0,0.8),lwd=3,col='steelblue', ylab=expression("Establishment rate (20yrs"^-1~")"), xlab='Year (CE)')
polygon(c(s1avg$years[18:37],rev(s1avg$years[18:37])),c((s1avg$upper[18:37]),rev((s1avg$lower[18:37]))),col = alpha("grey75",.75), border = FALSE)
lines(seq(1600,1980,by=20),s1avg$mean[18:37], type='l', ylim=c(0,0.8),lwd=3,col='steelblue')

par(new=T)

plot(seq(1600,1980,by=20),Temp, type='l',col='coral',lwd=3,axes = FALSE, bty = "n", xlab = "", ylab = "",ylim=c(14.5,16.5))
axis(side=4, col.axis="coral" ,col="coral" )
mtext(expression("Mean max temperature (\u00B0C 20yrs"^-1~")"),side=4,col="coral",line=2.5,cex=0.7) 
fig_label("A",region="plot", cex=2) 

cor(Temp,s1avg$mean[18:37])
Tmod<-lm(s1avg$mean[18:37]~Temp)
summary(Tmod)
newdata=as.data.frame(seq(14.5,16.5,length=20))
colnames(newdata)<-'Temp'
prednew<-predict.lm(Tmod, newdata=newdata,se.fit=T)

plot(Temp,s1avg$mean[18:37], pch=19,ylab=expression("Establishment rate (20yrs"^-1~")"),xlab=expression("Mean max temperature (\u00B0C 20yrs"^-1~")"))
lines(newdata$Temp,prednew$fit,lwd=3)
polygon(c(newdata$Temp,rev(newdata$Temp)),c((prednew$fit+prednew$se.fit),rev((prednew$fit-prednew$se.fit))),col = alpha("grey75",.75), border = FALSE)

text(15.5,0.45,labels=paste("Correlation=",round(cor(Temp,s1avg$mean[18:37]),2)))
text(15.5,0.47,labels=paste("Slope=",round(Tmod$coefficients[2],2)))
fig_label("B",region="plot", cex=2) 

dev.off() 

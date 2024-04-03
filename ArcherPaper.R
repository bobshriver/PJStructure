ArcherData<-read.csv('Archer1989.csv')


ArcherData$N<-cumsum(ArcherData$B)
years<-rev(seq(1985,1770,by=-20))
N<-ArcherData$N[match(years,ArcherData$Year)]/max(ArcherData$N)

Brate<-c(NA,exp(diff(log(N))))-1


par(mar = c(4.5,4.6,4,4.5))

plot(years-20,Brate, pch=19, lwd=3,ylab=expression("Establishment Rate (20yrs"^-1~")"), xlab="Year (CE)")
par(new=T)
plot(years-20,N, type='l',col='coral',lwd=3,axes = FALSE, bty = "n", xlab = "", ylab = "")
axis(side=4, col.axis="coral" ,col="coral" )
mtext(expression("Cumulative Establishment"),side=4,col="coral",line=2.5,cex=1) 

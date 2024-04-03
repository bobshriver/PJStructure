


#######################

par(mfrow=c(2,2))

plot((1*cumprod(c(rep(1.2,10),rep(1.6,2)))), type='l',lwd=3, col='steelblue')
lines((1*cumprod(c(rep(1.2,10),rep(1.2,3)))), col='coral',lwd=3)
lines((1*cumprod(c(rep(1.2,10),1.6,rep(1.2,2)))),lwd=3,col='steelblue')


plot((1*cumprod(c(rep(1.2,10),1.3,rep(1.2,2)))), type='l',lwd=3,col='steelblue')
lines((1*cumprod(c(rep(1.2,10),rep(1.2,3)))), col='coral',lwd=3)
lines((1*cumprod(c(rep(1.2,10),1.3,rep(1.2,2)))),lwd=3,col='steelblue')



plot(((0.1*1.4^(1:20))),lwd=5, col='steelblue',axes = FALSE, bty = "n", xlab = "", ylab = "")
box()
arrows(12,38,14.5,60,lwd=3, length=0.2)
arrows(12,35,14.5,18,lwd=3, length=0.2)
text(8,36,lab='Change point?')
par(new=T)
plot(rep(1.4,25),lwd=4, col='coral',type='l',axes = FALSE, bty = "n", xlab = "", ylab = "")
mtext("Per-capita establishment rate",side=4,col="Coral",line=1) 
mtext("Establishment",side=2,col="steelblue",line=1) 
mtext("Time",side=1,line=1) 



#lines(diff(1*cumprod(c(rep(1.2,10),rep(1.2,3)))), col40='coral',lwd=3)
lines(diff(.1*cumprod(c(rep(s1avg$mean[which(s1avg$years==1840)],length(yrs))+1))),lwd=3,col='coral',type='s')


plot(diff(1*cumprod(c(rep(1.2,15)))), type='l',lwd=3,col='steelblue')
lines(diff(1*cumprod(c(rep(1.2,10),rep(1.2,3)))), col='coral',lwd=3)
lines(diff(1*cumprod(c(rep(1.2,10),1.3,rep(1.2,2)))),lwd=3,col='steelblue')




#Extra

E<-s1avg$mean[18:37]
Years<-seq(1600,1980,by=20)

dat <- data.frame(Years, E, Temp, Precip)
ggplot(dat, aes(Years, E)) +
  geom_line(color="steelblue",size=2) +  geom_ribbon(aes(ymin=s1avg$lower[18:37], ymax=s1avg$upper[18:37]), alpha=0.2)+ylim(0,0.8)+
  geom_line(
    aes(y = Tmod), color='coral',size=2,
    data = ~ transform(., Tmod = scales::rescale(Temp, range(c(0,s1avg$upper[18:37],s1avg$lower[18:37])), range(Temp)))
  ) +
  scale_y_continuous(
    sec.axis = sec_axis(~ scales::rescale(., range(dat$Temp), range(c(0,s1avg$upper[18:37],s1avg$lower[18:37]))),
                        breaks = c(14.5,15,15.5,16,16.5))
  )




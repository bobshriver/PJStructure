

svariable<-function(x,sstart,backward=T,agint=20){ 
  
  TimeSh<-seq(min(x$Talign),2020,by=20)
  Time<-rep(TimeSh,max(x$Dataset))
  
  srateall<-rep(NA,length(TimeSh))
  s<-logit(sstart)
  for(i in 1:length(TimeSh)){s<-(s+rnorm(1,0,0.45)) ##Simulate random walk of survival values
  srateall[i]<-inv.logit(s)}
  if(backward==T){svec<-rev(srateall)}# forward or backward switched the direction of the random walk either present to past, or past to present
  if(backward==F){svec<-(srateall)}
  
  ###Create containers for output

  dataout<-as.data.frame(Time)
  dataout$dataset<-(rep(1:max(x$Dataset), each=length(TimeSh)))
  dataout$grate<-NA
  dataout$brate<-NA
  dataout$srate<-NA
  dataout$B<-NA
  dataout$Nstar<-NA
  dataout$Bstar<-NA
  dataout$Relpop<-NA
  dataout$Relt_1<-NA
  dataout$RelB<-NA
  dataout$Species<-NA
  dataout$deltabrate<-NA
  dataout$deltaB<-NA
  dataout$deltaN<-NA
  
  ####Loop over all datasets####
  for (i in 1:max(x$Dataset)){
    agedata<-x[which(x$Dataset==i),]#extract individual dataset
    agedata$firstsurvival<-min(agedata$CensusOffset)/agedata$Interval[1] ###Correction needed for offset between decade collected and last decade of observed establishment 
    i10<-(agedata$firstsurvival[1]*agedata$Interval[1])
    
    if (i10<=agint) {surv1<- svec[which(TimeSh==max(agedata$Talign)+agedata$Interval[1])]^(i10/agint)}
    
    if (i10>agint) { wsteps<-floor(i10/agint) ###Whole time steps missing in periods
                psteps<-(i10/agint)-wsteps ###partial time steps missing
                steps<-ceiling(i10/agint) ##Total number of steps, rounded up
       surv1<-prod(svec[(which(TimeSh==max(agedata$Talign)+agedata$Interval[1])):(which(TimeSh==max(agedata$Talign)+agedata$Interval[1])+wsteps-1)])*(svec[(which(TimeSh==max(agedata$Talign)+agedata$Interval[1])+steps-1)]^(psteps)) ##survival from the last establishment cohort to census
      }
    
    
    srateint<-numeric(length(agedata$B)+1)# Create vector of survival values for each interval. It's 1 longer than the number of intervals to include survival from the last cohort to census
    srateint[length(agedata$B)+1]<-surv1 ###last slot is survival from the last cohort to census
    for(ii in (length(agedata$B)):1){
      srateint[ii]<-svec[which(TimeSh>=agedata$Talign[ii] & TimeSh<agedata$Talign[ii]+agedata$Interval[1])] ###loop over and fill in the rest. Since survival rates are 5 year they need to be aggregated across the interval with the prodcuct
      
    }
    cumsurv<-rev(cumprod(rev(srateint))) #has to be double reversed for cumprod to do the accumulation in the right direction, seems like the function could use a reverse option
    agedata$Bstar<-agedata$B/cumsurv[-1] ###Correcting B based on cumulative survival of cohort from establishment to observation. Drop fist value which would be cumulative survival of cohort before the first est year
    
    N<-0 ###starting population size, loop over intervals
    for (t in 1:length(agedata$B)) {
      N<-agedata$Bstar[t]+N*srateint[t] # population size calculation
      agedata$Nstar[t]<-N
    }
    agedata$grate<-c(NA,diff(log(agedata$Nstar))/(agedata$Interval[-1]/agint)) #20 year growth rate on log scale corrected for interval length 
    agedata$brate<-exp(agedata$grate)-(srateint[-length(srateint)]^(1/(agedata$Interval[1]/agint)))#Establishment rate calculation. Remove last element in survival which is survival from the last cohort to the time of census
    
    ##Helpful for plugging results into containers
    yrs<-dataout$Time%in%agedata$Talign 
    loc<-which(dataout$dataset==i & yrs)
    
    ###Plug results into container###
    dataout$grate[loc]<-exp(agedata$grate)
    dataout$brate[loc]<-agedata$brate
    dataout$srate[loc]<-srateint[-length(srateint)]^(1/(agedata$Interval[1]/agint))
    dataout$B[loc]<-agedata$B
    dataout$Nstar[loc]<-agedata$Nstar
    dataout$Bstar[loc]<-agedata$Bstar
    dataout$Relpop[loc]<-agedata$Nstar/max(agedata$Nstar)
    dataout$Relt_1[loc]<-c(NA,agedata$Nstar[-length(agedata$B)])/max(agedata$Nstar)
    dataout$RelB[loc]<-agedata$Bstar/sum(agedata$Bstar)
    dataout$deltabrate[loc]<-c(NA,diff(agedata$brate))
    dataout$deltaB[loc]<-c(NA,diff(agedata$Bstar/sum(agedata$Bstar)))
    dataout$deltaN[loc]<-c(NA,diff(agedata$Nstar/max(agedata$Nstar)))
    dataout$Species[loc]<-agedata$Species
    
    
  }
  
  return(dataout)
  
  
  
}

blag<-function(x,lag, agint=20){
  
  ###Create containers for output
  TimeSh<-seq(min(x$Talign),max(x$Talign),by=agint)
  Time<-rep(TimeSh,max(x$Dataset))
  dataout<-as.data.frame(Time)
  dataout$dataset<-(rep(1:max(x$Dataset), each=length(TimeSh)))
  dataout$grate<-NA
  dataout$brate<-NA
  dataout$srate<-NA
  dataout$B<-NA
  dataout$Nstar<-NA
  dataout$Bstar<-NA
  dataout$Relpop<-NA
  dataout$Relt_1<-NA
  dataout$RelB<-NA
  dataout$Species<-NA
  dataout$deltabrate<-NA
  dataout$deltaB<-NA
  dataout$deltaN<-NA
  
  ####Loop over all datasets####
  for (i in 1:max(x$Dataset)){
    agedata<-x[which(x$Dataset==i),]#extract individual dataset

    agedata$Nstar<-cumsum(agedata$B)
      
    
        laguse<-lag/agedata$Interval[1]
        agedata$brate<-c(rep(NA,laguse),((agedata$B[-(1:laguse)])/(agedata$Nstar[-((length(agedata$Dataset)-(laguse-1)):length(agedata$Dataset))])))
      
    
    ##Helpful for plugging results into containers
    yrs<-dataout$Time%in%agedata$Talign 
    loc<-which(dataout$dataset==i & yrs)
    
    ###Plug results into container###
    dataout$grate[loc]<-NA
    dataout$brate[loc]<-agedata$brate
    dataout$srate[loc]<-1
    dataout$B[loc]<-agedata$B
    dataout$Nstar[loc]<-agedata$Nstar
    dataout$Bstar[loc]<-NA
    dataout$Relpop[loc]<-agedata$Nstar/max(agedata$Nstar)
    dataout$Relt_1[loc]<-c(NA,agedata$Nstar[-length(agedata$B)])/max(agedata$Nstar)

    dataout$Species[loc]<-agedata$Species
    
    
  }
  
  return(dataout)
  
  
  
}

decadeavg<-function(x,var='brate'){
  years<-unique(x$Time)
  out<-as.data.frame(years)
 
    out$mean<-by(x[,which(colnames(x)==var)],x$Time,mean,na.rm=T)
    out$upper<-out$mean+by(x[,which(colnames(x)==var)],x$Time,std.error,na.rm=T)*1.96
    out$lower<-out$mean-by(x[,which(colnames(x)==var)],x$Time,std.error,na.rm=T)*1.96
    
  
  return(out)
  
}

mvavg<-function(x,var='brate', smoothyears){
  years<-unique(x$Time)
  out<-as.data.frame(years)
  out$mean<-NA
  out$upper<-NA
  out$lower<-NA
  for(i in (smoothyears/5):length(years)){
    use<-which(x$Time<=years[i] & x$Time>years[i]-smoothyears)
    out$mean[i]<-mean(x[use,which(colnames(x)==var)],na.rm=T)
    out$upper[i]<-out$mean[i]+std.error(x[use,which(colnames(x)==var)],na.rm = T)*1.96
    out$lower[i]<-out$mean[i]-std.error(x[use,which(colnames(x)==var)],na.rm = T)*1.96
    
  }
  
  return(out)
  
}




sconstant<-function(x,srate,agint=20){
  
  ###Create containers for output
  TimeSh<-seq(min(x$Talign),max(x$Talign),by=agint)
  Time<-rep(TimeSh,max(x$Dataset))
  dataout<-as.data.frame(Time)
  dataout$dataset<-(rep(1:max(x$Dataset), each=length(TimeSh)))
  dataout$grate<-NA
  dataout$brate<-NA
  dataout$srate<-NA
  dataout$B<-NA
  dataout$Nstar<-NA
  dataout$Bstar<-NA
  dataout$Relpop<-NA
  dataout$Relt_1<-NA
  dataout$RelB<-NA
  dataout$Species<-NA
  dataout$deltabrate<-NA
  dataout$deltaB<-NA
  dataout$deltaN<-NA
  dataout$bratermean<-NA
  ####Loop over all datasets####
  for (i in 1:max(x$Dataset)){
    agedata<-x[which(x$Dataset==i),]#extract individual dataset
    agedata$firstsurvival<-min(agedata$CensusOffset)/agedata$Interval[1] ###Correction needed for offset between decade collected and last decade of observed establishment 
    srateint<-srate^(agedata$Interval[1]/agint) ###Adjust survival rate for interval length. 
    sratestart<-srateint^agedata$firstsurvival[1] ###calculating survival from last establishment cohort to when data was collected 
    agedata$Bstar<-agedata$B/c(sratestart*(srateint^((length(agedata$B)-1):1)),sratestart) ###Correcting B based on cumulative survival of cohort from establishment to observation. Survival is cumulative product of sratestart and each subsequent interval going back to the establishment cohort bin
    
    N<-0 ###starting population size, loop over intervals
    for (t in 1:length(agedata$B)) {
      N<-agedata$Bstar[t]+N*srateint # population size calculation
      agedata$Nstar[t]<-N
    }
    agedata$grate<-c(NA,diff(log(agedata$Nstar))/(agedata$Interval[-1]/agint)) #20 year growth rate on log scale corrected for interval length 
    agedata$brate<-exp(agedata$grate)-srate #Establishment rate calculation. Since grate is now on a 20 year time period for all, use 20 year survival rate not srateint
    
    ##Helpful for plugging results into containers
    yrs<-dataout$Time%in%agedata$Talign
    loc<-which(dataout$dataset==i & yrs)
    
    ###Plug results into container###
    dataout$grate[loc]<-agedata$grate
    dataout$brate[loc]<-agedata$brate
    dataout$srate[loc]<-srate
    dataout$B[loc]<-agedata$B
    dataout$Nstar[loc]<-agedata$Nstar
    dataout$Bstar[loc]<-agedata$Bstar
    dataout$Relpop[loc]<-agedata$Nstar/max(agedata$Nstar)
    dataout$Relt_1[loc]<-c(NA,agedata$Nstar[-length(agedata$B)])/max(agedata$Nstar)
    dataout$RelB[loc]<-agedata$Bstar/sum(agedata$Bstar)
    dataout$deltabrate[loc]<-c(NA,diff(agedata$brate))
    dataout$deltaB[loc]<-c(NA,diff(agedata$Bstar/sum(agedata$Bstar)))
    dataout$deltaN[loc]<-c(NA,diff(agedata$Nstar/max(agedata$Nstar)))
    dataout$Species[loc]<-agedata$Species
    dataout$bratermean[loc]<-agedata$brate-mean(agedata$brate[which(agedata$Talign>=1600 & agedata$Talign<1850)],na.rm=T)
    
  }
  
  return(dataout)
  
  
  
}


fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) { ### Code borrowed from https://gist.github.com/Pakillo/5712db8a54f3efb3b1583c52a2e2e270
  
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
  
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}

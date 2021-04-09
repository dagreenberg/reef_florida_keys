rm(list=ls())
setwd("C:/Users/14388/Desktop/reef_florida_keys_data")
library(rstan);library(loo);library(bayesplot);library(dplyr);library(tidyverse);library(stringr);library(lubridate);library(rlist)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

####Functions####
`%notin%`<- Negate(`%in%`)

abund_tranfs<- function(probs){
  sum<- probs[1]*0+probs[2]*1+probs[3]*2+probs[4]*11+probs[5]*101
  return(sum)
}

rvc_filter = function(x,GZ,sp){
  x$SSU_YEAR<- paste(x$PRIMARY_SAMPLE_UNIT,x$STATION_NR,x$YEAR,sep='_')
  x1= x %>% subset(region.id==GZ) %>% select(SSU_YEAR,SPECIES_CD,everything())
  x2= complete(x1,SSU_YEAR,nesting(SPECIES_CD),fill=list(NUM=0))
  zeros = anti_join(x2,x1)
  zeros[,3:29]<- x1[match(zeros$SSU_YEAR,x1$SSU_YEAR),3:29]
  x3= rbind(x1,zeros)
  x3$HAB_CD2<- gsub('\\_.*','',x3$HABITAT_CD)
  year_index<- data.frame(yr=seq(1993,2018),y.ind=seq(1,26))
  x3$year_index=x3$y.ind[match(x3$YEAR,year_index$yr)]
  
  rvc_occs<- list()
  for(i in 1:nrow(sp)){
    x4= subset(x3,SPECIES_CD==sp$rvc_code[i])
    x5=  x4 %>% dplyr::group_by(SSU_YEAR) %>%
      dplyr::summarise(NUM.total=sum(NUM),occ=NA) %>% #Sums up the number of counts
      mutate(occ=ifelse(NUM.total>0,1,0)) %>% arrange(SSU_YEAR) #Also scores presence/absence at the SSU level
    x5[,4:32]<- x4[match(x5$SSU_YEAR,x4$SSU_YEAR),2:30]
    x5<- transform(x5,psu_id=match(LAT_LON,unique(LAT_LON)))
    x5$NUM.total2<- ceiling(x5$NUM.total)
    rvc_occs[[i]]=x5
  }
  return(rvc_occs)
}

ts_rvc = function(x){ #Takes the output from the previous function
  ts<- list()
  
  for(i in 1:length(x)){
    x1 = x[[i]]
    x2= x1 %>% group_by(YEAR) %>% summarize(n.occ=sum(occ),n.surv=n(),p.occ=n.occ/n.surv,sp=unique(SPECIES_CD))
    x3= x1 %>% group_by(YEAR,PRIMARY_SAMPLE_UNIT) %>% summarize(psu_abund=mean(NUM.total2)) %>% group_by(YEAR) %>% summarize(mean_abund=mean(psu_abund),sd_abund=sd(psu_abund))
    x4=left_join(x2,x3) %>% complete(YEAR=seq(1993,2018))
    
    for(z in 1:nrow(x4)){
      if(is.na(x4$p.occ[z])==T){next}
      if(x4$p.occ[z]==0){
        x4$p.occ[z]=NA
        x4$mean_abund[z]=NA
      }
    }
    
    ts[[i]]=x4
  }
  return(ts)
}


reef_filter = function(R,GZ,sp,geog){ #function to trim dataframe and add in occurrence data by species
  occ_list<- list()
  TempDat<-R %>% subset(site4==GZ) %>% subset(year<=2018) %>% select('formid','speciesid','abundance',everything())
  TempDat$hab_class<- geog$hab_class[match(TempDat$geogr,geog$geogid)]
  Zeros<- complete(TempDat,formid,nesting(speciesid),
                   fill=list(abundance=0)) %>%
    anti_join(TempDat) #Creates a dataset where each species is represented in every survey (adding in zero observations)
  m<- match(Zeros$formid,TempDat$formid) #Match the zero data frame to the rest of the data
  Zeros[,4:ncol(Zeros)]<- TempDat[m,4:ncol(TempDat)] #Replicate the survey-level data (except abundance)
  TempDat2<- rbind(TempDat,Zeros) %>% arrange(formid,speciesid) #Combine presence and absence data for all species
  TempDat2$occ=ifelse(TempDat2$abundance>0,1,0) #Code presence/absence based on abundance
  TempDat2<- TempDat2 %>% mutate(exp_binary=ifelse(exp=='E',1,0))
  TempDat2$abundance2<- TempDat2$abundance+1 #For modelling abundance categories, has to start at 1
  TempDat2$abund_trans <-NA #Transform abundance categories into the minimum counts
  TempDat2<- TempDat2%>%
    mutate(
      abund_trans=ifelse(abundance==0,0,abund_trans),
      abund_trans=ifelse(abundance==1,1,abund_trans),
      abund_trans=ifelse(abundance==2,2,abund_trans),
      abund_trans=ifelse(abundance==3,11,abund_trans),
      abund_trans=ifelse(abundance==4,101,abund_trans)
    )
  surveyors<- TempDat2 %>% group_by(fish_memberid) %>% summarize(n=length(unique(formid))) #Determine survey effort of each diver in the region
  surveyors_trim<- subset(surveyors,n>=5) #only keep surveys by members with 5 or more dives
  
  TempDat3<- subset(TempDat2, fish_memberid %in% surveyors_trim$fish_memberid) #Trim out surveys from surveys with less than 5
  site_subsets<-TempDat3 %>% group_by(geogr) %>% summarize(n=n_distinct(formid),hab_class=unique(hab_class)) %>% subset(n>=5 & is.na(hab_class)==F) #Calculate surveys per site
  
  TempDat4<- TempDat3 %>% subset(fish_memberid %in% surveyors_trim$fish_memberid) %>% subset(geogr %in% site_subsets$geogr) #Only keep well surveyed sites (n>=5)
  TempDat4$hab_class2<- gsub('\\_.*','',TempDat4$hab_class)
  survs<- TempDat4 %>% distinct(formid,.keep_all = T) %>% group_by(geogr,day,month,year) %>% summarize(n=n()) %>% mutate(site_dmy=paste(geogr,day,month,year,sep='')) %>% subset(n>1)
  TempDat4<- TempDat4 %>% mutate(site_dmy=ifelse(paste(geogr,day,month,year,sep='') %in% survs$site_dmy,paste(geogr,day,month,year,sep=''),'baseline'))
  
  for(i in 1:nrow(sp)){
    occ_list[[i]]<- subset(TempDat4,speciesid==sp$speciesid[i]) #Subset out each species in the provided dataframe
  }
  return(occ_list)
}

ts_reef = function(X,sp){
  ts_list<- list()
  for(i in 1:length(X)){
    sp_x<-X[[i]] #Subset out each species in the provided dataframe
    occ_by_year<- sp_x %>% group_by(year) %>% summarize(n.occ=sum(occ),n.surv=n(),p.occ=n.occ/n.surv,sp=sp$commonname[match(unique(sp_x$speciesid),sp$speciesid)])
    abun_by_year<- sp_x %>% group_by(year,geogr) %>% summarize(site_abun=mean(abund_trans),sd_abund=sd(abund_trans),n.surv=n()) %>% group_by(year) %>% summarize(mean_abund=mean(site_abun),n.survs=sum(n.surv),n.sites=n())
    total_sd<- sp_x %>% group_by(year) %>% summarize(sd=sd(abund_trans))
    
    comb<- left_join(occ_by_year,abun_by_year)
    comb2<- left_join(comb,total_sd)
    ts_list[[i]]<- comb2
  }
  return(ts_list)
}


TS_occ_plot_MARSS<- function(ts1,ts2,sp,GZ,mod){
  pdf(paste(paste(sp,GZ,sep='_'),'.pdf',sep=''),width=8,height=6)
  par(xpd=T)
  plot(plogis(ts1[1,])~c(seq(1993,2018)),type='n',ylim=c(min(na.omit(c(ts1,plogis(ts2)))),max(na.omit(c(ts1,plogis(ts2))))),col='darkblue',bty='l',ylab=expression('Sighting Probability'),xlab='Year',main=paste(sp,GZ,sep=', '))
  
  if(mod==1){
    lines(plogis(fit_1$states[1,])~c(seq(1993,2018)),lty=5,lwd=2,col='slategray')
    lines(plogis(fit_1$states[1,]+fit_1$par$A[,1])~c(seq(1993,2018)),lty=5,lwd=2,col='slategray')
    x<- c(c(seq(1993,2018)), rev(c(seq(1993,2018))))
    y1<- c(plogis(fit_1$states[1,]-fit_1$states.se[1,]*2), rev(plogis(fit_1$states[1,]+fit_1$states.se[1,]*2)))
    y2<- c(plogis(fit_1$states[1,]+fit_1$par$A[,1]-fit_1$states.se[1,]*2), rev(plogis(fit_1$states[1,]+fit_1$par$A[,1]+fit_1$states.se[1,]*2)))
    polygon(x, y1, col = adjustcolor('slategray', alpha = 0.1), border=NA) # Add uncertainty polygon
    polygon(x, y2, col = adjustcolor('slategray', alpha = 0.1), border=NA) # Add uncertainty polygon
  }
  if(mod==2){
    lines(plogis(fit_3$states[1,])~c(seq(1993,2018)),lty=5,lwd=2,col='darkcyan')
    lines(plogis(fit_3$states[2,])~c(seq(1993,2018)),lty=5,lwd=2,col='darksalmon')
    x<- c(c(seq(1993,2018)), rev(c(seq(1993,2018))))
    y1<- c(plogis(fit_3$states[1,]-fit_3$states.se[1,]*2), rev(plogis(fit_3$states[1,]+fit_3$states.se[1,]*2)))
    y2<- c(plogis(fit_3$states[2,]-fit_3$states.se[2,]*2), rev(plogis(fit_3$states[2,]+fit_3$states.se[2,]*2)))
    polygon(x, y1, col = adjustcolor('darkcyan', alpha = 0.1), border=NA) # Add uncertainty polygon
    polygon(x, y2, col = adjustcolor('darksalmon', alpha = 0.2), border=NA) # Add uncertainty polygon
  }
  lines(ts1[1,]~c(seq(1993,2018)),col=adjustcolor('navy',alpha.f=0.5),lwd=2)
  points(ts1[1,]~c(seq(1993,2018)),col='white',pch=21,bg=adjustcolor('navy',alpha.f=0.5),cex=1.5)
  lines(plogis(ts2[1,])~c(seq(1993,2018)),col='dodgerblue4',lwd=2)
  points(plogis(ts2[1,])~c(seq(1993,2018)),col='white',pch=21,bg='dodgerblue4',cex=1.5)
  
  lines(ts1[2,]~c(seq(1993,2018)),col=adjustcolor('darkred',alpha.f=0.5),lwd=2)
  points(ts1[2,]~c(seq(1993,2018)),col='white',pch=21,bg=adjustcolor('darkred',alpha.f=0.5),cex=1.5)
  lines(plogis(ts2[2,])~c(seq(1993,2018)),col='firebrick4',lwd=2)
  points(plogis(ts2[2,])~c(seq(1993,2018)),col='white',pch=21,bg='firebrick4',cex=1.5)
  
  legend(2013,c(max(na.omit(c(ts1,plogis(ts2))))*1.05),c('Obs. NOAA surveys','Obs. REEF surveys','Est. NOAA surveys','Est. REEF surveys'),text.col=c(adjustcolor('navy',alpha.f=0.5),adjustcolor('darkred',alpha.f=0.5),'dodgerblue4','firebrick4'),bty='n')
  dev.off(paste(paste(sp,GZ,sep='_'),'.pdf',sep=''))
}

####Stan occupancy plot function####
TS_stan_abund_plot_MARSS<- function(ts1,ts2,sp,GZ,params){
 # pdf(paste(paste(i,sp,GZ,sep='_'),'.pdf',sep=''),width=8,height=6)
  par(xpd=T)
  plot(ts1$mean_abund~c(seq(1993,2018)),type='n',ylim=c(min(na.omit(c(ts1$mean_abund,ts2$mean_abund))),max(na.omit(c(ts1$mean_abund,ts2$mean_abund)))),col='darkblue',bty='l',ylab=expression('Mean count per survey'),xlab='Year',main=paste(sp,GZ,sep=', '))
  
  lambda_mat<- list()
    
  reef_coef<- data.frame(p_0.y=NA,p_1.y=NA,p_2.y=NA,p_11.y=NA,p_101.y=NA,p_0.x=NA,p_1.x=NA,p_2.x=NA,p_11.x=NA,p_101.x=NA,lambda.y=NA,lambda.x=NA,iter=seq(1,1200))
  for(i in 1:26){
    if(ncol(params$c)==2){
      reef_coef[,1]<- plogis(params$c[,1]-params$a_yr2[,i])
      reef_coef[,2]<-plogis(params$c[,2]-params$a_yr2[,i])-plogis(params$c[,1]-params$a_yr2[,i])
      reef_coef[,3]<-1-plogis(params$c[,2]-params$a_yr2[,i])
      reef_coef[,4]<- 0
      reef_coef[,5]<- 0
      reef_coef[,6]=plogis(params$c[,1]-params$x2[,i])
      reef_coef[,7]<-plogis(params$c[,2]-params$x2[,i])-plogis(params$c[,1]-params$x2[,i])
      reef_coef[,8]<-1-plogis(params$c[,2]-params$x2[,i])
      reef_coef[,9]<- 0
      reef_coef[,10]<- 0
      reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
      reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
          }
    if(ncol(params$c)==3){
      reef_coef[,1]<- plogis(params$c[,1]-params$a_yr2[,i])
      reef_coef[,2]<-plogis(params$c[,2]-params$a_yr2[,i])-plogis(params$c[,1]-params$a_yr2[,i])
      reef_coef[,3]<-plogis(params$c[,3]-params$a_yr2[,i])-plogis(params$c[,2]-params$a_yr2[,i])
      reef_coef[,4]<- 1-plogis(params$c[,3]-params$a_yr2[,i])
      reef_coef[,5]<- 0
      reef_coef[,6]=plogis(params$c[,1]-params$x2[,i])
      reef_coef[,7]<-plogis(params$c[,2]-params$x2[,i])-plogis(params$c[,1]-params$x2[,i])
      reef_coef[,8]<-plogis(params$c[,3]-params$x2[,i])-plogis(params$c[,2]-params$x2[,i])
      reef_coef[,9]<- 1-plogis(params$c[,3]-params$x2[,i])
      reef_coef[,10]<- 0
      reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
      reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
    }
    if(ncol(params$c)==4){
      reef_coef[,1]<- plogis(params$c[,1]-params$a_yr2[,i])
      reef_coef[,2]<-plogis(params$c[,2]-params$a_yr2[,i])-plogis(params$c[,1]-params$a_yr2[,i])
      reef_coef[,3]<-plogis(params$c[,3]-params$a_yr2[,i])-plogis(params$c[,2]-params$a_yr2[,i])
      reef_coef[,4]<- plogis(params$c[,4]-params$a_yr2[,i])-plogis(params$c[,3]-params$a_yr2[,i])
      reef_coef[,5]<- 1- plogis(params$c[,4]-params$a_yr2[,i])
      reef_coef[,6]=plogis(params$c[,1]-params$x2[,i])
      reef_coef[,7]<-plogis(params$c[,2]-params$x2[,i])-plogis(params$c[,1]-params$x2[,i])
      reef_coef[,8]<-plogis(params$c[,3]-params$x2[,i])-plogis(params$c[,2]-params$x2[,i])
      reef_coef[,9]<- plogis(params$c[,4]-params$x2[,i])-plogis(params$c[,3]-params$x2[,i])
      reef_coef[,10]<- 1-plogis(params$c[,4]-params$x2[,i])
      reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
      reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
    }
    
    lambda_mat[[i]]<- reef_coef
    
  }
  
    x_mat<- data.frame(median.rvc=NA,l.95=NA,u.95=NA,median.reef=NA,l.95.rf=NA,u.95.rf=NA)
    for(i in 1:26){
      x_mat[i,1]=exp(median(params$x1[,i]))
      x_mat[i,2]=exp(quantile(params$x1[,i],0.1))
      x_mat[i,3]=exp(quantile(params$x1[,i],0.9))
      x_mat[i,4]=median(lambda_mat[[i]]$lambda.x)
      x_mat[i,5]=quantile(lambda_mat[[i]]$lambda.x,0.1)
      x_mat[i,6]=quantile(lambda_mat[[i]]$lambda.x,0.9)
    }
    
    y_mat_rvc<- data.frame(year=ts1$YEAR[complete.cases(ts1)],median.rvc=NA,l.95.rvc=NA,u.95.rvc=NA)
    y_mat_reef<- data.frame(year=seq(1993,2018),median.reef=NA,l.95.rf=NA,u.95.rf=NA)
    for(i in 1:23){
      y_mat_rvc[i,2]=exp(median(params$a_yr1[,i]))
      y_mat_rvc[i,3]=exp(quantile(params$a_yr1[,i],0.1))
      y_mat_rvc[i,4]=exp(quantile(params$a_yr1[,i],0.9))
    }
    
    for(i in 1:26){
      y_mat_reef[i,2]=median(lambda_mat[[i]]$lambda.y)
      y_mat_reef[i,3]=quantile(lambda_mat[[i]]$lambda.y,0.1)
      y_mat_reef[i,4]=quantile(lambda_mat[[i]]$lambda.y,0.9)
    }
    y_mat<- full_join(y_mat_reef,y_mat_rvc)
    
    lines(x_mat[,1]~c(seq(1993,2018)),lty=5,lwd=2,col='darkcyan')
    lines(x_mat[,4]~c(seq(1993,2018)),lty=5,lwd=2,col='darksalmon')
    x<- c(c(seq(1993,2018)), rev(c(seq(1993,2018))))
    y1<- c(x_mat[,2], rev(x_mat[,3]))
    y2<-  c(x_mat[,5], rev(x_mat[,6]))
    polygon(x, y1, col = adjustcolor('darkcyan', alpha = 0.1), border=NA) # Add uncertainty polygon
    polygon(x, y2, col = adjustcolor('darksalmon', alpha = 0.2), border=NA) # Add uncertainty polygon

  lines(ts1$mean_abund~ts1$YEAR,col=adjustcolor('navy',alpha.f=0.5),lwd=2)
  points(ts1$mean_abund~ts1$YEAR,col='white',pch=21,bg=adjustcolor('navy',alpha.f=0.5),cex=1.5)
  lines(y_mat$median.rvc~c(seq(1993,2018)),col='dodgerblue4',lwd=2)
  points(y_mat$median.rvc~c(seq(1993,2018)),col='white',pch=21,bg='dodgerblue4',cex=1.5)
  
  lines(ts2$mean_abund~c(seq(1993,2018)),col=adjustcolor('darkred',alpha.f=0.5),lwd=2)
  points(ts2$mean_abund~c(seq(1993,2018)),col='white',pch=21,bg=adjustcolor('darkred',alpha.f=0.5),cex=1.5)
  lines(y_mat$median.reef~c(seq(1993,2018)),col='firebrick4',lwd=2)
  points(y_mat$median.reef~c(seq(1993,2018)),col='white',pch=21,bg='firebrick4',cex=1.5)
  
  legend(2013,c(max(na.omit(c(ts1$mean_abund,ts2$mean_abund)))*1.15),c('Obs. RVC surveys','Obs. REEF surveys','Est. RVC surveys','Est. REEF surveys'),text.col=c(adjustcolor('navy',alpha.f=0.5),adjustcolor('darkred',alpha.f=0.5),'dodgerblue4','firebrick4'),bty='n')
#  dev.off(paste(paste(sp,GZ,sep='_'),'.pdf',sep=''))
}


####Stan occupancy plot function####
TS_stan_plot_MARSS<- function(ts1,ts2,sp,GZ,params){
  
  x_mat<- data.frame(median.rvc=NA,l.95=NA,u.95=NA,median.reef=NA,l.95.rf=NA,u.95.rf=NA)
  for(i in 1:26){
    x_mat[i,1]=plogis(median(params$x1[,i]))
    x_mat[i,2]=plogis(quantile(params$x1[,i],0.1))
    x_mat[i,3]=plogis(quantile(params$x1[,i],0.9))
    x_mat[i,4]=plogis(median(params$x2[,i]))
    x_mat[i,5]=plogis(quantile(params$x2[,i],0.1))
    x_mat[i,6]=plogis(quantile(params$x2[,i],0.9))
  }
  
  y_mat_rvc<- data.frame(year=ts1$YEAR[complete.cases(ts1)],median.rvc=NA,l.95.rvc=NA,u.95.rvc=NA)
  y_mat_reef<- data.frame(year=seq(1993,2018),median.reef=NA,l.95.rf=NA,u.95.rf=NA)
  for(i in 1:23){
    y_mat_rvc[i,2]=plogis(median(params$a_yr1[,i]))
    y_mat_rvc[i,3]=plogis(quantile(params$a_yr1[,i],0.1))
    y_mat_rvc[i,4]=plogis(quantile(params$a_yr1[,i],0.9))
  }
  
  for(i in 1:26){
    y_mat_reef[i,2]=plogis(median(params$a_yr2[,i]))
    y_mat_reef[i,3]=plogis(quantile(params$a_yr2[,i],0.1))
    y_mat_reef[i,4]=plogis(quantile(params$a_yr2[,i],0.9))
  }
  y_mat<- full_join(y_mat_reef,y_mat_rvc)
  
  
  pdf(paste(paste(i,sp,GZ,sep='_'),'.pdf',sep=''),width=8,height=6)
  par(xpd=T)
  plot(ts1$p.occ~c(seq(1993,2018)),type='n',ylim=c(min(na.omit(c(y_mat$l.95.rf,y_mat$l.95.rvc))),max(na.omit(c(y_mat$u.95.rf,y_mat$u.95.rvc)))),col='darkblue',bty='l',ylab=expression('Mean count per survey'),xlab='Year',main=paste(sp,GZ,sep=', '))
  lines(x_mat[,1]~c(seq(1993,2018)),lty=5,lwd=2,col='darkcyan')
  lines(x_mat[,4]~c(seq(1993,2018)),lty=5,lwd=2,col='darksalmon')
  x<- c(c(seq(1993,2018)), rev(c(seq(1993,2018))))
  y1<- c(x_mat[,2], rev(x_mat[,3]))
  y2<-  c(x_mat[,5], rev(x_mat[,6]))
  polygon(x, y1, col = adjustcolor('darkcyan', alpha = 0.1), border=NA) # Add uncertainty polygon
  polygon(x, y2, col = adjustcolor('darksalmon', alpha = 0.2), border=NA) # Add uncertainty polygon
  
  lines(ts1$p.occ~ts1$YEAR,col=adjustcolor('navy',alpha.f=0.5),lwd=2)
  points(ts1$p.occ~ts1$YEAR,col='white',pch=21,bg=adjustcolor('navy',alpha.f=0.5),cex=1.5)
  lines(y_mat$median.rvc~c(seq(1993,2018)),col='dodgerblue4',lwd=2)
  points(y_mat$median.rvc~c(seq(1993,2018)),col='white',pch=21,bg='dodgerblue4',cex=1.5)
  
  lines(ts2$p.occ~c(seq(1993,2018)),col=adjustcolor('darkred',alpha.f=0.5),lwd=2)
  points(ts2$p.occ~c(seq(1993,2018)),col='white',pch=21,bg=adjustcolor('darkred',alpha.f=0.5),cex=1.5)
  lines(y_mat$median.reef~c(seq(1993,2018)),col='firebrick4',lwd=2)
  points(y_mat$median.reef~c(seq(1993,2018)),col='white',pch=21,bg='firebrick4',cex=1.5)
  
  legend(2013,c(max(na.omit(c(y_mat$u.95.rf,y_mat$u.95.rvc)))*1.15),c('Obs. RVC surveys','Obs. REEF surveys','Est. RVC surveys','Est. REEF surveys'),text.col=c(adjustcolor('navy',alpha.f=0.5),adjustcolor('darkred',alpha.f=0.5),'dodgerblue4','firebrick4'),bty='n')
  dev.off(paste(paste(sp,GZ,sep='_'),'.pdf',sep=''))
}

####Data####
###REEF data
REEF<- read.csv("./REEF Tropical Western Atlantic/TWA.csv") #Reef sightings 
REEF_survey<- read.csv("./REEF Tropical Western Atlantic/TWAsurveys.csv") #Survey metadata
m<- match(REEF$formid,REEF_survey$formid) #match sightings to surveys
REEF[,14:30]<- REEF_survey[m,4:20] #merge survey level data

###REEF geog
library(sp)
reef_geog<- read.csv("./REEF Tropical Western Atlantic/TWAgeog.csv")
site_surveys<-REEF %>% group_by(geogr) %>% summarize(n=n_distinct(formid))
m<- match(reef_geog$geogid,site_surveys$geogr)
reef_geog$no.surveys<- site_surveys$n[m]
reef_geog$region.id<- substr(reef_geog$geogid, 1,4) #Get the region id (first four digits)
reef_geog$lat_full<- reef_geog$lat #Copy of the full latitude
reef_geog$lon_full<- reef_geog$lon #Copy of the full longitude
reef_geog<-reef_geog%>% #Separate out degrees and minutes
  separate(lat,into=c("lat_deg","lat_min"),sep=" ")%>%
  separate(lon,into=c("lon_deg","lon_min"),sep=" ")

#convert desired columns numeric
col.num<-c("lat_deg","lat_min","lon_deg","lon_min")
reef_geog[col.num]<-sapply(reef_geog[col.num],as.numeric)
reef_geog<-reef_geog[complete.cases(reef_geog),] #Discard sites with no coordinates

#Convert to decimal degrees (matches RVC lat/lons)
reef_geog$lat_dd<- reef_geog$lat_deg+reef_geog$lat_min/60
reef_geog$lon_dd<- reef_geog$lon_deg-reef_geog$lon_min/60

#Remove these surveys
R<- subset(REEF, type==1) #Remove species-only surveys

# Get rid of  rows of data with no date reported
R<-R[-which(R$date=="0000-00-00"),]

# let's get the year & month from the Date field (requires lubridate)
R$date<-ymd(R$date)#put into proper date format
R<-cbind(R,year=year(R$date))
R<-cbind(R,month=month(R$date))
R<-cbind(R,day=day(R$date))

# Remove sites without lat-long (b/c those sites are probably janky anyway)
R<-filter(R, R$geogr %in% reef_geog$geogid)

# Thin  raw data to exclude dives shorter than 20min, longer than 120min
R<-R[as.numeric(R$btime)>20,]
R<-R[as.numeric(R$btime)<120,]

# Thin the data to remove night dives (start before 5am or after 8pm)
R<-R[as.numeric(R$start)>5,]
R<-R[as.numeric(R$start)<20,]

#RVC
fk_99_18<- read.csv("./RVC/Florida_Keys_RVC.csv") #RVC survey data from 1999 to 2018
fk_79_98<- read.csv("./RVC/Pre1999_edit.txt") #RVC survey dating from 1979 to 1998

#Synonymize column names
colnames(fk_79_98)<- c(colnames(fk_99_18)[1],colnames(fk_99_18)[3:5],colnames(fk_99_18)[2],colnames(fk_99_18)[6],colnames(fk_99_18)[7:18],colnames(fk_99_18)[19],colnames(fk_99_18)[20]) 
#Join together datasets
fk_79_18<- full_join(fk_99_18,fk_79_98)


library(sp)
reef_geog<- read.csv("./REEF Tropical Western Atlantic/TWAgeog.csv")
site_surveys<-REEF %>% group_by(geogr) %>% summarize(n=n_distinct(formid))
m<- match(reef_geog$geogid,site_surveys$geogr)
reef_geog$no.surveys<- site_surveys$n[m]
reef_geog$region.id<- substr(reef_geog$geogid, 1,4) #Get the region id (first four digits)
reef_geog$lat_full<- reef_geog$lat #Copy of the full latitude
reef_geog$lon_full<- reef_geog$lon #Copy of the full longitude
reef_geog<-reef_geog%>% #Separate out degrees and minutes
  separate(lat,into=c("lat_deg","lat_min"),sep=" ")%>%
  separate(lon,into=c("lon_deg","lon_min"),sep=" ")

#convert desired columns numeric
col.num<-c("lat_deg","lat_min","lon_deg","lon_min")
reef_geog[col.num]<-sapply(reef_geog[col.num],as.numeric)
reef_geog<-reef_geog[complete.cases(reef_geog),] #Discard sites with no coordinates

#Convert to decimal degrees (matches RVC lat/lons)
reef_geog$lat_dd<- reef_geog$lat_deg+reef_geog$lat_min/60
reef_geog$lon_dd<- reef_geog$lon_deg-reef_geog$lon_min/60

#Extract out florida key subregions for select zones
reef_3403_sites<- subset(reef_geog,region.id=='3403')

#RVC sites
fk_79_18$LAT_LON<- paste(fk_79_18$LAT_DEGREES,fk_79_18$LON_DEGREES,sep='_') #Find unique geographic position
rvc_sites<- data.frame(lat=fk_79_18$LAT_DEGREES,lon=fk_79_18$LON_DEGREES,lat_lon=fk_79_18$LAT_LON) 
rvc_sites<- distinct(rvc_sites,lat_lon,.keep_all = T) #keep unique sample sites from the RVC surveys

#Match up RVC sites to REEF sites
library(rgeos)
set1 <- SpatialPoints(cbind(rvc_sites$lat,rvc_sites$lon)) #Set of RVC site points
set2 <- SpatialPoints(cbind(reef_geog$lat_dd,reef_geog$lon_dd)) #Set of REEF site points
matched_reef_site<- apply(gDistance(set1, set2, byid=TRUE), 2, which.min) #Find closest REEF site by geographic distance between point sets
rvc_sites$REEF_site <- reef_geog$geog[matched_reef_site]
rvc_sites$REEF_lat<- reef_geog$lat_dd[matched_reef_site]
rvc_sites$REEF_lon<- reef_geog$lon_dd[matched_reef_site]
rvc_sites$region.id<- reef_geog$region.id[matched_reef_site]
fk_79_18[,25:28]<- rvc_sites[match(fk_79_18$LAT_LON,rvc_sites$lat_lon),4:7]


###3. Geographic filtering ###
### Intersect REEF & RVC spatial data ###
library(sf)
rvc_grid<-st_read(dsn='RVC Grid', layer='FlaKeys_Grid') #Read in Florida Keys sampling grid
rvc_grid_84<- st_transform(rvc_grid,4326) #Set to WGS84

#3403 regions
reef_pts<-st_as_sf(x = reef_geog, 
                   coords = c("lon_dd", "lat_dd"),
                   crs = 4326)

rvc_pts<-st_as_sf(x = rvc_sites, 
                  coords = c("lon", "lat"),
                  crs = 4326)

grid_match<- st_intersects(reef_pts,rvc_grid_84,sparse=T)
reef_geog$grid_match<- NA
for(i in 1:nrow(reef_pts)){
  if(length(grid_match[[i]])==1){
    reef_geog$grid_match[i]=grid_match[[i]]  
  }else{
    reef_geog$grid_match[i]=NA
  }
}
reef_geog$hab_class<- as.factor(rvc_grid$habclass[reef_geog$grid_match])

grid_match_rvc<- st_intersects(rvc_pts,rvc_grid_84,sparse=T)
rvc_sites$grid_match<- NA
for(i in 1:nrow(rvc_sites)){
  if(length(grid_match_rvc[[i]])==1){
    rvc_sites$grid_match[i]=grid_match_rvc[[i]]  
  }else{
    rvc_sites$grid_match[i]=NA
  }
}
rvc_sites$hab_class<- as.factor(rvc_grid$habclass[rvc_sites$grid_match])

reef_geog_3403<- reef_geog %>% subset(is.na(grid_match)==F & region.id==3403)

### 4. Creating RVC time-series ###
#Fish data from REEF - remove ultra rare and basket species designations
fish_reef<- read.csv("Caribbean_fish_trait_matrix.csv") #fish species found in the Tropical Western Atlantic
fish_reef<- subset(fish_reef,expert_sighting_freq>1) #Take out the very rare species
fish_rvc<- read.csv("Florida_keys_taxonomic_data.csv")
fish_rvc<- subset(fish_rvc,gsub('.*\\ ', '', fish_rvc$SCINAME)!='sp.') #remove unknown species
fish_rvc<- subset(fish_rvc, SCINAME %in% fish_reef$sciname2)
m<- match(fish_reef$sciname2,fish_rvc$SCINAME)
fish_reef$rvc_code<- fish_rvc$SPECIES_CD[m]

fk_93_18<- subset(fk_79_18,YEAR>=1993) #Subset for the dataset from 1993 to match the first year of REEF surveys

rvc_occs_1<- rvc_filter(fk_93_18,GZ='3403',sp=fish_reef)
rvc_ts<- ts_rvc(rvc_occs_1)
rvc_ts_filter<- rlist::list.filter(rvc_ts,length(na.omit(p.occ))>18)

rvc.green<- do.call(rbind, lapply(rvc_ts_filter, data.frame, stringsAsFactors=FALSE))
rvc.green.sp<- unique(rvc.green$sp)
fish_reef_trim<- subset(fish_reef, rvc_code %in% rvc.green.sp)

rvc_occs<- rvc_filter(fk_93_18,GZ='3403',sp=fish_reef_trim)
reef_occs<- reef_filter(R,GZ='3403',sp=fish_reef_trim,geog=reef_geog_3403)
reef_ts<- ts_reef(reef_occs,sp=fish_reef_trim)
rvc_ts<- ts_rvc(rvc_occs)

####Stan models####
nb_test_rvc<-"data{
  int<lower=1> N;//number of observations (SSU surveys)
  int y[N]; //presence or absence on each survey
  int<lower=0> N_hab; //number of habitat classes
  int<lower=1,upper=N_hab> hab_class[N]; // vector of habitat class identities
  int<lower=0> N_yr; //number of years
  int<lower=1,upper=N_yr> year[N]; // vector of year
  int K; // columns in the covariate matrix
  matrix[N,K] X; // design matrix X
}
parameters{
  real alpha; 
  //global intercept
  real<lower = 0> recip_phi;
  //dispersion parameter
  vector[K] beta;
  //depth covariate
  
  //deviations from intercept
  vector[N_hab] z_hab; //deviation between habitats
  vector[N_yr] z_yr; //deviation between years
 
  //st dev on the deviations
  real<lower = 0> tau_hab; //sigma on habitat
  real<lower = 0> tau_yr; //sigma on year
}
transformed parameters{
  real phi;

  phi = 1/recip_phi;
}  
model{
  //priors
  alpha ~ normal(0,3);
  beta ~ normal(0,2);
  
  //shape parameter
  recip_phi ~ cauchy(0, 5);
  
  //varying intercepts
  z_hab ~ normal(0,3);
  z_yr ~ student_t(5, 0, 3);
  tau_yr ~ cauchy(0,3);
  tau_hab ~ cauchy(0,3);
  
  y ~ neg_binomial_2_log(alpha + z_yr[year]*tau_yr + z_hab[hab_class]*tau_hab + X*beta,phi);
}
"



nb_test_SS_rvc<-"data{
  int TT; //timespan
  int N;//number of observations (SSU surveys)
  int y[N]; //presence or absence on each survey
  int N_hab; //number of habitat classes
  int<lower=1,upper=N_hab> hab_class[N]; // vector of habitat class identities
  int N_yr; //number of years
  int yr_index[N_yr]; //index of years
  int<lower=1,upper=N_yr> year_id[N]; // vector of year ids
  int K; // columns in the covariate matrix
  matrix[N,K] X; // design matrix X
}
parameters {
  //global intercept
  real x0;
  vector[K] beta;
  real<lower = 0> recip_phi;
  
  //deviations from intercept
  vector[N_hab] a_hab; //deviation between habitats
  vector[TT] pro_dev; //state deviations
  vector[TT] obs_dev; //state deviations

  
  //st dev on the deviations
  real<lower = 0> sd_hab; //sigma on habitat
  real<lower = 0> sd_r; //sigma for observation error
  real<lower = 0> sd_q; //sigma on process error
}

transformed parameters{
  real phi;
  vector[TT] x;
  vector[N_yr] a_yr;
  
  phi=1/recip_phi;
  
  x[1] = x0 + pro_dev[1];
   
  for(t in 2:TT){
    x[t] = x[t-1] + pro_dev[t];
  }
   
  for(i in 1:N_yr){
    a_yr[i] = x[yr_index[i]] + obs_dev[i]; 
  }
  
}  
  
model{
  //priors
  beta ~ normal(0,2);
 
  //shape parameter
  recip_phi ~ cauchy(0, 10);
 
  //standard deviations
  sd_hab ~ cauchy(0,1);
  sd_r ~ cauchy(0,0.5);
  sd_q ~ cauchy(0,0.5);
  
  //varying intercepts
  a_hab ~ normal(0,sd_hab);
  
  x0 ~ normal(0,3);
  
  for(t in 1:TT){
    pro_dev[t] ~ normal(0, sd_q);
    obs_dev[t] ~ normal(0,sd_r);
  }
 
  //draw presences
    y ~ neg_binomial_2_log(a_yr[year_id] + a_hab[hab_class] + X*beta,phi);
}
"

nb_test_SS_rvc2<-"data{
  int TT; //timespan
  int N;//number of observations (SSU surveys)
  int y[N]; //presence or absence on each survey
  int N_hab; //number of habitat classes
  int<lower=1,upper=N_hab> hab_class[N]; // vector of habitat class identities
  int N_yr; //number of years
  int yr_index[N_yr]; //index of years
  int<lower=1,upper=N_yr> year_id[N]; // vector of year ids
  int K; // columns in the covariate matrix
  matrix[N,K] X; // design matrix X
}
parameters {
  //global intercept
  real x0;
  vector[K] beta;
  vector[T]
  real<lower = 0> recip_phi;
  
  //deviations from intercept
  vector[N_hab] z_hab; //deviation between habitats
  vector[TT] pro_dev; //state deviations
  vector[TT] obs_dev; //state deviations

  
  //st dev on the deviations
  real<lower = 0> sd_hab; //sigma on habitat
  real<lower = 0> sd_r; //sigma for observation error
  real<lower = 0> sd_q; //sigma on process error
}

transformed parameters{
  real phi;
  vector[TT] x;
  vector[N_yr] a_yr;
  
  phi=1/recip_phi;
  
  x[1] = x0 + pro_dev[1];
   
  for(t in 2:TT){
    x[t] = x[t-1] + pro_dev[t];
  }
   
  for(i in 1:N_yr){
    a_yr[i] = x[yr_index[i]] + obs_dev[i]; 
  }
  
}  
  
model{
  //priors
  beta ~ normal(0,2);
 
  //shape parameter
  recip_phi ~ cauchy(0, 10);
 
  //standard deviations
  sd_hab ~ cauchy(0,1);
  sd_r ~ cauchy(0,1);
  sd_q ~ cauchy(0,0.5);
  
  //varying intercepts
  z_hab ~ normal(0,5);
  
  x0 ~ normal(0,5);
  
  for(t in 1:TT){
    pro_dev[t] ~ normal(0, sd_q);
    obs_dev[t] ~ normal(0,sd_r);
  }
 
  //draw presences
    y ~ neg_binomial_2_log(a_yr[year_id] + z_hab[hab_class]*sd_hab + X*beta,phi);
}
"

ord_test_reef<-"
functions {
  real induced_dirichlet_lpdf(vector c, vector alpha, real phi) {
    int K = num_elements(c) + 1;
    vector[K - 1] sigma = inv_logit(phi - c);
    vector[K] p;
    matrix[K, K] J = rep_matrix(0, K, K);
    
    // Induced ordinal probabilities
    p[1] = 1 - sigma[1];
    for (k in 2:(K - 1))
      p[k] = sigma[k - 1] - sigma[k];
    p[K] = sigma[K - 1];
    
    // Baseline column of Jacobian
    for (k in 1:K) J[k, 1] = 1;
    
    // Diagonal entries of Jacobian
    for (k in 2:K) {
      real rho = sigma[k - 1] * (1 - sigma[k - 1]);
      J[k, k] = - rho;
      J[k - 1, k] = rho;
    }
    
    return   dirichlet_lpdf(p | alpha)
           + log_determinant(J);
  }
}
data{
  int<lower=1> N;//number of observations (SSU surveys)
  int y[N]; //abundance category for each survey
  int<lower=0> N_hab; //number of habitat classes
  int<lower=1,upper=N_hab> hab_class[N]; // vector of habitat class identities
  int<lower=0> N_yr; //number of years
  int<lower=1,upper=N_yr> year[N]; // vector of year
  int<lower=0> N_site; //number of sites
  int<lower=1,upper=N_site> site[N]; // vector of site identities
  int<lower=0> N_dv; //number of divers
  int<lower=1,upper=N_dv> diver[N]; // vector of diver identities
  int<lower=0> N_dmy; //number of site day clusters
  int<lower=1,upper=N_dmy> dmy[N]; // vector of site day cluster identities
  int Z; // columns in the covariate matrix
  matrix[N,Z] X; // design matrix X
  int K; //ordinal levels
}
parameters {
  ordered[K-1] c; //cutpoints
 
  //deviations from intercept
  vector[Z] beta; //effort coefficients
  vector[N_hab] a_hab; //deviation between habitats
  vector[N_yr] a_yr; //deviation between years
  vector[N_site] a_site; //deviation between sites
  vector[N_dv] a_dv; //deviation between divers
  vector[N_dmy] a_dmy; //deviation between site day clusters
 
  //st dev on the deviations
  real<lower = 0> sigma_hab;
  real<lower = 0> sigma_yr;
  real<lower = 0> sigma_site;
  real<lower = 0> sigma_dv;
  real<lower = 0> sigma_dmy;
}

model{
  //priors
  c ~ induced_dirichlet(rep_vector(1, K), 0);
  beta ~ normal(0,2);
  
  //standard deviations

  sigma_hab ~ cauchy(0, 3);
  sigma_yr ~ cauchy(0, 3);
  sigma_site ~ cauchy(0, 3);
  sigma_dv ~ cauchy(0, 3);
  sigma_dmy ~ cauchy(0, 3);
  
  //varying intercepts
  a_hab ~ normal(0, sigma_hab);
  a_yr ~ normal(0, sigma_yr);
  a_site ~ normal(0, sigma_site);
  a_dv ~ normal(0, sigma_dv);
  a_dmy ~ normal(0, sigma_dmy);

  y ~ ordered_logistic(a_hab[hab_class]+a_yr[year]+a_site[site]+a_dv[diver]+a_dmy[dmy]+X*beta,c);
  
}
"

ord_test_SS_reef<-"
functions {
  real induced_dirichlet_lpdf(vector c, vector alpha, real phi) {
    int K = num_elements(c) + 1;
    vector[K - 1] sigma = inv_logit(phi - c);
    vector[K] p;
    matrix[K, K] J = rep_matrix(0, K, K);
    
    // Induced ordinal probabilities
    p[1] = 1 - sigma[1];
    for (k in 2:(K - 1))
      p[k] = sigma[k - 1] - sigma[k];
    p[K] = sigma[K - 1];
    
    // Baseline column of Jacobian
    for (k in 1:K) J[k, 1] = 1;
    
    // Diagonal entries of Jacobian
    for (k in 2:K) {
      real rho = sigma[k - 1] * (1 - sigma[k - 1]);
      J[k, k] = - rho;
      J[k - 1, k] = rho;
    }
    
    return   dirichlet_lpdf(p | alpha)
           + log_determinant(J);
  }
}
data{
  int<lower=1> N;//number of observations (SSU surveys)
  int y[N]; //abundance category for each survey
  int<lower=0> N_hab; //number of habitat classes
  int<lower=1,upper=N_hab> hab_class[N]; // vector of habitat class identities
  int<lower=0> N_site; //number of sites
  int<lower=1,upper=N_site> site[N]; // vector of site identities
  int<lower=0> N_dv; //number of divers
  int<lower=1,upper=N_dv> diver[N]; // vector of diver identities
  int<lower=0> N_dmy; //number of site day clusters
  int<lower=1,upper=N_dmy> dmy[N]; // vector of site day cluster identities
  int Z; // columns in the covariate matrix
  matrix[N,Z] X; // design matrix X
  int K; //ordinal levels
  int TT; // timespan
  int<lower=0> N_yr; //number of years
  int yr_index[N_yr]; //index of years
  int<lower=1,upper=N_yr> year_id[N]; // vector of year
}
parameters {
  ordered[K-1] c; //cutpoints
  real x0; //initial popn size
  
  //deviations from intercept
  vector[Z] beta; //effort coefficients
  vector[N_hab] a_hab; //deviation between habitats
  vector[N_site] a_site; //deviation between sites
  vector[N_dv] a_dv; //deviation between divers
  vector[N_dmy] a_dmy; //deviation between site day clusters
 
  //st dev on the deviations
  real<lower = 0> sd_hab;
  real<lower = 0> sd_site;
  real<lower = 0> sd_dv;
  real<lower = 0> sd_dmy;
  real<lower = 0> sd_r;
  real<lower = 0> sd_q;
  
  vector[TT] pro_dev; //process deviations
  vector[N_yr] obs_dev; //observation deviations 
  
}

transformed parameters{
  vector[TT] x;
  vector[N_yr] a_yr;
  x[1] = x0 + pro_dev[1];
   
  for(t in 2:TT){
    x[t] = x[t-1] + pro_dev[t];
  }
   
  for(i in 1:N_yr){
    a_yr[i] = x[yr_index[i]] + obs_dev[i]; 
  }
  
}  

model{
  //priors
  c ~ induced_dirichlet(rep_vector(1, K), 0);
  beta ~ normal(0,2);
  x0 ~ normal(0,5);
  
  //standard deviations

  sd_hab ~ cauchy(0, 3);
  sd_q ~ cauchy(0, 3);
  sd_r ~ cauchy(0, 3);
  sd_site ~ cauchy(0, 3);
  sd_dv ~ cauchy(0, 3);
  sd_dmy ~ cauchy(0, 3);
  
  //varying intercepts
  a_hab ~ normal(0, sd_hab);
  a_site ~ normal(0, sd_site);
  a_dv ~ normal(0, sd_dv);
  a_dmy ~ normal(0, sd_dmy);
  
  for(t in 1:TT){
    pro_dev[t] ~ normal(0, sd_q);
    obs_dev[t] ~ normal(0,sd_r);
  }

  y ~ ordered_logistic(a_hab[hab_class]+a_yr[year_id]+a_site[site]+a_dv[diver]+a_dmy[dmy]+X*beta,c);
  
}
"

ord_test_SS_reef2<-"
functions {
  real induced_dirichlet_lpdf(vector c, vector alpha, real phi) {
    int K = num_elements(c) + 1;
    vector[K - 1] sigma = inv_logit(phi - c);
    vector[K] p;
    matrix[K, K] J = rep_matrix(0, K, K);
    
    // Induced ordinal probabilities
    p[1] = 1 - sigma[1];
    for (k in 2:(K - 1))
      p[k] = sigma[k - 1] - sigma[k];
    p[K] = sigma[K - 1];
    
    // Baseline column of Jacobian
    for (k in 1:K) J[k, 1] = 1;
    
    // Diagonal entries of Jacobian
    for (k in 2:K) {
      real rho = sigma[k - 1] * (1 - sigma[k - 1]);
      J[k, k] = - rho;
      J[k - 1, k] = rho;
    }
    
    return   dirichlet_lpdf(p | alpha)
           + log_determinant(J);
  }
  
  transf_abundance(c,beta_yr,N_yr){
    matrix[N_yr,num_elements(c)+1] p;
    vector[N_yr] lambda;
    vector[num_elements(c)+1] abunds;
    
    if(num_elements(abunds)==3){
    abunds[1]=0
    abunds[2]=1
    abunds[3]=2
    }
    if(num_elements(abunds)==4){
    abunds[1]=0
    abunds[2]=1
    abunds[3]=2
    abunds[4]=11
    }
     if(num_elements(abunds)==5){
    abunds[1]=0
    abunds[2]=1
    abunds[3]=2
    abunds[4]=11
    abunds[5]=101
    }
    
    for(t in 1:N_yr){
      if(num_elements(p)==3){
        p[t,1] = inv_logit(c[1] + a_yr[t])
        p[t,2] = inv_logit(c[2] + a_yr[t]) - inv_logit(c[1] + a_yr[t])
        p[t,3] = 1-inv_logit(c[2] + a_yr[t])
      }
      if(num_elements(p)==4){
        p[t,1] = inv_logit(c[1] + a_yr[t])
        p[t,2] = inv_logit(c[2] + a_yr[t]) - inv_logit(c[1] + a_yr[t])
        p[t,3] = inv_logit(c[3] + a_yr[t]) - inv_logit(c[2] + a_yr[t])
        p[t,4] = 1-inv_logit(c[3] + a_yr[t])
      }
      
    }
    if(num_elements(p)==5){
        p[t,1] = inv_logit(c[1] + a_yr[t])
        p[t,2] = inv_logit(c[2] + a_yr[t]) - inv_logit(c[1] + a_yr[t])
        p[t,3] = inv_logit(c[3] + a_yr[t]) - inv_logit(c[2] + a_yr[t])
        p[t,4] = inv_logit(c[4] + a_yr[t]) - inv_logit(c[3] + a_yr[t])
        p[t,5] = 1-inv_logit(c[4] + a_yr[t])
      }
  for(z in 1:N_yr){
    lambda[z]=p[z,]*abunds
  }
  return lambda
  }
}
data{
  int<lower=1> N;//number of observations (SSU surveys)
  int y[N]; //abundance category for each survey
  int<lower=0> N_hab; //number of habitat classes
  int<lower=1,upper=N_hab> hab_class[N]; // vector of habitat class identities
  int<lower=0> N_site; //number of sites
  int<lower=1,upper=N_site> site[N]; // vector of site identities
  int<lower=0> N_dv; //number of divers
  int<lower=1,upper=N_dv> diver[N]; // vector of diver identities
  int<lower=0> N_dmy; //number of site day clusters
  int<lower=1,upper=N_dmy> dmy[N]; // vector of site day cluster identities
  int Z; // columns in the covariate matrix
  matrix[N,Z] X; // design matrix X
  int K; //ordinal levels
  int TT; // timespan
  int<lower=0> N_yr; //number of years
  int yr_index[N_yr]; //index of years
  matrix[N,N_yr]; //year matrix - fixed effect
}
parameters {
  ordered[K-1] c; //cutpoints
  real x0; //initial popn size
  
  //deviations from intercept
  vector[Z] beta; //effort coefficients
  vector[N_yr] beta_yr; //year coefficients
  vector[N_hab] a_hab; //deviation between habitats
  vector[N_site] a_site; //deviation between sites
  vector[N_dv] a_dv; //deviation between divers
  vector[N_dmy] a_dmy; //deviation between site day clusters
 
  //st dev on the deviations
  real<lower = 0> sd_hab;
  real<lower = 0> sd_site;
  real<lower = 0> sd_dv;
  real<lower = 0> sd_dmy;
  real<lower = 0> sd_r;
  real<lower = 0> sd_q;
  
  vector[TT] pro_dev; //process deviations
  vector[N_yr] obs_dev; //observation deviations 
  
}

transformed parameters{
  vector[N_yr] lambda_y; // estimate of annual abundance
  vector[TT] x; // unobserved state for abundance
  vector[TT] pro_dev; //process deviations
  vector[N_yr] obs_dev; //observation deviations 
 
  lambda_y = transf_abundance(c,a_yr,N_yr)
 
  x[1] = x0 + pro_dev[1];
   
  for(t in 2:TT){
    x[t] = x[t-1] + pro_dev[t];
  }
   
  for(t in 1:TT){
    lambda_y[t] = x[t] + obs_dev[t];
  }
   
  
}  

model{
  //priors
  c ~ induced_dirichlet(rep_vector(1, K), 0);
  beta ~ normal(0,2);
  beta_yr ~ normal(0,5);
  x0 ~ normal(0,5);
  
  //standard deviations

  sd_hab ~ cauchy(0, 3);
  sd_q ~ cauchy(0, 3);
  sd_r ~ cauchy(0, 3);
  sd_site ~ cauchy(0, 3);
  sd_dv ~ cauchy(0, 3);
  sd_dmy ~ cauchy(0, 3);
  
  //varying intercepts
  a_hab ~ normal(0, sd_hab);
  a_site ~ normal(0, sd_site);
  a_dv ~ normal(0, sd_dv);
  a_dmy ~ normal(0, sd_dmy);
  
  for(t in 1:TT){
    pro_dev[t] ~ normal(0, sd_q);
    obs_dev[t] ~ normal(0,sd_r);
  }

  y ~ ordered_logistic(a_yr[year_id]+a_hab[hab_class]+a_site[site]+a_dv[diver]+a_dmy[dmy]+X*beta+Y*beta_yr,c);
  
  
}
"


abund_test_SS_comb<-"functions {
  real induced_dirichlet_lpdf(vector c, vector alpha, real phi) {
    int K = num_elements(c) + 1;
    vector[K - 1] sigma = inv_logit(phi - c);
    vector[K] p;
    matrix[K, K] J = rep_matrix(0, K, K);
    
    // Induced ordinal probabilities
    p[1] = 1 - sigma[1];
    for (k in 2:(K - 1))
      p[k] = sigma[k - 1] - sigma[k];
    p[K] = sigma[K - 1];
    
    // Baseline column of Jacobian
    for (k in 1:K) J[k, 1] = 1;
    
    // Diagonal entries of Jacobian
    for (k in 2:K) {
      real rho = sigma[k - 1] * (1 - sigma[k - 1]);
      J[k, k] = - rho;
      J[k - 1, k] = rho;
    }
    
    return   dirichlet_lpdf(p | alpha)
           + log_determinant(J);
  }
}
data{
  int<lower=1> N1;//number of observations (REEF surveys)
  int<lower=1> N2;//number of observations (REEF surveys)
  int y1[N1]; //abundance category for each survey
  int y2[N2]; //abundance category for each survey
  int<lower=0> N_hab1; //number of habitat classes
  int<lower=1,upper=N_hab1> hab_class1[N1]; // vector of habitat class identities
  int<lower=0> N_hab2; //number of habitat classes
  int<lower=1,upper=N_hab2> hab_class2[N2]; // vector of habitat class identities
  int<lower=0> N_site; //number of sites
  int<lower=1,upper=N_site> site[N2]; // vector of site identities
  int<lower=0> N_dv; //number of divers
  int<lower=1,upper=N_dv> diver[N2]; // vector of diver identities
  int<lower=0> N_dmy; //number of site day clusters
  int<lower=1,upper=N_dmy> dmy[N2]; // vector of site day cluster identities
  int Z1; // columns in the covariate matrix
  int Z2; // columns in the covariate matrix
  matrix[N1,Z1] X1; // design matrix X
  matrix[N2,Z2] X2; // design matrix X
  int K; //ordinal levels
  int TT; // timespan
  int<lower=0> N_yr1; //number of years
  int yr_index1[N_yr1]; //index of years
  int<lower=1,upper=N_yr1> year_id1[N1]; // vector of year
  int<lower=0> N_yr2; //number of years
  int yr_index2[N_yr2]; //index of years
  int<lower=1,upper=N_yr2> year_id2[N2]; // vector of year

}
parameters {
  ordered[K-1] c; //cutpoints
  real x0; //initial popn size
  real<lower = 0> recip_phi; //overdispersion parameter
  real a; //scalar for time-series 2

  //deviations from intercept
  vector[Z1] beta1; //effort coefficients - RVC
  vector[Z2] beta2; //effort coefficients - REEF
  vector[N_hab1] a_hab1; //deviation between habitats - RVC
  vector[N_hab2] a_hab2; //deviation between habitats - REEF
  vector[N_site] a_site; //deviation between sites
  vector[N_dv] a_dv; //deviation between divers
  vector[N_dmy] a_dmy; //deviation between site day clusters
 
  //variance on the deviance components
  real<lower = 0> sd_hab1;
  real<lower = 0> sd_hab2;
  real<lower = 0> sd_site;
  real<lower = 0> sd_dv;
  real<lower = 0> sd_dmy;
  real<lower = 0> sd_r1;
  real<lower = 0> sd_r2;
  real<lower = 0> sd_q;
  
  vector[TT] pro_dev; //process deviations
  vector[N_yr1] obs_dev1; //observation deviations 
  vector[N_yr2] obs_dev2; //observation deviations 
  
}

transformed parameters{
  vector[TT] x;
  vector[N_yr1] a_yr1;
  vector[N_yr2] a_yr2;
  real phi;
 
  phi=1/recip_phi;
  
  x[1] = x0 + pro_dev[1];
   
  for(t in 2:TT){
    x[t] = x[t-1] + pro_dev[t];
  }
   
  for(i in 1:N_yr1){
    a_yr1[i] = x[yr_index1[i]] + obs_dev1[i]; 
  }
    for(i in 1:N_yr2){
      a_yr2[i] = x[yr_index2[i]] + a + obs_dev2[i]; 
}
  
}  

model{
  //priors
  c ~ induced_dirichlet(rep_vector(1, K), 0); //prior on ordered cut-points
  beta1 ~ normal(0,2); //covariates - rvc
  beta2 ~ normal(0,2); //covariates - reef
  x0 ~ normal(0,5); //initial state
  a ~ normal(0,5); //scalar
  
  //variance terms
  sd_hab1 ~ inv_gamma(2, 1);
  sd_hab2 ~ inv_gamma(2, 1);
  sd_q ~inv_gamma(2,0.25);
  sd_r1 ~ inv_gamma(2,0.25);
  sd_r2 ~ inv_gamma(2,0.25);
  sd_site ~ inv_gamma(2, 1);
  sd_dv ~ inv_gamma(2, 1);
  sd_dmy ~ inv_gamma(2, 1);
  
  //varying intercepts
  a_hab1 ~ normal(0, sd_hab1);
  a_hab2 ~ normal(0, sd_hab2);
  a_site ~ normal(0, sd_site);
  a_dv ~ normal(0, sd_dv);
  a_dmy ~ normal(0, sd_dmy);
  
  for(t in 1:N_yr1){
   obs_dev1[t] ~ normal(0,sd_r1);
  }
  
  for(t in 1:TT){
    pro_dev[t] ~ normal(0, sd_q);
    obs_dev2[t] ~ normal(0,sd_r2);
  }

  y1 ~ neg_binomial_2_log(a_yr1[year_id1] + a_hab1[hab_class1] + X1*beta1,phi);
  
  y2 ~ ordered_logistic(a_yr2[year_id2]+a_hab2[hab_class2]+a_site[site]+a_dv[diver]+a_dmy[dmy]+X2*beta2,c);
  
}
  generated quantities{
  vector[N1+N2] log_lik;
  for (i in 1:N1) log_lik[i] = neg_binomial_2_log_lpmf(y1[i]|a_yr1[year_id1[i]] + a_hab1[hab_class1[i]] + X1[i,]*beta1,phi);
  for (z in 1:N2) log_lik[N1+z] = ordered_logistic_lpmf(y2[z]|a_hab2[hab_class2[z]]+a_yr2[year_id2[z]]+a_site[site[z]]+a_dv[diver[z]]+a_dmy[dmy[z]]+X2[z,]*beta2, c);
} 
"

abund_test_SS_sep<-"functions {
  real induced_dirichlet_lpdf(vector c, vector alpha, real phi) {
    int K = num_elements(c) + 1;
    vector[K - 1] sigma = inv_logit(phi - c);
    vector[K] p;
    matrix[K, K] J = rep_matrix(0, K, K);
    
    // Induced ordinal probabilities
    p[1] = 1 - sigma[1];
    for (k in 2:(K - 1))
      p[k] = sigma[k - 1] - sigma[k];
    p[K] = sigma[K - 1];
    
    // Baseline column of Jacobian
    for (k in 1:K) J[k, 1] = 1;
    
    // Diagonal entries of Jacobian
    for (k in 2:K) {
      real rho = sigma[k - 1] * (1 - sigma[k - 1]);
      J[k, k] = - rho;
      J[k - 1, k] = rho;
    }
    
    return   dirichlet_lpdf(p | alpha)
           + log_determinant(J);
  }
}
data{
  int<lower=1> N1;//number of observations (REEF surveys)
  int<lower=1> N2;//number of observations (REEF surveys)
  int y1[N1]; //abundance category for each survey
  int y2[N2]; //abundance category for each survey
  int<lower=0> N_hab1; //number of habitat classes
  int<lower=1,upper=N_hab1> hab_class1[N1]; // vector of habitat class identities
  int<lower=0> N_hab2; //number of habitat classes
  int<lower=1,upper=N_hab2> hab_class2[N2]; // vector of habitat class identities
  int<lower=0> N_site; //number of sites
  int<lower=1,upper=N_site> site[N2]; // vector of site identities
  int<lower=0> N_dv; //number of divers
  int<lower=1,upper=N_dv> diver[N2]; // vector of diver identities
  int<lower=0> N_dmy; //number of site day clusters
  int<lower=1,upper=N_dmy> dmy[N2]; // vector of site day cluster identities
  int Z1; // columns in the covariate matrix
  int Z2; // columns in the covariate matrix
  matrix[N1,Z1] X1; // design matrix X
  matrix[N2,Z2] X2; // design matrix X
  int K; //ordinal levels
  int TT; // timespan
  int<lower=0> N_yr1; //number of years
  int yr_index1[N_yr1]; //index of years
  int<lower=1,upper=N_yr1> year_id1[N1]; // vector of year
  int<lower=0> N_yr2; //number of years
  int yr_index2[N_yr2]; //index of years
  int<lower=1,upper=N_yr2> year_id2[N2]; // vector of year

}
parameters {
  ordered[K-1] c; //cutpoints
  real x01; //initial popn size - rvc
  real x02; //initial popn size - reef
  real<lower = 0> recip_phi; //overdispersion parameter
 
  //deviations from intercept
  vector[Z1] beta1; //effort coefficients - RVC
  vector[Z2] beta2; //effort coefficients - REEF
  vector[N_hab1] a_hab1; //deviation between habitats - RVC
  vector[N_hab2] a_hab2; //deviation between habitats - REEF
  vector[N_site] a_site; //deviation between sites
  vector[N_dv] a_dv; //deviation between divers
  vector[N_dmy] a_dmy; //deviation between site day clusters
 
  //variance on the deviance components
  real<lower = 0> sd_hab1;
  real<lower = 0> sd_hab2;
  real<lower = 0> sd_site;
  real<lower = 0> sd_dv;
  real<lower = 0> sd_dmy;
  real<lower = 0> sd_r1;
  real<lower = 0> sd_r2;
  real<lower = 0> sd_q1;
  real<lower = 0> sd_q2;
  
  vector[TT] pro_dev1; //process deviations
  vector[TT] pro_dev2; //process deviations
  vector[N_yr1] obs_dev1; //observation deviations 
  vector[N_yr2] obs_dev2; //observation deviations 
  
}

transformed parameters{
  vector[TT] x1;
  vector[TT] x2;
  vector[N_yr1] a_yr1;
  vector[N_yr2] a_yr2;
  real phi;
 
  phi=1/recip_phi;
  
  x1[1] = x01 + pro_dev1[1];
  x2[1] = x02 + pro_dev2[1];
   
  for(t in 2:TT){
    x1[t] = x1[t-1]  + pro_dev1[t];
    x2[t] = x2[t-1]  + pro_dev2[t];
  }
   
  for(i in 1:N_yr1){
    a_yr1[i] = x1[yr_index1[i]] + obs_dev1[i]; 
  }
  for(i in 1:N_yr2){
    a_yr2[i] = x2[yr_index2[i]] + obs_dev2[i]; 
  }
  
}  

model{
  //priors
  c ~ induced_dirichlet(rep_vector(1, K), 0); //prior on ordered cut-points
  beta1 ~ normal(0,2); //covariates - rvc
  beta2 ~ normal(0,2); //covariates - reef
  x01 ~ normal(0,5); //initial state - rvc
  x02 ~ normal(0,5); //initial state - reef

  //variance terms
  sd_hab1 ~ cauchy(0, 1);
  sd_hab2 ~ cauchy(0, 1);
  sd_q1 ~ inv_gamma(2,0.25);
  sd_q2 ~ inv_gamma(2,0.25);
  sd_r1 ~ inv_gamma(2,0.25);
  sd_r2 ~ inv_gamma(2,0.25);
  sd_site ~ cauchy(0, 1);
  sd_dv ~ cauchy(0, 1);
  sd_dmy ~ cauchy(0, 1);
  
  //varying intercepts
  a_hab1 ~ normal(0, sd_hab1);
  a_hab2 ~ normal(0, sd_hab2);
  a_site ~ normal(0, sd_site);
  a_dv ~ normal(0, sd_dv);
  a_dmy ~ normal(0, sd_dmy);
  
  for(t in 1:N_yr1){
   obs_dev1[t] ~ normal(0,sd_r1);
  }
  
  for(t in 1:TT){
    pro_dev1[t] ~ normal(0, sd_q1);
    pro_dev2[t] ~ normal(0, sd_q2);
    obs_dev2[t] ~ normal(0,sd_r2);
  }
  
  y1 ~ neg_binomial_2_log(a_yr1[year_id1] + a_hab1[hab_class1] + X1*beta1,phi);
  
  y2 ~ ordered_logistic(a_hab2[hab_class2]+a_yr2[year_id2]+a_site[site]+a_dv[diver]+a_dmy[dmy]+X2*beta2,c);
  
}
 generated quantities{
  vector[N1+N2] log_lik;
  for (i in 1:N1) log_lik[i] = neg_binomial_2_log_lpmf(y1[i]|a_yr1[year_id1[i]] + a_hab1[hab_class1[i]] + X1[i,]*beta1,phi);
  for (z in 1:N2) log_lik[N1+z] = ordered_logistic_lpmf(y2[z]|a_hab2[hab_class2[z]]+a_yr2[year_id2[z]]+a_site[site[z]]+a_dv[diver[z]]+a_dmy[dmy[z]]+X2[z,]*beta2, c);
}  
"


abund_test_SS_sep_trend<-"functions {
  real induced_dirichlet_lpdf(vector c, vector alpha, real phi) {
    int K = num_elements(c) + 1;
    vector[K - 1] sigma = inv_logit(phi - c);
    vector[K] p;
    matrix[K, K] J = rep_matrix(0, K, K);
    
    // Induced ordinal probabilities
    p[1] = 1 - sigma[1];
    for (k in 2:(K - 1))
      p[k] = sigma[k - 1] - sigma[k];
    p[K] = sigma[K - 1];
    
    // Baseline column of Jacobian
    for (k in 1:K) J[k, 1] = 1;
    
    // Diagonal entries of Jacobian
    for (k in 2:K) {
      real rho = sigma[k - 1] * (1 - sigma[k - 1]);
      J[k, k] = - rho;
      J[k - 1, k] = rho;
    }
    
    return   dirichlet_lpdf(p | alpha)
           + log_determinant(J);
  }
}
data{
  int<lower=1> N1;//number of observations (REEF surveys)
  int<lower=1> N2;//number of observations (REEF surveys)
  int y1[N1]; //abundance category for each survey
  int y2[N2]; //abundance category for each survey
  int<lower=0> N_hab1; //number of habitat classes
  int<lower=1,upper=N_hab1> hab_class1[N1]; // vector of habitat class identities
  int<lower=0> N_hab2; //number of habitat classes
  int<lower=1,upper=N_hab2> hab_class2[N2]; // vector of habitat class identities
  int<lower=0> N_site; //number of sites
  int<lower=1,upper=N_site> site[N2]; // vector of site identities
  int<lower=0> N_dv; //number of divers
  int<lower=1,upper=N_dv> diver[N2]; // vector of diver identities
  int<lower=0> N_dmy; //number of site day clusters
  int<lower=1,upper=N_dmy> dmy[N2]; // vector of site day cluster identities
  int Z1; // columns in the covariate matrix
  int Z2; // columns in the covariate matrix
  matrix[N1,Z1] X1; // design matrix X
  matrix[N2,Z2] X2; // design matrix X
  int K; //ordinal levels
  int TT; // timespan
  int<lower=0> N_yr1; //number of years
  int yr_index1[N_yr1]; //index of years
  int<lower=1,upper=N_yr1> year_id1[N1]; // vector of year
  int<lower=0> N_yr2; //number of years
  int yr_index2[N_yr2]; //index of years
  int<lower=1,upper=N_yr2> year_id2[N2]; // vector of year

}
parameters {
  ordered[K-1] c; //cutpoints
  real x01; //initial popn size - rvc
  real x02; //initial popn size - reef
  real<lower = 0> recip_phi; //overdispersion parameter
  real u1; //trend constant - rvc
  real u2; //trend constant - reef
  
  //deviations from intercept
  vector[Z1] beta1; //effort coefficients - RVC
  vector[Z2] beta2; //effort coefficients - REEF
  vector[N_hab1] a_hab1; //deviation between habitats - RVC
  vector[N_hab2] a_hab2; //deviation between habitats - REEF
  vector[N_site] a_site; //deviation between sites
  vector[N_dv] a_dv; //deviation between divers
  vector[N_dmy] a_dmy; //deviation between site day clusters
 
  //variance on the deviance components
  real<lower = 0> sd_hab1;
  real<lower = 0> sd_hab2;
  real<lower = 0> sd_site;
  real<lower = 0> sd_dv;
  real<lower = 0> sd_dmy;
  real<lower = 0> sd_r1;
  real<lower = 0> sd_r2;
  real<lower = 0> sd_q1;
  real<lower = 0> sd_q2;
  
  vector[TT] pro_dev1; //process deviations
  vector[TT] pro_dev2; //process deviations
  vector[N_yr1] obs_dev1; //observation deviations 
  vector[N_yr2] obs_dev2; //observation deviations 
  
}

transformed parameters{
  vector[TT] x1;
  vector[TT] x2;
  vector[N_yr1] a_yr1;
  vector[N_yr2] a_yr2;
  real phi;
 
  phi=1/recip_phi;
  
  x1[1] = x01 + pro_dev1[1];
  x2[1] = x02 + pro_dev2[1];
   
  for(t in 2:TT){
    x1[t] = x1[t-1] + u1 + pro_dev1[t];
    x2[t] = x2[t-1] + u2 + pro_dev2[t];
  }
   
  for(i in 1:N_yr1){
    a_yr1[i] = x1[yr_index1[i]] + obs_dev1[i]; 
  }
  for(i in 1:N_yr2){
    a_yr2[i] = x2[yr_index2[i]] + obs_dev2[i]; 
  }
  
}  

model{
  //priors
  c ~ induced_dirichlet(rep_vector(1, K), 0); //prior on ordered cut-points
  beta1 ~ normal(0,2); //covariates - rvc
  beta2 ~ normal(0,2); //covariates - reef
  x01 ~ normal(0,5); //initial state - rvc
  x02 ~ normal(0,5); //initial state - reef
  u1 ~ normal(0,0.5); //trend constant - rvc
  u2 ~ normal(0,0.5); //trend constant - reef
 
  //variance terms
  sd_hab1 ~ cauchy(0, 1);
  sd_hab2 ~ cauchy(0, 1);
  sd_q1 ~ cauchy(0, 1);
  sd_q2 ~ cauchy(0, 1);
  sd_r1 ~ cauchy(0, 1);
  sd_r2 ~ cauchy(0, 1);
  sd_site ~ cauchy(0, 1);
  sd_dv ~ cauchy(0, 1);
  sd_dmy ~ cauchy(0, 1);
  
  //varying intercepts
  a_hab1 ~ normal(0, sd_hab1);
  a_hab2 ~ normal(0, sd_hab2);
  a_site ~ normal(0, sd_site);
  a_dv ~ normal(0, sd_dv);
  a_dmy ~ normal(0, sd_dmy);
  
  for(t in 1:N_yr1){
   obs_dev1[t] ~ normal(0,sd_r1);
  }
  
  for(t in 1:TT){
    pro_dev1[t] ~ normal(0, sd_q1);
    pro_dev2[t] ~ normal(0, sd_q2);
    obs_dev2[t] ~ normal(0,sd_r2);
  }
  
  y1 ~ neg_binomial_2_log(a_yr1[year_id1] + a_hab1[hab_class1] + X1*beta1,phi);
  
  y2 ~ ordered_logistic(a_hab2[hab_class2]+a_yr2[year_id2]+a_site[site]+a_dv[diver]+a_dmy[dmy]+X2*beta2,c);
  
}
 generated quantities{
  vector[N1+N2] log_lik;
  for (i in 1:N1) log_lik[i] = neg_binomial_2_log_lpmf(y1[i]|a_yr1[year_id1[i]] + a_hab1[hab_class1[i]] + X1[i,]*beta1,phi);
  for (z in 1:N2) log_lik[N1+z] = ordered_logistic_lpmf(y2[z]|a_hab2[hab_class2[z]]+a_yr2[year_id2[z]]+a_site[site[z]]+a_dv[diver[z]]+a_dmy[dmy[z]]+X2[z,]*beta2, c);
}  
"

occ_test_SS_comb<-"
data{
  int<lower=1> N1;//number of observations - rvc
  int<lower=1> N2;//number of observations - reef
  int y1[N1]; //occurence for each survey - rvc
  int y2[N2]; //occurence category for each survey
  int<lower=0> N_hab1; //number of habitat classes
  int<lower=1,upper=N_hab1> hab_class1[N1]; // vector of habitat class identities
  int<lower=0> N_hab2; //number of habitat classes
  int<lower=1,upper=N_hab2> hab_class2[N2]; // vector of habitat class identities
  int<lower=0> N_site; //number of sites
  int<lower=1,upper=N_site> site[N2]; // vector of site identities
  int<lower=0> N_dv; //number of divers
  int<lower=1,upper=N_dv> diver[N2]; // vector of diver identities
  int<lower=0> N_dmy; //number of site day clusters
  int<lower=1,upper=N_dmy> dmy[N2]; // vector of site day cluster identities
  int Z1; // columns in the covariate matrix
  int Z2; // columns in the covariate matrix
  matrix[N1,Z1] X1; // design matrix X
  matrix[N2,Z2] X2; // design matrix X
  int TT; // timespan
  int<lower=0> N_yr1; //number of years
  int yr_index1[N_yr1]; //index of years
  int<lower=1,upper=N_yr1> year_id1[N1]; // vector of year
  int<lower=0> N_yr2; //number of years
  int yr_index2[N_yr2]; //index of years
  int<lower=1,upper=N_yr2> year_id2[N2]; // vector of year

}
parameters {
  real x0; //initial popn size
  real a; //scalar

  //deviations from intercept
  vector[Z1] beta1; //effort coefficients - RVC
  vector[Z2] beta2; //effort coefficients - REEF
  vector[N_hab1] a_hab1; //deviation between habitats - RVC
  vector[N_hab2] a_hab2; //deviation between habitats - REEF
  vector[N_site] a_site; //deviation between sites
  vector[N_dv] a_dv; //deviation between divers
  vector[N_dmy] a_dmy; //deviation between site day clusters
 
  //variance on the deviance components
  real<lower = 0> sd_hab1;
  real<lower = 0> sd_hab2;
  real<lower = 0> sd_site;
  real<lower = 0> sd_dv;
  real<lower = 0> sd_dmy;
  real<lower = 0> sd_r1;
  real<lower = 0> sd_r2;
  real<lower = 0> sd_q;

  
  vector[TT] pro_dev; //process deviations
  vector[N_yr1] obs_dev1; //observation deviations 
  vector[N_yr2] obs_dev2; //observation deviations 
  
}

transformed parameters{
  vector[TT] x;
  vector[N_yr1] a_yr1;
  vector[N_yr2] a_yr2;

  x[1] = x0 + pro_dev[1];
   
  for(t in 2:TT){
    x[t] = x[t-1] + pro_dev[t];
  }
   
  for(i in 1:N_yr1){
    a_yr1[i] = x[yr_index1[i]] + obs_dev1[i]; 
  }
  for(i in 1:N_yr2){
    a_yr2[i] = x[yr_index2[i]] + a + obs_dev2[i]; 
  }
  
}  

model{
  //priors
  beta1 ~ normal(0,2); //covariates - rvc
  beta2 ~ normal(0,2); //covariates - reef
  x0 ~ normal(0,5); //initial state
  a ~ normal(0,3); //scalar

  //variance terms
  sd_hab1 ~ cauchy(0, 1);
  sd_hab2 ~ cauchy(0, 1);
  sd_q ~ cauchy(0, 1);
  sd_r1 ~ cauchy(0, 1);
  sd_r2 ~ cauchy(0, 1);
  sd_site ~ cauchy(0, 1);
  sd_dv ~ cauchy(0, 1);
  sd_dmy ~ cauchy(0, 1);
  
  //varying intercepts
  a_hab1 ~ normal(0, sd_hab1);
  a_hab2 ~ normal(0, sd_hab2);
  a_site ~ normal(0, sd_site);
  a_dv ~ normal(0, sd_dv);
  a_dmy ~ normal(0, sd_dmy);
  
  for(t in 1:N_yr1){
   obs_dev1[t] ~ normal(0,sd_r1);
  }
  
  for(t in 1:TT){
    pro_dev[t] ~ normal(0, sd_q);
    obs_dev2[t] ~ normal(0,sd_r2);
  }
  
  y1 ~ bernoulli_logit(a_yr1[year_id1] + a_hab1[hab_class1] + X1*beta1);
  
  y2 ~ bernoulli_logit(a_yr2[year_id2]+ a_hab2[hab_class2]+a_site[site]+a_dv[diver]+a_dmy[dmy]+X2*beta2);
  
}
 generated quantities{
  vector[N1+N2] log_lik;
  for (i in 1:N1) log_lik[i] = bernoulli_logit_lpmf(y1[i]|a_yr1[year_id1[i]] + a_hab1[hab_class1[i]] + X1[i,]*beta1);
  for (z in 1:N2) log_lik[N1+z] = bernoulli_logit_lpmf(y2[z]|a_hab2[hab_class2[z]]+a_yr2[year_id2[z]]+a_site[site[z]]+a_dv[diver[z]]+a_dmy[dmy[z]]+X2[z,]*beta2);
}  
"

occ_test_SS_sep<-"
data{
  int<lower=1> N1;//number of observations - rvc
  int<lower=1> N2;//number of observations - reef
  int y1[N1]; //occurence for each survey - rvc
  int y2[N2]; //occurence category for each survey
  int<lower=0> N_hab1; //number of habitat classes
  int<lower=1,upper=N_hab1> hab_class1[N1]; // vector of habitat class identities
  int<lower=0> N_hab2; //number of habitat classes
  int<lower=1,upper=N_hab2> hab_class2[N2]; // vector of habitat class identities
  int<lower=0> N_site; //number of sites
  int<lower=1,upper=N_site> site[N2]; // vector of site identities
  int<lower=0> N_dv; //number of divers
  int<lower=1,upper=N_dv> diver[N2]; // vector of diver identities
  int<lower=0> N_dmy; //number of site day clusters
  int<lower=1,upper=N_dmy> dmy[N2]; // vector of site day cluster identities
  int Z1; // columns in the covariate matrix
  int Z2; // columns in the covariate matrix
  matrix[N1,Z1] X1; // design matrix X
  matrix[N2,Z2] X2; // design matrix X
  int TT; // timespan
  int<lower=0> N_yr1; //number of years
  int yr_index1[N_yr1]; //index of years
  int<lower=1,upper=N_yr1> year_id1[N1]; // vector of year
  int<lower=0> N_yr2; //number of years
  int yr_index2[N_yr2]; //index of years
  int<lower=1,upper=N_yr2> year_id2[N2]; // vector of year

}
parameters {
  real x01; //initial popn size - rvc
  real x02; //initial popn size - reef

  //deviations from intercept
  vector[Z1] beta1; //effort coefficients - RVC
  vector[Z2] beta2; //effort coefficients - REEF
  vector[N_hab1] a_hab1; //deviation between habitats - RVC
  vector[N_hab2] a_hab2; //deviation between habitats - REEF
  vector[N_site] a_site; //deviation between sites
  vector[N_dv] a_dv; //deviation between divers
  vector[N_dmy] a_dmy; //deviation between site day clusters
 
  //variance on the deviance components
  real<lower = 0> sd_hab1;
  real<lower = 0> sd_hab2;
  real<lower = 0> sd_site;
  real<lower = 0> sd_dv;
  real<lower = 0> sd_dmy;
  real<lower = 0> sd_r1;
  real<lower = 0> sd_r2;
  real<lower = 0> sd_q1;
  real<lower = 0> sd_q2;
  
  vector[TT] pro_dev1; //process deviations
  vector[TT] pro_dev2; //process deviations
  vector[N_yr1] obs_dev1; //observation deviations 
  vector[N_yr2] obs_dev2; //observation deviations 
  
}

transformed parameters{
  vector[TT] x1;
  vector[TT] x2;
  vector[N_yr1] a_yr1;
  vector[N_yr2] a_yr2;

  x1[1] = x01 + pro_dev1[1];
  x2[1] = x02 + pro_dev2[1];
   
  for(t in 2:TT){
    x1[t] = x1[t-1] + pro_dev1[t];
    x2[t] = x2[t-1] + pro_dev2[t];
  }
   
  for(i in 1:N_yr1){
    a_yr1[i] = x1[yr_index1[i]] + obs_dev1[i]; 
  }
  for(i in 1:N_yr2){
    a_yr2[i] = x2[yr_index2[i]] + obs_dev2[i]; 
  }
  
}  

model{
  //priors
  beta1 ~ normal(0,2); //covariates - rvc
  beta2 ~ normal(0,2); //covariates - reef
  x01 ~ normal(0,5); //initial state - rvc
  x02 ~ normal(0,5); //initial state - reef

  //variance terms
  sd_hab1 ~ cauchy(0, 1);
  sd_hab2 ~ cauchy(0, 1);
  sd_q1 ~ cauchy(0, 1);
  sd_q2 ~ cauchy(0, 1);
  sd_r1 ~ cauchy(0, 1);
  sd_r2 ~ cauchy(0, 1);
  sd_site ~ cauchy(0, 1);
  sd_dv ~ cauchy(0, 1);
  sd_dmy ~ cauchy(0, 1);
  
  //varying intercepts
  a_hab1 ~ normal(0, sd_hab1);
  a_hab2 ~ normal(0, sd_hab2);
  a_site ~ normal(0, sd_site);
  a_dv ~ normal(0, sd_dv);
  a_dmy ~ normal(0, sd_dmy);
  
  for(t in 1:N_yr1){
   obs_dev1[t] ~ normal(0,sd_r1);
  }
  
  for(t in 1:TT){
    pro_dev1[t] ~ normal(0, sd_q1);
    pro_dev2[t] ~ normal(0, sd_q2);
    obs_dev2[t] ~ normal(0,sd_r2);
  }
  
  y1 ~ bernoulli_logit(a_yr1[year_id1] + a_hab1[hab_class1] + X1*beta1);
  
  y2 ~ bernoulli_logit(a_yr2[year_id2]+ a_hab2[hab_class2]+a_site[site]+a_dv[diver]+a_dmy[dmy]+X2*beta2);
  
}
 generated quantities{
  vector[N1+N2] log_lik;
  for (i in 1:N1) log_lik[i] = bernoulli_logit_lpmf(y1[i]|a_yr1[year_id1[i]] + a_hab1[hab_class1[i]] + X1[i,]*beta1);
  for (z in 1:N2) log_lik[N1+z] = bernoulli_logit_lpmf(y2[z]|a_hab2[hab_class2[z]]+a_yr2[year_id2[z]]+a_site[site[z]]+a_dv[diver[z]]+a_dmy[dmy[z]]+X2[z,]*beta2);
}  
"


####Combined state - blue angelfish

blue_angel<- rvc_occs[[3]]

year_index<- data.frame(yr=seq(1993,2018),y.ind=seq(1,26))
blue_angel$year_index=year_index$y.ind[match(blue_angel$YEAR,year_index$yr)]
X_rvc<- matrix(data=c(scale(as.numeric(blue_angel$DEPTH))),ncol=1,nrow=nrow(blue_angel))


test_nb_rvc<- rstan::stan(model_code = nb_test_rvc, data = list(y = blue_angel$NUM.total2, 
                                                                 N = nrow(blue_angel),
                                                                 N_hab = length(unique(blue_angel$HAB_CD2)),
                                                                 hab_class=as.numeric(factor(blue_angel$HAB_CD2)),
                                                                 N_yr =length(unique(blue_angel$YEAR)),
                                                                 year=as.numeric(factor(blue_angel$YEAR)),
                                                                 X=X_rvc,
                                                                 K=ncol(X_rvc)),
                        pars = c("alpha", "tau_hab",'z_hab','tau_yr','z_yr','phi','beta'),
                        control = list(adapt_delta = 0.995,max_treedepth = 15), warmup = 200, chains = 4, iter = 400, thin = 1)

shinystan::launch_shinystan(test_nb_rvc)

posterior <- as.array(test_nb_rvc)
dimnames(posterior)
color_scheme_set("viridis")
mcmc_areas(
  posterior,
  pars = c('alpha',paste('z_hab[',seq(1:4),']',sep='')),
  prob = 0.95, # 80% intervals
  prob_outer = 1, # 99%
  point_est = "median"
)
mcmc_dens_overlay(posterior, c("alpha", paste('a_yr[',seq(1:23),']',sep=''),'sigma_hab','sigma_yr'))
mcmc_pairs(posterior, pars = c('alpha',paste('z_hab[',seq(1:4),']',sep=''),'tau_hab'),
           off_diag_args = list(size = 1.5))
mcmc_trace(posterior, pars = c("tau_hab"))
mcmc_trace(posterior, pars = c("alpha"))
mcmc_trace(posterior, pars = c("sd_r"))


ba_params<- rstan::extract(test_1_comb)
x_mat<- data.frame(median=NA,l.95=NA,u.95=NA)
for(i in 1:26){
  x_mat[i,1]=median(plogis(ba_params$x[,i]))
  x_mat[i,2]=quantile(plogis(ba_params$x[,i]),0.975)
  x_mat[i,3]=quantile(plogis(ba_params$x[,i]),0.025)
}
y_mat_1<- data.frame(median=NA,l.95=NA,u.95=NA)
for(i in 1:23){
  y_mat_1[i,1]=median(plogis(ba_params$a_yr1[,i]))
  y_mat_1[i,2]=quantile(plogis(ba_params$a_yr1[,i]),0.975)
  y_mat_1[i,3]=quantile(plogis(ba_params$a_yr1[,i]),0.025)
}

y_mat_2<- data.frame(median=NA,l.95=NA,u.95=NA)
for(i in 1:26){
  y_mat_2[i,1]=median(plogis(ba_params$a_yr2[,i]))
  y_mat_2[i,2]=quantile(plogis(ba_params$a_yr2[,i]+ba_params$s),0.975)
  y_mat_2[i,3]=quantile(plogis(ba_params$a_yr2[,i]+ba_params$s),0.025)
}

comp_plot_SS_comb(ts1=rvc_ts[[1]],ts2=reef_ts[[1]],x_mat = x_mat,y_mat1=y_mat_1,y_mat2=y_mat_2,sp='Blue Angelfish',GZ='Key Largo')


blue_angel<- rvc_occs[[3]]

year_index<- data.frame(yr=seq(1993,2018),y.ind=seq(1,26))
blue_angel$year_index=year_index$y.ind[match(blue_angel$YEAR,year_index$yr)]
X_rvc<- matrix(data=c(scale(as.numeric(blue_angel$DEPTH))),ncol=1,nrow=nrow(blue_angel))

test_nb_rvc<- rstan::stan(model_code = nb_test_SS_rvc, data = list(y = blue_angel$NUM.total2, 
                                                                N = nrow(blue_angel),
                                                                N_hab = length(unique(blue_angel$HAB_CD2)),
                                                                hab_class=as.numeric(factor(blue_angel$HAB_CD2)),
                                                               
                                                                N_yr =length(unique(blue_angel$YEAR)),
                                                                year_id=as.numeric(factor(blue_angel$YEAR)),
                                                                X=X_rvc,
                                                                K=ncol(X_rvc),
                                                                TT=26,
                                                                yr_index=sort(unique(blue_angel$year_index))),
                          pars = c('x',"sd_hab",'a_hab','sd_r','sd_q','a_yr','phi','beta'),
                          control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 200, chains = 4, iter = 600, thin = 1)

shinystan::launch_shinystan(test_nb_rvc)


##Reef ordinal model

library(brms)
test<- brm(ordered(abundance)~1+(1|geogr)+(1|hab_class2),
               data=reef_occs[[3]],
               family=cumulative("logit"),
               prior = c(prior(normal(0,10), "Intercept")),
               control = list(adapt_delta = 0.98),
               cores=6,
               save_all_pars = T,
               warmup = 100, iter = 3000, chains = 2)

brms::make_stancode(ordered(abundance)~1+(1|geogr)+(1|hab_class2),
              data=reef_occs[[3]],
              family=cumulative("logit"))

reef_occs[[3]]<- reef_occs[[3]] %>% mutate(exp_binary=ifelse(exp=='E',1,0))

X_rvc<- matrix(data=c(scale(as.numeric(blue_angel$DEPTH))),ncol=1,nrow=nrow(blue_angel))
X<- matrix(data=c(scale(as.numeric(reef_occs[[3]]$btime)),scale(as.numeric(reef_occs[[3]]$averagedepth)),scale(as.numeric(reef_occs[[3]]$visibility)),scale(as.numeric(reef_occs[[3]]$current)),reef_occs[[3]]$exp_binary),ncol=5,nrow=nrow(reef_occs[[3]]))

reef_occs[[3]]$abundance2<- reef_occs[[3]]$abundance+1
test_ord_reef<- rstan::stan(model_code = ord_test_reef, data = list(y =reef_occs[[3]]$abundance2, 
                                                                    N = nrow(reef_occs[[3]]),
                                                                    N_hab = length(unique(reef_occs[[3]]$hab_class2)),
                                                                    hab_class=as.numeric(factor(reef_occs[[3]]$hab_class2)),
                                                                    N_yr =length(unique(reef_occs[[3]]$year)),
                                                                    year=as.numeric(factor(reef_occs[[3]]$year)),
                                                                    site=as.numeric(factor(reef_occs[[3]]$geogr)),
                                                                    N_site=length(unique(reef_occs[[3]]$geogr)),
                                                                    diver=as.numeric(factor(reef_occs[[3]]$fish_memberid)),
                                                                    N_dv=length(unique(reef_occs[[3]]$fish_memberid)),
                                                                    Z=ncol(X),
                                                                    X=X,
                                                                    K=length(unique(reef_occs[[3]]$abundance))),
                          pars = c('a_yr',"sigma_hab",'sigma_yr','sigma_site','sigma_dv','c'),
                          control = list(adapt_delta = 0.995,max_treedepth = 15), warmup = 200, chains = 4, iter = 500, thin = 1)

test_ord_reef<- rstan::stan(model_code = ord_test_reef1, data = list(y =reef_occs[[3]]$abundance2, 
                                                                    N = nrow(reef_occs[[3]]),
                                                                    K=length(unique(reef_occs[[3]]$abundance))),
                            pars = c('c','gamma'),
                            control = list(adapt_delta = 0.995,max_treedepth = 15), warmup = 200, chains = 4, iter = 500, thin = 1)

test_ord_reef<- rstan::stan(model_code = ord_test_reef, data = list(y =reef_occs[[3]]$abundance2, 
                                                                     N = nrow(reef_occs[[3]]),
                                                                     N_hab = length(unique(reef_occs[[3]]$hab_class2)),
                                                                     hab_class=as.numeric(factor(reef_occs[[3]]$hab_class2)),
                                                                     N_yr =length(unique(reef_occs[[3]]$year)),
                                                                     year=as.numeric(factor(reef_occs[[3]]$year)),
                                                                     site=as.numeric(factor(reef_occs[[3]]$geogr)),
                                                                     N_site=length(unique(reef_occs[[3]]$geogr)),
                                                                     diver=as.numeric(factor(reef_occs[[3]]$fish_memberid)),
                                                                     N_dv=length(unique(reef_occs[[3]]$fish_memberid)),
                                                                     dmy=as.numeric(factor(reef_occs[[3]]$site_dmy)),
                                                                     N_dmy=length(unique(reef_occs[[3]]$site_dmy)),
                                                                     K=length(unique(reef_occs[[3]]$abundance)),
                                                                     X=X,
                                                                    Z=ncol(X)),
                            pars = c('c','a_hab','sigma_hab','sigma_yr','sigma_site','sigma_dv','sigma_dmy'),
                            control = list(adapt_delta = 0.995,max_treedepth = 15), warmup = 200, chains = 4, iter = 500, thin = 1)

test_ord_reef2<- rstan::stan(model_code = ord_test_SS_reef, data = list(y =reef_occs[[3]]$abundance2, 
                                                                     N = nrow(reef_occs[[3]]),
                                                                     N_hab = length(unique(reef_occs[[3]]$hab_class2)),
                                                                     hab_class=as.numeric(factor(reef_occs[[3]]$hab_class2)),
                                                                     site=as.numeric(factor(reef_occs[[3]]$geogr)),
                                                                     N_site=length(unique(reef_occs[[3]]$geogr)),
                                                                     diver=as.numeric(factor(reef_occs[[3]]$fish_memberid)),
                                                                     N_dv=length(unique(reef_occs[[3]]$fish_memberid)),
                                                                     dmy=as.numeric(factor(reef_occs[[3]]$site_dmy)),
                                                                     N_dmy=length(unique(reef_occs[[3]]$site_dmy)),
                                                                     K=length(unique(reef_occs[[3]]$abundance)),
                                                                     X=X,
                                                                     Z=ncol(X),
                                                                     TT=26,
                                                                     N_yr=length(unique(reef_occs[[3]]$year)),
                                                                     yr_index=sort(unique(as.numeric(factor(reef_occs[[3]]$year)))),
                                                                     year_id=as.numeric(factor(reef_occs[[3]]$year))),
                             pars = c('c','a_hab','sd_hab','sd_site','sd_dv','sd_dmy','sd_r','sd_q','beta','x'),
                             control = list(adapt_delta = 0.995,max_treedepth = 15), warmup = 200, chains = 4, iter = 500, thin = 1)

X1<- matrix(data=c(scale(as.numeric(rvc_occs[[3]]$DEPTH))),ncol=1,nrow=nrow(rvc_occs[[3]]))
X2<- matrix(data=c(scale(as.numeric(reef_occs[[3]]$btime)),scale(as.numeric(reef_occs[[3]]$averagedepth)),scale(as.numeric(reef_occs[[3]]$visibility)),scale(as.numeric(reef_occs[[3]]$current)),reef_occs[[3]]$exp_binary),ncol=5,nrow=nrow(reef_occs[[3]]))


test_comb<- rstan::stan(model_code = abund_test_SS_comb, data = list(y1 = rvc_occs[[3]]$NUM.total2,
                                                                        y2 =reef_occs[[3]]$abundance2,
                                                                        N1 = nrow(rvc_occs[[3]]),
                                                                        N2 = nrow(reef_occs[[3]]),
                                                                        N_hab1 = length(unique(rvc_occs[[3]]$HAB_CD2)),
                                                                        hab_class1=as.numeric(factor(rvc_occs[[3]]$HAB_CD2)),
                                                                        N_hab2 = length(unique(reef_occs[[3]]$hab_class2)),
                                                                        hab_class2=as.numeric(factor(reef_occs[[3]]$hab_class2)),
                                                                        site=as.numeric(factor(reef_occs[[3]]$geogr)),
                                                                        N_site=length(unique(reef_occs[[3]]$geogr)),
                                                                        diver=as.numeric(factor(reef_occs[[3]]$fish_memberid)),
                                                                        N_dv=length(unique(reef_occs[[3]]$fish_memberid)),
                                                                        dmy=as.numeric(factor(reef_occs[[3]]$site_dmy)),
                                                                        N_dmy=length(unique(reef_occs[[3]]$site_dmy)),
                                                                        K=length(unique(reef_occs[[3]]$abundance)),
                                                                        X1=X1,
                                                                        Z1=ncol(X1),
                                                                        X2=X2,
                                                                        Z2=ncol(X2),
                                                                        TT=26,
                                                                        N_yr1=length(unique(rvc_occs[[3]]$YEAR)),
                                                                        yr_index1=sort(unique(as.numeric(factor(rvc_occs[[3]]$YEAR)))),
                                                                        year_id1=as.numeric(factor(rvc_occs[[3]]$YEAR)),
                                                                        N_yr2=length(unique(reef_occs[[3]]$year)),
                                                                        yr_index2=sort(unique(as.numeric(factor(reef_occs[[3]]$year)))),
                                                                        year_id2=as.numeric(factor(reef_occs[[3]]$year))),
                             pars = c('c','a_hab1','a_hab2','sd_hab1','sd_hab2','sd_site','sd_dv','sd_dmy','sd_r1','sd_r2','sd_q','x','a','log_lik'),
                             control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 150, chains = 4, iter = 450, thin = 1)

test_sep<- rstan::stan(model_code = abund_test_SS_sep, data = list(y1 = rvc_occs[[3]]$NUM.total2,
                                                                     y2 =reef_occs[[3]]$abundance2,
                                                                     N1 = nrow(rvc_occs[[3]]),
                                                                     N2 = nrow(reef_occs[[3]]),
                                                                     N_hab1 = length(unique(rvc_occs[[3]]$HAB_CD2)),
                                                                     hab_class1=as.numeric(factor(rvc_occs[[3]]$HAB_CD2)),
                                                                     N_hab2 = length(unique(reef_occs[[3]]$hab_class2)),
                                                                     hab_class2=as.numeric(factor(reef_occs[[3]]$hab_class2)),
                                                                     site=as.numeric(factor(reef_occs[[3]]$geogr)),
                                                                     N_site=length(unique(reef_occs[[3]]$geogr)),
                                                                     diver=as.numeric(factor(reef_occs[[3]]$fish_memberid)),
                                                                     N_dv=length(unique(reef_occs[[3]]$fish_memberid)),
                                                                     dmy=as.numeric(factor(reef_occs[[3]]$site_dmy)),
                                                                     N_dmy=length(unique(reef_occs[[3]]$site_dmy)),
                                                                     K=length(unique(reef_occs[[3]]$abundance)),
                                                                     X1=X1,
                                                                     Z1=ncol(X1),
                                                                     X2=X2,
                                                                     Z2=ncol(X2),
                                                                     TT=26,
                                                                     N_yr1=length(unique(rvc_occs[[3]]$YEAR)),
                                                                     yr_index1=sort(unique(as.numeric(factor(rvc_occs[[3]]$YEAR)))),
                                                                     year_id1=as.numeric(factor(rvc_occs[[3]]$YEAR)),
                                                                     N_yr2=length(unique(reef_occs[[3]]$year)),
                                                                     yr_index2=sort(unique(as.numeric(factor(reef_occs[[3]]$year)))),
                                                                     year_id2=as.numeric(factor(reef_occs[[3]]$year))),
                        pars = c('c','a_hab1','a_hab2','sd_hab1','sd_hab2','sd_site','sd_dv','sd_dmy','sd_r1','sd_r2','sd_q1','sd_q2','x1','x2','a_yr1','a_yr2','log_lik'),
                        control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 150, chains = 4, iter = 450, thin = 1)
system.time(loo1= loo::loo(test_comb))
start.time <- Sys.time()
system.time(loo2<- loo::loo(test_sep_rf))


shinystan::launch_shinystan(test_sep)

params<- rstan::extract(test_sep)

test_occ_sep<- rstan::stan(model_code = occ_test_SS_sep, data = list(y1 = rvc_occs[[3]]$occ,
                                                                   y2 =reef_occs[[3]]$occ,
                                                                   N1 = nrow(rvc_occs[[3]]),
                                                                   N2 = nrow(reef_occs[[3]]),
                                                                   N_hab1 = length(unique(rvc_occs[[3]]$HAB_CD2)),
                                                                   hab_class1=as.numeric(factor(rvc_occs[[3]]$HAB_CD2)),
                                                                   N_hab2 = length(unique(reef_occs[[3]]$hab_class2)),
                                                                   hab_class2=as.numeric(factor(reef_occs[[3]]$hab_class2)),
                                                                   site=as.numeric(factor(reef_occs[[3]]$geogr)),
                                                                   N_site=length(unique(reef_occs[[3]]$geogr)),
                                                                   diver=as.numeric(factor(reef_occs[[3]]$fish_memberid)),
                                                                   N_dv=length(unique(reef_occs[[3]]$fish_memberid)),
                                                                   dmy=as.numeric(factor(reef_occs[[3]]$site_dmy)),
                                                                   N_dmy=length(unique(reef_occs[[3]]$site_dmy)),
                                                                   X1=X1,
                                                                   Z1=ncol(X1),
                                                                   X2=X2,
                                                                   Z2=ncol(X2),
                                                                   TT=26,
                                                                   N_yr1=length(unique(rvc_occs[[3]]$YEAR)),
                                                                   yr_index1=sort(unique(as.numeric(factor(rvc_occs[[3]]$YEAR)))),
                                                                   year_id1=as.numeric(factor(rvc_occs[[3]]$YEAR)),
                                                                   N_yr2=length(unique(reef_occs[[3]]$year)),
                                                                   yr_index2=sort(unique(as.numeric(factor(reef_occs[[3]]$year)))),
                                                                   year_id2=as.numeric(factor(reef_occs[[3]]$year))),
                       pars = c('a_hab1','a_hab2','sd_hab1','sd_hab2','sd_site','sd_dv','sd_dmy','sd_r1','sd_r2','sd_q1','sd_q2','x1','x2','a_yr1','a_yr2',
                                'log_lik'),
                       control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 150, chains = 4, iter = 450, thin = 1)

test_occ_comb<- rstan::stan(model_code = occ_test_SS_comb, data = list(y1 = rvc_occs[[3]]$occ,
                                                                     y2 =reef_occs[[3]]$occ,
                                                                     N1 = nrow(rvc_occs[[3]]),
                                                                     N2 = nrow(reef_occs[[3]]),
                                                                     N_hab1 = length(unique(rvc_occs[[3]]$HAB_CD2)),
                                                                     hab_class1=as.numeric(factor(rvc_occs[[3]]$HAB_CD2)),
                                                                     N_hab2 = length(unique(reef_occs[[3]]$hab_class2)),
                                                                     hab_class2=as.numeric(factor(reef_occs[[3]]$hab_class2)),
                                                                     site=as.numeric(factor(reef_occs[[3]]$geogr)),
                                                                     N_site=length(unique(reef_occs[[3]]$geogr)),
                                                                     diver=as.numeric(factor(reef_occs[[3]]$fish_memberid)),
                                                                     N_dv=length(unique(reef_occs[[3]]$fish_memberid)),
                                                                     dmy=as.numeric(factor(reef_occs[[3]]$site_dmy)),
                                                                     N_dmy=length(unique(reef_occs[[3]]$site_dmy)),
                                                                     X1=X1,
                                                                     Z1=ncol(X1),
                                                                     X2=X2,
                                                                     Z2=ncol(X2),
                                                                     TT=26,
                                                                     N_yr1=length(unique(rvc_occs[[3]]$YEAR)),
                                                                     yr_index1=sort(unique(as.numeric(factor(rvc_occs[[3]]$YEAR)))),
                                                                     year_id1=as.numeric(factor(rvc_occs[[3]]$YEAR)),
                                                                     N_yr2=length(unique(reef_occs[[3]]$year)),
                                                                     yr_index2=sort(unique(as.numeric(factor(reef_occs[[3]]$year)))),
                                                                     year_id2=as.numeric(factor(reef_occs[[3]]$year))),
                           pars = c('a_hab1','a_hab2','sd_hab1','sd_hab2','sd_site','sd_dv','sd_dmy','sd_r1','sd_r2','sd_q','a','x','a_yr1','a_yr2'),
                           control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 150, chains = 4, iter = 450, thin = 1)


X1<- matrix(data=c(scale(as.numeric(rvc_occs[[5]]$DEPTH))),ncol=1,nrow=nrow(rvc_occs[[5]]))
X2<- matrix(data=c(scale(as.numeric(reef_occs[[5]]$btime)),scale(as.numeric(reef_occs[[5]]$averagedepth)),scale(as.numeric(reef_occs[[5]]$visibility)),scale(as.numeric(reef_occs[[5]]$current)),reef_occs[[5]]$exp_binary),ncol=5,nrow=nrow(reef_occs[[5]]))


test_comb_rf<- rstan::stan(model_code = abund_test_SS_comb, data = list(y1 = rvc_occs[[5]]$NUM.total2,
                                                                     y2 =reef_occs[[5]]$abundance2,
                                                                     N1 = nrow(rvc_occs[[5]]),
                                                                     N2 = nrow(reef_occs[[5]]),
                                                                     N_hab1 = length(unique(rvc_occs[[5]]$HAB_CD2)),
                                                                     hab_class1=as.numeric(factor(rvc_occs[[5]]$HAB_CD2)),
                                                                     N_hab2 = length(unique(reef_occs[[5]]$hab_class2)),
                                                                     hab_class2=as.numeric(factor(reef_occs[[5]]$hab_class2)),
                                                                     site=as.numeric(factor(reef_occs[[5]]$geogr)),
                                                                     N_site=length(unique(reef_occs[[5]]$geogr)),
                                                                     diver=as.numeric(factor(reef_occs[[5]]$fish_memberid)),
                                                                     N_dv=length(unique(reef_occs[[5]]$fish_memberid)),
                                                                     dmy=as.numeric(factor(reef_occs[[5]]$site_dmy)),
                                                                     N_dmy=length(unique(reef_occs[[5]]$site_dmy)),
                                                                     K=length(unique(reef_occs[[5]]$abundance)),
                                                                     X1=X1,
                                                                     Z1=ncol(X1),
                                                                     X2=X2,
                                                                     Z2=ncol(X2),
                                                                     TT=26,
                                                                     N_yr1=length(unique(rvc_occs[[5]]$YEAR)),
                                                                     yr_index1=sort(unique(as.numeric(factor(rvc_occs[[5]]$YEAR)))),
                                                                     year_id1=as.numeric(factor(rvc_occs[[5]]$YEAR)),
                                                                     N_yr2=length(unique(reef_occs[[5]]$year)),
                                                                     yr_index2=sort(unique(as.numeric(factor(reef_occs[[5]]$year)))),
                                                                     year_id2=as.numeric(factor(reef_occs[[5]]$year))),
                        pars = c('c','a_hab1','a_hab2','sd_hab1','sd_hab2','sd_site','sd_dv','sd_dmy','sd_r1','sd_r2','sd_q','x','a','log_lik'),
                        control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 150, chains = 2, iter = 300, thin = 1)

test_sep_rf<- rstan::stan(model_code = abund_test_SS_sep, data = list(y1 = rvc_occs[[5]]$NUM.total2,
                                                                   y2 =reef_occs[[5]]$abundance2,
                                                                   N1 = nrow(rvc_occs[[5]]),
                                                                   N2 = nrow(reef_occs[[5]]),
                                                                   N_hab1 = length(unique(rvc_occs[[5]]$HAB_CD2)),
                                                                   hab_class1=as.numeric(factor(rvc_occs[[5]]$HAB_CD2)),
                                                                   N_hab2 = length(unique(reef_occs[[5]]$hab_class2)),
                                                                   hab_class2=as.numeric(factor(reef_occs[[5]]$hab_class2)),
                                                                   site=as.numeric(factor(reef_occs[[5]]$geogr)),
                                                                   N_site=length(unique(reef_occs[[5]]$geogr)),
                                                                   diver=as.numeric(factor(reef_occs[[5]]$fish_memberid)),
                                                                   N_dv=length(unique(reef_occs[[5]]$fish_memberid)),
                                                                   dmy=as.numeric(factor(reef_occs[[5]]$site_dmy)),
                                                                   N_dmy=length(unique(reef_occs[[5]]$site_dmy)),
                                                                   K=length(unique(reef_occs[[5]]$abundance)),
                                                                   X1=X1,
                                                                   Z1=ncol(X1),
                                                                   X2=X2,
                                                                   Z2=ncol(X2),
                                                                   TT=26,
                                                                   N_yr1=length(unique(rvc_occs[[5]]$YEAR)),
                                                                   yr_index1=sort(unique(as.numeric(factor(rvc_occs[[5]]$YEAR)))),
                                                                   year_id1=as.numeric(factor(rvc_occs[[5]]$YEAR)),
                                                                   N_yr2=length(unique(reef_occs[[5]]$year)),
                                                                   yr_index2=sort(unique(as.numeric(factor(reef_occs[[5]]$year)))),
                                                                   year_id2=as.numeric(factor(reef_occs[[5]]$year))),
                       pars = c('c','a_hab1','a_hab2','sd_hab1','sd_hab2','sd_site','sd_dv','sd_dmy','sd_r1','sd_r2','sd_q1','sd_q2','x1','x2','a_yr1','a_yr2','log_lik'),
                       control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 150, chains = 2, iter = 600, thin = 1)

shinystan::launch_shinystan(test_comb_rf)

params_1<- rstan::extract(test_comb_rf)
params_2<- rstan::extract(test_sep_rf)

loo1= loo::loo(test_comb_rf)
start.time <- Sys.time()
system.time(loo2<- loo::loo(test_sep_rf))

time.taken <- end.time - start.time
time.taken
start.time <- Sys.time()
loo_comp<- loo::loo_compare(loo1,loo2)
time.taken <- end.time - start.time
time.taken

params_sep<- rstan::extract(test_sep_rf)

TS_stan_abund_plot_MARSS(ts1=rvc_ts[[5]],ts2=reef_ts[[5]],sp='Rock Beauty',GZ='Key Largo',params=params_sep)


test_sep_rf_trend<- rstan::stan(model_code = abund_test_SS_sep_trend, data = list(y1 = rvc_occs[[5]]$NUM.total2,
                                                                      y2 =reef_occs[[5]]$abundance2,
                                                                      N1 = nrow(rvc_occs[[5]]),
                                                                      N2 = nrow(reef_occs[[5]]),
                                                                      N_hab1 = length(unique(rvc_occs[[5]]$HAB_CD2)),
                                                                      hab_class1=as.numeric(factor(rvc_occs[[5]]$HAB_CD2)),
                                                                      N_hab2 = length(unique(reef_occs[[5]]$hab_class2)),
                                                                      hab_class2=as.numeric(factor(reef_occs[[5]]$hab_class2)),
                                                                      site=as.numeric(factor(reef_occs[[5]]$geogr)),
                                                                      N_site=length(unique(reef_occs[[5]]$geogr)),
                                                                      diver=as.numeric(factor(reef_occs[[5]]$fish_memberid)),
                                                                      N_dv=length(unique(reef_occs[[5]]$fish_memberid)),
                                                                      dmy=as.numeric(factor(reef_occs[[5]]$site_dmy)),
                                                                      N_dmy=length(unique(reef_occs[[5]]$site_dmy)),
                                                                      K=length(unique(reef_occs[[5]]$abundance)),
                                                                      X1=X1,
                                                                      Z1=ncol(X1),
                                                                      X2=X2,
                                                                      Z2=ncol(X2),
                                                                      TT=26,
                                                                      N_yr1=length(unique(rvc_occs[[5]]$YEAR)),
                                                                      yr_index1=sort(unique(as.numeric(factor(rvc_occs[[5]]$YEAR)))),
                                                                      year_id1=as.numeric(factor(rvc_occs[[5]]$YEAR)),
                                                                      N_yr2=length(unique(reef_occs[[5]]$year)),
                                                                      yr_index2=sort(unique(as.numeric(factor(reef_occs[[5]]$year)))),
                                                                      year_id2=as.numeric(factor(reef_occs[[5]]$year))),
                          pars = c('c','a_hab1','a_hab2','sd_hab1','sd_hab2','sd_site','sd_dv','sd_dmy','sd_r1','sd_r2','sd_q1','sd_q2','u1','u2','x1','x2','a_yr1','a_yr2'),
                          control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 150, chains = 4, iter = 450, thin = 1)

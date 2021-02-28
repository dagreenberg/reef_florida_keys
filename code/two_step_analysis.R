rm(list=ls())
setwd("C:/Users/14388/Desktop/reef_florida_keys_data")
library(rstan);library(loo);library(bayesplot);library(dplyr);library(tidyverse);library(stringr);library(lubridate);library(rlist)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

####Functions####
rvc_filter = function(x,GZ,sp){
  x$SSU_YEAR<- paste(x$PRIMARY_SAMPLE_UNIT,x$STATION_NR,x$YEAR,sep='_')
  x1= x %>% subset(region.id==GZ) %>% select(SSU_YEAR,SPECIES_CD,everything())
  x2= complete(x1,SSU_YEAR,nesting(SPECIES_CD),fill=list(NUM=0))
  zeros = anti_join(x2,x1)
  zeros[,3:29]<- x1[match(zeros$SSU_YEAR,x1$SSU_YEAR),3:29]
  x3= rbind(x1,zeros)
  x3$HAB_CD2<- gsub('\\_.*','',x3$HABITAT_CD)
  rvc_occs<- list()
  for(i in 1:nrow(sp)){
    x4= subset(x3,SPECIES_CD==sp$rvc_code[i])
    x5=  x4 %>% dplyr::group_by(SSU_YEAR) %>%
      dplyr::summarise(NUM.total=sum(NUM),occ=NA) %>% #Sums up the number of counts
      mutate(occ=ifelse(NUM.total>0,1,0)) %>% arrange(SSU_YEAR) #Also scores presence/absence at the SSU level
    x5[,4:32]<- x4[match(x5$SSU_YEAR,x4$SSU_YEAR),2:30]
    x5<- transform(x5,psu_id=match(LAT_LON,unique(LAT_LON)))
    x5$NUM.total<- round(x5$NUM.total)
    rvc_occs[[i]]=x5
  }
  return(rvc_occs)
}

ts_rvc = function(x){ #Takes the output from the previous function
  ts<- list()
  for(i in 1:length(x)){
    x1 = x[[i]]
    x2= x1 %>% group_by(YEAR) %>% summarize(n.occ=sum(occ),n.surv=n(),p.occ=n.occ/n.surv,sp=unique(SPECIES_CD))
    x3= x1 %>% group_by(YEAR,PRIMARY_SAMPLE_UNIT) %>% summarize(psu_abund=mean(NUM.total)) %>% group_by(YEAR) %>% summarize(mean_abund=mean(psu_abund),sd_abund=sd(psu_abund))
    for(z in 1:nrow(x2)){
      if(x2$p.occ[z]==0){
        x2$p.occ[z]=NA
      }
    }
    x4=left_join(x2,x3) %>% complete(YEAR=seq(1993,2018))
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
  TempDat4$site_dmy<- paste(TempDat4$geogr,TempDat4$day,TempDat4$month,TempDat4$year,sep='')
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


###3. Geographic filtering ####
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

#### 4. Creating RVC time-series ####
#Fish data from REEF - remove ultra rare and basket species designations
fish_reef<- read.csv("Caribbean_fish_trait_matrix.csv") #fish species found in the Tropical Western Atlantic
fish_reef<- subset(fish_reef,expert_sighting_freq>1) #Take out the very rare species
fish_rvc<- read.csv("Florida_keys_taxonomic_data.csv")
fish_rvc<- subset(fish_rvc,gsub('.*\\ ', '', fish_rvc$SCINAME)!='sp.') #remove unknown species
fish_rvc<- subset(fish_rvc, SCINAME %in% fish_reef$sciname2)
m<- match(fish_reef$sciname2,fish_rvc$SCINAME)
fish_reef$rvc_code<- fish_rvc$SPECIES_CD[m]

fk_93_18<- subset(fk_79_18,YEAR>=1993) #Subset for the dataset from 1993 to match the first year of REEF surveys


rvc_occs<- rvc_filter(fk_93_18,GZ='3403',sp=fish_reef)
rvc_ts<- ts_rvc(rvc_occs)
rvc_ts_filter<- rlist::list.filter(rvc_ts,length(na.omit(p.occ))>16)
rvc.green<- do.call(rbind, lapply(rvc_ts_filter, data.frame, stringsAsFactors=FALSE))
rvc.green.sp<- unique(rvc.green$sp)
fish_reef_trim<- subset(fish_reef, rvc_code %in% rvc.green.sp)

rvc_occs<- rvc_filter(fk_93_18,GZ='3403',sp=fish_reef_trim)
reef_occs<- reef_filter(R,GZ='3403',sp=fish_reef_trim,geog=reef_geog_3403)
reef_ts<- ts_reef(reef_occs,sp=fish_reef_trim)

reef_mod<- vector(mode='list',length=nrow(fish_reef))
rvc_mod<- vector(mode='list',length=nrow(fish_reef))
ts_comp<- vector(mode='list',length=nrow(fish_reef))
ts_orig<- vector(mode='list',length=nrow(fish_reef))
library(glmmTMB); library(MARSS)
setwd("C:/Users/14388/Desktop/reef_florida_keys_data/occ ts plots feb")
mars_3403<- data.frame(SP=NA,conv.1=NA,conv.2=NA,tier1.AICc=NA,tier2.AICc=NA,Q1=NA,Q2.rvc=NA,Q2.reef=NA,R1.rvc=NA,R1.reef=NA,R2.rvc=NA,R2.reef=NA,dAIC1=NA,dAIC2=NA,mod=NA,years.rvc=NA,years.reef=NA,rvc.b0=NA,reef.b0=NA,rvc.var.hab=NA,reef.var.hab=NA,rvc.var.yr=NA,reef.var.year=NA,reef.var.site=NA,reef.var.dmy=NA,rvc.b.depth=NA,reef.b.depth=NA,reef.b.btime=NA,reef.b.vis=NA,reef.b.curr=NA,rvc.b.cont=NA,reef.b.cont=NA,rvc.b.isol=NA,reef.b.isol=NA,rvc.b.rubb=NA,reef.b.hubb=NA,rvc.b.sg=NA,reef.b.sg=NA)
for(i in 1:nrow(fish_reef)){
  spp_rvc<- rvc_occs[[i]]
  spp_reef<- reef_occs[[i]]
  
  rvc_mod[[i]]<-glmmTMB(occ~scale(DEPTH) +(1|HAB_CD2)+(1|YEAR),family = binomial,data=spp_rvc)
  reef_mod[[i]]<-glmmTMB(occ~scale(as.numeric(btime))+scale(as.numeric(averagedepth))+scale(as.numeric(visibility))+scale(as.numeric(current))+(1|hab_class2)+(1|hab_class2:geogr)+(1|year)+(1|site_dmy),family = binomial,data=spp_reef)
  
  rvc_years<- data.frame(year=as.numeric(rownames(coef(rvc_mod[[i]])$cond$YEAR)),year_prob=coef(rvc_mod[[i]])$cond$YEAR[,1])
  rvc_years<- as.data.frame(complete(rvc_years,year=seq(1993,2018)))
  
  ts_comp[[i]]<- t(rvc_years$year_prob)
  ts_comp[[i]]<- rbind(ts_comp[[i]],t(coef(reef_mod[[i]])$cond$year[,1]))
  colnames(ts_comp[[i]])<- rvc_years$year
  
  ts_orig[[i]]<- t(rvc_ts[[i]]$p.occ)
  ts_orig[[i]]<- rbind(ts_orig[[i]],t(reef_ts[[i]]$p.occ))
  
  Z_1=factor(c('rvc_reef','rvc_reef')) #Common state process for both time-series
  Z_2=factor(c('rvc','reef')) #Different state process for each time-series
  
  #tier 1 - common state process
  tier_1<- list(Z=Z_1,
                Q = 'diagonal and equal',
                R = 'diagonal and unequal',
                U = 'zero',
                A = 'scaling',
                x0 = 'equal',
                tinitx=0)
  #tier 3
  tier_3<- list(Z=Z_2,
                Q = 'diagonal and unequal',
                R = 'diagonal and unequal',
                U = 'unequal',
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
 
  fit_1<-  MARSS(ts_comp[[i]], model = tier_1,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1,allow.degen=FALSE),method='kem')

  fit_3<-  MARSS(ts_comp[[i]], model = tier_3,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1,allow.degen=FALSE),method='kem')
  params.1<- MARSSparamCIs(fit_1)
  params.3<- MARSSparamCIs(fit_3)
  
  mars_3403[i,1]=fish_reef$commonname[i]
  mars_3403[i,2]=ifelse(is.null(fit_1$errors)==T,1,0)
  mars_3403[i,3]=ifelse(is.null(fit_3$errors)==T,1,0)  
  mars_3403[i,4]=fit_1$AICc
  mars_3403[i,5]=fit_3$AICc  
  mars_3403[i,6]=params.1$parMean[4]
  mars_3403[i,7]=params.3$parMean[3]
  mars_3403[i,8]=params.3$parMean[4]
  mars_3403[i,9]=params.1$parMean[2]
  mars_3403[i,10]=params.1$parMean[3] 
  mars_3403[i,11]=params.3$parMean[1]
  mars_3403[i,12]=params.3$parMean[2]
  mars_3403[i,13]=mars_3403[i,4]-min(mars_3403[i,4:5])
  mars_3403[i,14]=mars_3403[i,5]-min(mars_3403[i,4:5])
  if(mars_3403[i,13]==0){mars_3403[i,15]=1}
  if(mars_3403[i,14]==0){mars_3403[i,15]=2}
  mars_3403[i,16]=length(na.omit(ts_comp[[i]][1,]))
  mars_3403[i,17]=length(na.omit(ts_comp[[i]][2,]))
  mars_3403[i,18]=rvc_mod[[i]]$fit$par[1]
  mars_3403[i,19]=reef_mod[[i]]$fit$par[1]
  mars_3403[i,20]=VarCorr(rvc_mod[[i]])$cond$HAB_CD2[1]
  mars_3403[i,21]=VarCorr(reef_mod[[i]])$cond$hab_class2[1]
  mars_3403[i,22]=VarCorr(rvc_mod[[i]])$cond$YEAR[1]
  mars_3403[i,23]=VarCorr(reef_mod[[i]])$cond$year[1]
  mars_3403[i,24]=VarCorr(reef_mod[[i]])$cond$`hab_class2:geogr`[1]
  mars_3403[i,25]=VarCorr(reef_mod[[i]])$cond$site_dmy[1]
  mars_3403[i,26]=fixef(rvc_mod[[i]])$cond[2]
  mars_3403[i,27]=fixef(reef_mod[[i]])$cond[3]
  mars_3403[i,28]=fixef(reef_mod[[i]])$cond[2]
  mars_3403[i,29]=fixef(reef_mod[[i]])$cond[4]
  mars_3403[i,30]=fixef(reef_mod[[i]])$cond[5]
  mars_3403[i,31]=ranef(rvc_mod[[i]])$cond$HAB_CD2[1,]
  mars_3403[i,32]=ranef(reef_mod[[i]])$cond$hab_class2[1,]
  mars_3403[i,33]=ranef(rvc_mod[[i]])$cond$HAB_CD2[2,]
  mars_3403[i,34]=ranef(reef_mod[[i]])$cond$hab_class2[2,]
  mars_3403[i,35]=ranef(rvc_mod[[i]])$cond$HAB_CD2[3,]
  mars_3403[i,36]=ranef(reef_mod[[i]])$cond$hab_class2[3,]
  mars_3403[i,37]=ranef(rvc_mod[[i]])$cond$HAB_CD2[4,]
  mars_3403[i,38]=ranef(reef_mod[[i]])$cond$hab_class2[4,]
  
  TS_occ_plot_MARSS(ts1=ts_orig[[i]],ts2=ts_comp[[i]],sp=fish_reef$commonname[i],GZ='Key Largo',mod=mars_3403[i,15])
  dev.off()
  print(i)
  
}
write.csv(mars_3403,'GLMM_smoothed_occ_ts_3403.csv')

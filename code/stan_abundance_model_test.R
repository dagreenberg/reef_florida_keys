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
    x4=left_join(x2,x3)
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
  for(i in 1:nrow(sp)){
    occ_list[[i]]<- subset(TempDat4,speciesid==sp$speciesid[i]) #Subset out each species in the provided dataframe
  }
  return(occ_list)
}

pool_sd = function (sd, n, m){
  dat<- data.frame(sd=sd,n=n,m=m)
  dat<- dat[is.na(dat$sd)==F,]
  q<- NA
  for(i in 1:nrow(dat)){
    q[i] = (dat$n[i]-1)*dat$sd[i]^2 + dat$n[i]*dat$m^2
  }
  qc = sum(q)
  m_a=mean(dat$m)
  sc = sqrt((qc - (sum(dat$n))*dat$m_a^2)/(sum(dat$n)-1) )
  return(sc)  
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

comp_plot_rvc = function(ts,mod_mat,sp,GZ){
  plot(ts$p.occ~ts$YEAR,type='n',xlab='Year',ylab='Probability of occurrence',ylim=c(0,1),bty='l',main=paste(sp,GZ,sep=' - '))
  lines(ts$p.occ~ts$YEAR,lwd=2,col='darkblue')
  points(ts$p.occ~ts$YEAR,pch=21,col='darkblue',bg='white',cex=1.5,lwd=1.5)
  points(mod_mat$median~ts$YEAR,pch=21,col=adjustcolor('dodgerblue3',alpha.f = 0.8),cex=1.2,lwd=1.2)
  lines(mod_mat$median~ts$YEAR,lwd=1.5,col=adjustcolor('dodgerblue3',alpha.f = 0.8))
  x.polygon <- c(ts$YEAR, rev(ts$YEAR)) # Define a polygon x value for adding to a plot
  y.polygon <- c(mod_mat$l.95, rev(mod_mat$u.95)) # Define a polygon y value for adding to a plot
  polygon(x.polygon, y.polygon, col = adjustcolor('dodgerblue3', alpha = 0.2), border=NA) # Add uncertainty polygon
  
}

comp_plot_SS_rvc = function(ts,x_mat,y_mat,sp,GZ){
  plot(ts$p.occ~ts$YEAR,type='n',xlab='Year',ylab='Probability of occurrence',ylim=c(0,1),bty='l',main=paste(sp,GZ,sep=' - '))
  lines(ts$p.occ~ts$YEAR,lwd=2,col='darkblue')
  points(ts$p.occ~ts$YEAR,pch=21,col='darkblue',bg='white',cex=1.5,lwd=1.5)
  points(x_mat$median~seq(min(ts$YEAR),max(ts$YEAR)),pch=21,col=adjustcolor('dodgerblue3',alpha.f = 0.8),cex=1.2,lwd=1.2)
  lines(x_mat$median~seq(min(ts$YEAR),max(ts$YEAR)),lwd=1.5,col=adjustcolor('dodgerblue3',alpha.f = 0.8))
  points(y_mat$median~ts$YEAR,pch=21,col=adjustcolor('cyan',alpha.f = 0.8),cex=1.2,lwd=1.2)
  lines(y_mat$median~ts$YEAR,lwd=1.5,col=adjustcolor('cyan',alpha.f = 0.8),lty=5)
  x.polygon <- c(seq(min(ts$YEAR),max(ts$YEAR)), rev(seq(min(ts$YEAR),max(ts$YEAR)))) # Define a polygon x value for adding to a plot
  y.polygon <- c(x_mat$l.95, rev(x_mat$u.95)) # Define a polygon y value for adding to a plot
  polygon(x.polygon, y.polygon, col = adjustcolor('dodgerblue3', alpha = 0.1), border=NA) # Add uncertainty polygon
}

comp_plot_SS_reef = function(ts,x_mat,y_mat,sp,GZ){
  plot(ts$p.occ~ts$year,type='n',xlab='Year',ylab='Probability of occurrence',ylim=c(0,1),bty='l',main=paste(sp,GZ,sep=' - '))
  lines(ts$p.occ~ts$year,lwd=2,col='darkblue')
  points(ts$p.occ~ts$year,pch=21,col='darkblue',bg='white',cex=1.5,lwd=1.5)
  points(x_mat$median~seq(min(ts$year),max(ts$year)),pch=21,col=adjustcolor('dodgerblue3',alpha.f = 0.8),cex=1.2,lwd=1.2)
  lines(x_mat$median~seq(min(ts$year),max(ts$year)),lwd=1.5,col=adjustcolor('dodgerblue3',alpha.f = 0.8))
  points(y_mat$median~ts$year,pch=21,col=adjustcolor('cyan',alpha.f = 0.8),cex=1.2,lwd=1.2)
  lines(y_mat$median~ts$year,lwd=1.5,col=adjustcolor('cyan',alpha.f = 0.8),lty=5)
  x.polygon <- c(seq(min(ts$year),max(ts$year)), rev(seq(min(ts$year),max(ts$year)))) # Define a polygon x value for adding to a plot
  y.polygon <- c(x_mat$l.95, rev(x_mat$u.95)) # Define a polygon y value for adding to a plot
  polygon(x.polygon, y.polygon, col = adjustcolor('dodgerblue3', alpha = 0.1), border=NA) # Add uncertainty polygon
  
}

comp_plot_SS_comb = function(ts1,ts2,x_mat,y_mat1,y_mat2,sp,GZ){
  plot(ts1$p.occ~ts1$YEAR,type='n',xlab='Year',ylab='Probability of occurrence',ylim=c(0,1),bty='l',main=paste(sp,GZ,sep=' - '))
  lines(ts1$p.occ~ts1$YEAR,lwd=2,col='darkblue')
  points(ts1$p.occ~ts1$YEAR,pch=21,col='darkblue',bg='white',cex=1.5,lwd=1.5)
  lines(ts2$p.occ~ts2$year,lwd=2,col='darkred')
  points(ts2$p.occ~ts2$year,pch=21,col='darkred',bg='white',cex=1.5,lwd=1.5)
  
  
  points(x_mat$median~seq(min(ts1$YEAR),max(ts1$YEAR)),pch=21,col=adjustcolor('dodgerblue3',alpha.f = 0.8),cex=1.2,lwd=1.2)
  lines(x_mat$median~seq(min(ts1$YEAR),max(ts1$YEAR)),lwd=1.5,col=adjustcolor('dodgerblue3',alpha.f = 0.8))
  points(y_mat1$median~ts1$YEAR,pch=21,col=adjustcolor('darkblue',alpha.f = 0.8),cex=1.2,lwd=1.2)
  lines(y_mat1$median~ts1$YEAR,lwd=1.5,col=adjustcolor('darkblue',alpha.f = 0.8),lty=5)
  points(y_mat2$median~ts2$year,pch=21,col=adjustcolor('darkred',alpha.f = 0.8),cex=1.2,lwd=1.2)
  lines(y_mat2$median~ts2$year,lwd=1.5,col=adjustcolor('darkred',alpha.f = 0.8),lty=5)
  
  
  x.polygon <- c(seq(min(ts1$YEAR),max(ts1$YEAR)), rev(seq(min(ts1$YEAR),max(ts1$YEAR)))) # Define a polygon x value for adding to a plot
  y.polygon <- c(x_mat$l.95, rev(x_mat$u.95)) # Define a polygon y value for adding to a plot
  polygon(x.polygon, y.polygon, col = adjustcolor('dodgerblue3', alpha = 0.1), border=NA) # Add uncertainty polygon
  
}

comp_plot_SS_sep = function(ts1,ts2,x_mat1,x_mat2,y_mat1,y_mat2,sp,GZ){
  plot(ts1$p.occ~ts1$YEAR,type='n',xlab='Year',ylab='Probability of occurrence',ylim=c(0,1),bty='l',main=paste(sp,GZ,sep=' - '))
  lines(ts1$p.occ~ts1$YEAR,lwd=2,col='darkblue')
  points(ts1$p.occ~ts1$YEAR,pch=21,col='darkblue',bg='white',cex=1.5,lwd=1.5)
  lines(ts2$p.occ~ts2$year,lwd=2,col='darkred')
  points(ts2$p.occ~ts2$year,pch=21,col='darkred',bg='white',cex=1.5,lwd=1.5)
  
  
  points(x_mat1$median~seq(min(ts1$YEAR),max(ts1$YEAR)),pch=21,col=adjustcolor('dodgerblue3',alpha.f = 0.8),cex=1.2,lwd=1.2)
  lines(x_mat1$median~seq(min(ts1$YEAR),max(ts1$YEAR)),lwd=1.5,col=adjustcolor('dodgerblue3',alpha.f = 0.8))
  
  points(x_mat2$median~seq(min(ts1$YEAR),max(ts1$YEAR)),pch=21,col=adjustcolor('dodgerblue3
                                                                               ',alpha.f = 0.8),cex=1.2,lwd=1.2)
  lines(x_mat2$median~seq(min(ts1$YEAR),max(ts1$YEAR)),lwd=1.5,col=adjustcolor('dodgerblue3',alpha.f = 0.8))
  points(y_mat1$median~ts1$YEAR,pch=21,col=adjustcolor('cyan',alpha.f = 0.8),cex=1.2,lwd=1.2)
  lines(y_mat1$median~ts1$YEAR,lwd=1.5,col=adjustcolor('cyan',alpha.f = 0.8),lty=5)
  points(y_mat2$median~ts2$year,pch=21,col=adjustcolor('cyan',alpha.f = 0.8),cex=1.2,lwd=1.2)
  lines(y_mat2$median~ts2$year,lwd=1.5,col=adjustcolor('cyan',alpha.f = 0.8),lty=5)
  
  
  x.polygon <- c(seq(min(ts1$YEAR),max(ts1$YEAR)), rev(seq(min(ts1$YEAR),max(ts1$YEAR)))) # Define a polygon x value for adding to a plot
  y1.polygon <- c(x_mat1$l.95, rev(x_mat1$u.95)) # Define a polygon y value for adding to a plot
  y2.polygon <- c(x_mat2$l.95, rev(x_mat2$u.95)) # Define a polygon y value for adding to a plot
  polygon(x.polygon, y1.polygon, col = adjustcolor('dodgerblue3', alpha = 0.1), border=NA) # Add uncertainty polygon
  polygon(x.polygon, y2.polygon, col = adjustcolor('dodgerblue3', alpha = 0.1), border=NA) # Add uncertainty polygon
  
}

comp_plot_reef = function(ts,mod_mat,sp,GZ){
  plot(ts$p.occ~ts$year,type='n',xlab='Year',ylab='Probability of occurrence',ylim=c(0,1),bty='l',main=paste(sp,GZ,sep=' - '))
  lines(ts$p.occ~ts$year,lwd=2,col='darkblue')
  points(ts$p.occ~ts$year,pch=21,col='darkblue',bg='white',cex=1.5,lwd=1.5)
  points(mod_mat$median~ts$year,pch=21,col=adjustcolor('dodgerblue3',alpha.f = 0.8),cex=1.2,lwd=1.2)
  lines(mod_mat$median~ts$year,lwd=1.5,col=adjustcolor('dodgerblue3',alpha.f = 0.8))
  x.polygon <- c(ts$year, rev(ts$year)) # Define a polygon x value for adding to a plot
  y.polygon <- c(mod_mat$l.95, rev(mod_mat$u.95)) # Define a polygon y value for adding to a plot
  polygon(x.polygon, y.polygon, col = adjustcolor('dodgerblue3', alpha = 0.2), border=NA) # Add uncertainty polygon
  
}



hab_breakdown = function(x){
  hab = x %>% group_by(YEAR) %>% count(HABITAT_CD) %>% arrange(YEAR,HABITAT_CD)
  sum_yr = hab %>% group_by(YEAR) %>% summarize(n.tot=sum(n))
  comb<- right_join(hab,sum_yr,by='YEAR')
  comb$prop<- comb$n/comb$n.tot
  return(comb)
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
reef_occs<- reef_filter(R,GZ='3403',sp=fish_reef,geog=reef_geog_3403)
reef_ts<- ts_reef(reef_occs,sp=fish_reef)


library(glmmTMB)
glmmTMB(NUM.total~1+offset(log(n.assess)),data=rvc_occs[[1]],family=nbinom2)

rep_survs<- R %>% group_by(geogr,month,day,year) %>% summarize(n_survs=n_distinct(formid))

####Stan model
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

zip_test_rvc<-"data{
  int<lower=1> N;//number of observations (SSU surveys)
  int y[N]; //rounded counts from each survey
  int<lower=0> N_hab; //number of habitat classes
  int<lower=1,upper=N_hab> hab_class[N]; // vector of habitat class identities
  int<lower=0> N_yr; //number of years
  int<lower=1,upper=N_yr> year[N]; // vector of year
  int K; // columns in the covariate matrix
  matrix[N,K] X; // design matrix X
}
functions {
  int num_zero(int[] y) {
    int sum = 0;
    for (n in 1:size(y))
      sum += (y[n] == 0);
    return sum;
  }
}
transformed data {
  int<lower=0, upper=N> N0 = num_zero(y);
  int<lower=0, upper=N> Ngt0 = N - N0;
  int<lower=1> y_nz[N - num_zero(y)];
  {
    int pos = 1;
    for (n in 1:N) {
      if (y[n] != 0) {
        y_nz[pos] = y[n];
        pos += 1;
      }
    }
  }
}
transformed parameters{
  vector[Ngt0] lambda;
  
  lambda = alpha + z_yr[year[n]]*tau_yr + z_hab[hab_class[n]]*tau_hab + X[n,]*beta;

}  
parameters {
  real alpha; 
  //global intercept
  real<lower = 0, upper = 1> theta;
  //dispersion parameter
  real beta;
  //depth covariate
  
  //deviations from intercept
  real z_hab[N_hab]; //deviation between habitats
  real z_yr[N_yr]; //deviation between years
 
  
  //st dev on the deviations
  real<lower = 0> tau_hab; //sigma on habitat
  real<lower = 0> tau_yr; //sigma on year
}
model{
  //priors
  alpha ~ normal(0,3);
  beta ~ normal(0,2);
  
  //shape parameter
  theta ~ beta(2,2);
  
  //varying intercepts
  z_hab ~ normal(0,3);
  z_yr ~ student_t(5, 0, 3);
  tau_yr ~ cauchy(0,3);
  tau_hab ~ cauchy(0,3);
  
  N0 ~ binomial(N, theta);
  
  y_nz ~ poisson(lambda);
  target += -Ngt0 * log1m_exp(-lambda);
  }
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
  vector[N_hab] z_hab; //deviation between habitats
  vector[N_yr] z_yr; //deviation between years
  vector[TT] z_dev; //state deviations

  
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
  
  x[1] = x0 + z_dev[1]*sd_q;
   
  for(t in 2:TT){
    x[t] = x[t-1] + z_dev[t]*sd_q;
  }
   
  for(i in 1:N_yr){
    a_yr[i] = x[yr_index[i]] + sd_r*z_yr[year_id[i]]; 
  }
  
}  
  
model{
  //priors
  x0 ~ normal(0,10);
  beta ~ normal(0,2);
 
  //shape parameter
  recip_phi ~ cauchy(0, 5);
 
  //standard deviations
  sd_hab ~ cauchy(0,5);
  sd_r ~ cauchy(0,5);
  sd_q ~ cauchy(0,5);
  
  //varying intercepts
  z_hab ~ normal(0,5);
  z_yr ~ normal(0,5);
 
 
  //draw presences
    y ~ neg_binomial_2_log(a_yr[year_id] + z_hab[hab_class]*sd_hab + X*beta,phi);
}
"

ord_test_reef<-"data{
  int<lower=1> N;//number of observations (SSU surveys)
  int y[N]; //abundance category for each survey
  int<lower=0> N_hab; //number of habitat classes
  int<lower=1,upper=N_hab> hab_class[N]; // vector of habitat class identities
  int<lower=0> N_site; //number of sites
  int<lower=1,upper=N_site> site[N]; // vector of site identities
  int<lower=0> N_dv; //number of divers
  int<lower=1,upper=N_dv> diver[N]; // vector of diver identities
  int<lower=0> N_yr; //number of years
  int<lower=1,upper=N_yr> year[N]; // vector of year
  int<lower=0> N_dmy; //number of day-month-year site samples
  int<lower=1,upper=N_yr> dmy[N]; // vector of day-month-year site samples
  int K; // columns in the covariate matrix
  matrix[N,K] X; // design matrix X
  int J; //ordinal levels
}
parameters {
  ordered[J-1] c; //cutpoints
  real<lower = 0> c_sigma; //sigma for cutpoints
  real c_mu; //mean for cutpoints
  
  //deviations from intercept
  real a_hab[N_hab]; //deviation between habitats
  real a_yr[N_yr]; //deviation between years
  real a_site[N_site]; //deviation between sites
  real a_dv[N_dv]; //deviation between divers
  
  //covariates for effort variables
  vector[K] beta;
  
  //st dev on the deviations
  real<lower = 0> sigma_hab;
  real<lower = 0> sigma_yr;
  real<lower = 0> sigma_site;
  real<lower = 0> sigma_dv;
  real<lower = 0> sigma_dmy;
  
}

transformed parameters{
  vector[N] eta;
  
  for(n in 1:N){
    eta[n] = alpha + a_hab[hab_class[n]] + a_yr[year[n]] + a_site[site[n]]+ a_dv[diver[n]] + X[n,]*beta;
  }
}

model{
  //priors
  alpha ~ normal(0,10);
  beta ~ normal(0,2);
  
  //standard deviations
  sigma_hab ~ cauchy(0, 5);
  sigma_yr ~ cauchy(0, 5);
  sigma_site ~ cauchy(0, 5);
  sigma_dv ~ cauchy(0, 5);
  
  //varying intercepts
  a_hab ~ normal(0, sigma_hab);
  a_yr ~ student_t(5, 0, sigma_yr);
  a_site ~ normal(0, sigma_site);
  a_dv ~ normal(0, sigma_dv);
 
 
  y ~ ordered_logistic(mu,c);
  
}
"

logit_test_SS_reef<-"data{
  int<lower=1> N;//number of observations (SSU surveys)
  int y[N]; //presence or absence on each survey
  int<lower=0> N_hab; //number of habitat classes
  int<lower=1,upper=N_hab> hab_class[N]; // vector of habitat class identities
  int<lower=0> N_site; //number of sites
  int<lower=1,upper=N_site> site[N]; // vector of site identities
  int<lower=0> N_dv; //number of divers
  int<lower=1,upper=N_dv> diver[N]; // vector of diver identities
  int<lower=0> N_yr; //number of years
  int yr_index[N_yr]; //index of years
  int<lower=1,upper=N_yr> year_id[N]; // vector of year
  int K; // columns in the covariate matrix
  matrix[N,K] X; // design matrix X
  int TT; //timespan
  
}
parameters {
  //global intercept
  real x[TT];
  real x0;
  
  //deviations from intercept
  real a_hab[N_hab]; //deviation between habitats
  real a_yr[N_yr]; //deviation between years
  real a_site[N_site]; //deviation between sites
  real a_dv[N_dv]; //deviation between divers
  
  //covariates for effort variables
  vector[K] beta;
  
  //st dev on the deviations
  real<lower = 0> sigma_hab;
  real<lower = 0> sd_q;
  real<lower = 0> sd_r;
  real<lower = 0> sigma_site;
  real<lower = 0> sigma_dv;
}


model{
  //priors
  x0 ~ normal(0, 10);
  beta ~ normal(0,2);
  
  //standard deviations
  sigma_hab ~ cauchy(0, 5);
  sd_q ~ cauchy(0, 5);
  sd_r ~ cauchy(0, 5);
  sigma_site ~ cauchy(0, 5);
  sigma_dv ~ cauchy(0, 5);
  
  //varying intercepts
  a_hab ~ normal(0, sigma_hab);
  a_site ~ normal(0, sigma_site);
  a_dv ~ normal(0, sigma_dv);
 
 
 //state process 
  x[1] ~ normal(x0,sd_q);
  
  for(t in 2:TT){
    x[t] ~ normal(x[t-1],sd_q);
  }
  
  //observation process for year
  for(i in 1:N_yr){
    a_yr[i]~ normal(x[yr_index[i]],sd_r); 
  }
 
  for(i in 1:N){
    y[i] ~ bernoulli_logit(a_yr[year_id[i]] + a_hab[hab_class[i]]+ a_site[site[i]]+ a_dv[diver[i]] + X[i,]*beta);
  }
  
 target +=  multinomial_lpmf(int[] y | vector theta)
}
"


logit_test_SS_comb<-"data{
  int TT; //timespan
  int N_rvc;//number of observations (SSU surveys)
  int y1[N_rvc]; //presence or absence on each survey
  int N_reef;//number of observations (REEF surveys)
  int y2[N_reef]; //presence or absence on each survey
  int N_hab1; //number of habitat classes
  int<lower=1,upper=N_hab1> hab_class1[N_rvc]; // vector of habitat class identities
  int N_hab2; //number of habitat classes
  int<lower=1,upper=N_hab2> hab_class2[N_reef]; // vector of habitat class identities
  int N_yr1; //number of years rvc
  int N_yr2; //number of years reef
  int yr_index1[N_yr1]; //index of years - RVC
  int<lower=1,upper=N_yr1> year_id1[N_rvc]; // vector of year ids
  int yr_index2[N_yr2]; //index of years - REEF
  int<lower=1,upper=N_yr2> year_id2[N_reef]; // vector of year ids
  int<lower=0> N_site; //number of sites in reef
  int<lower=1,upper=N_site> site[N_reef]; // vector of site identities for reef
  int<lower=0> N_dv; //number of divers in reef
  int<lower=1,upper=N_dv> diver[N_reef]; // vector of diver identities for reef
  int K; // columns in the covariate matrix
  matrix[N_reef,K] X; // design matrix X
}
parameters {
  //global intercept
  real x0;
  real x[TT];
  real beta;
  real s;
  
  //deviations from intercept
  real a_hab1[N_hab1]; //deviation between habitats
  real a_hab2[N_hab2]; //deviation between habitats
  real a_yr1[N_yr1]; //deviation between years
  real a_yr2[N_yr2]; //deviation between years
  real a_site[N_site]; //deviation between sites
  real a_dv[N_dv]; //deviation between divers
  
  
  //st dev on the deviations
  real<lower = 0> sd_hab1; //sigma on habitat
  real<lower = 0> sd_hab2; //sigma on habitat
  real<lower = 0> sd_r1; //sigma for observation error
  real<lower = 0> sd_r2; //sigma for observation error
  real<lower = 0> sd_q; //sigma on process error
  real<lower = 0> sd_site; //sigma on sites
  real<lower = 0> sd_dv; //sigma on divers
}

model{
  //priors
  x0 ~ normal(0,10);
  beta ~ normal(0,2);
  s~normal(0,5);
   
  //standard deviations
  sd_hab1 ~ cauchy(0,5);
  sd_hab2 ~ cauchy(0,5);
  sd_r1 ~ cauchy(0,5);
  sd_r2 ~ cauchy(0,5);
  sd_q ~ cauchy(0,5);
  sd_dv ~ cauchy(0,5);
  sd_site ~ cauchy(0,5);
  
  //varying intercepts
  a_hab1 ~ normal(0, sd_hab1);
  a_hab2 ~ normal(0, sd_hab2);
  a_site ~ normal(0, sd_site);
  a_dv ~ normal(0, sd_dv);
  
  //state process 
  x[1] ~ normal(x0,sd_q);
  
  for(t in 2:TT){
    x[t] ~ normal(x[t-1],sd_q);
  }
  
  //observation process for year
  
  for(i in 1:N_yr1){
    a_yr1[i]~ normal(x[yr_index1[i]],sd_r1); 
  }
  
  for(i in 1:N_yr2){
    a_yr2[i]~ normal(x[yr_index2[i]]+s,sd_r2); 
  }
 
  //draw presences
  for(i in 1:N_rvc){
    y1[i] ~ bernoulli_logit(a_yr1[year_id1[i]] + a_hab1[hab_class1[i]]);
  }
  
   for(n in 1:N_reef){
    y2[n] ~ bernoulli_logit(a_yr2[year_id2[n]] + a_hab2[hab_class2[n]] + a_site[site[n]]+ a_dv[diver[n]] + X[n,]*beta);
  }
}
"

logit_test_SS_sep<-"data{
  int TT; //timespan
  int N_rvc;//number of observations (SSU surveys)
  int y1[N_rvc]; //presence or absence on each survey
  int N_reef;//number of observations (REEF surveys)
  int y2[N_reef]; //presence or absence on each survey
  int N_hab1; //number of habitat classes
  int<lower=1,upper=N_hab1> hab_class1[N_rvc]; // vector of habitat class identities
  int N_hab2; //number of habitat classes
  int<lower=1,upper=N_hab2> hab_class2[N_reef]; // vector of habitat class identities
  int N_yr1; //number of years rvc
  int N_yr2; //number of years reef
  int yr_index1[N_yr1]; //index of years - RVC
  int<lower=1,upper=N_yr1> year_id1[N_rvc]; // vector of year ids
  int yr_index2[N_yr2]; //index of years - REEF
  int<lower=1,upper=N_yr2> year_id2[N_reef]; // vector of year ids
  int<lower=0> N_site; //number of sites in reef
  int<lower=1,upper=N_site> site[N_reef]; // vector of site identities for reef
  int<lower=0> N_dv; //number of divers in reef
  int<lower=1,upper=N_dv> diver[N_reef]; // vector of diver identities for reef
  int K; // columns in the covariate matrix
  matrix[N_reef,K] X; // design matrix X
}
parameters {
  //global intercept
  real x01;
  real x02;
  real x1[TT];
  real x2[TT];
  real beta;
  
  //deviations from intercept
  real a_hab1[N_hab1]; //deviation between habitats
  real a_hab2[N_hab2]; //deviation between habitats
  real a_yr1[N_yr1]; //deviation between years
  real a_yr2[N_yr2]; //deviation between years
  real a_site[N_site]; //deviation between sites
  real a_dv[N_dv]; //deviation between divers
  
  
  //st dev on the deviations
  real<lower = 0> sd_hab1; //sigma on habitat
  real<lower = 0> sd_hab2; //sigma on habitat
  real<lower = 0> sd_r1; //sigma for observation error
  real<lower = 0> sd_r2; //sigma for observation error
  real<lower = 0> sd_q1; //sigma on process error
  real<lower = 0> sd_q2; //sigma on process error
  real<lower = 0> sd_site; //sigma on sites
  real<lower = 0> sd_dv; //sigma on divers
}

model{
  //priors
  x01 ~ normal(0,10);
  x02 ~ normal(0,10);
  beta ~ normal(0,2);
   
  //standard deviations
  sd_hab1 ~ cauchy(0,5);
  sd_hab2 ~ cauchy(0,5);
  sd_r1 ~ cauchy(0,5);
  sd_r2 ~ cauchy(0,5);
  sd_q1 ~ cauchy(0,5);
  sd_q2 ~ cauchy(0,5);
  sd_dv ~ cauchy(0,5);
  sd_site ~ cauchy(0,5);
  
  //varying intercepts
  a_hab1 ~ normal(0, sd_hab1);
  a_hab2 ~ normal(0, sd_hab2);
  a_site ~ normal(0, sd_site);
  a_dv ~ normal(0, sd_dv);
  
  //state process 
  x1[1] ~ normal(x01,sd_q1);
  x2[1] ~ normal(x02,sd_q2);
  
  for(t in 2:TT){
    x1[t] ~ normal(x1[t-1],sd_q1);
    x2[t] ~ normal(x2[t-1],sd_q2);
  }
  
  //observation process for year
  
  for(i in 1:N_yr1){
    a_yr1[i]~ normal(x1[yr_index1[i]],sd_r1); 
  }
  
  for(i in 1:N_yr2){
    a_yr2[i]~ normal(x2[yr_index2[i]],sd_r2); 
  }
 
  //draw presences
  for(i in 1:N_rvc){
    y1[i] ~ bernoulli_logit(a_yr1[year_id1[i]] + a_hab1[hab_class1[i]]);
  }
  
   for(n in 1:N_reef){
    y2[n] ~ bernoulli_logit(a_yr2[year_id2[n]] + a_hab2[hab_class2[n]] + a_site[site[n]]+ a_dv[diver[n]] + X[n,]*beta);
  }
}
"


####Combined state - blue angelfish

blue_angel<- rvc_occs[[3]]

year_index<- data.frame(yr=seq(1993,2018),y.ind=seq(1,26))
blue_angel$year_index=year_index$y.ind[match(blue_angel$YEAR,year_index$yr)]
X_rvc<- matrix(data=c(scale(as.numeric(blue_angel$DEPTH))),ncol=1,nrow=nrow(blue_angel))


test_nb_rvc<- rstan::stan(model_code = nb_test_rvc, data = list(y = blue_angel$NUM.total, 
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


blue_angel<- rvc_occs[[1]]

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
                          pars = c('x',"sd_hab",'z_hab','sd_r','sd_q','z_yr','phi','beta'),
                          control = list(adapt_delta = 0.995,max_treedepth = 15), warmup = 200, chains = 4, iter = 1000, thin = 1)

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
               warmup = 400, iter = 1000, chains = 2)

make_stancode(ordered(abundance)~1+(1|geogr)+(1|hab_class2),
              data=reef_occs[[3]],
              family=cumulative("logit"))


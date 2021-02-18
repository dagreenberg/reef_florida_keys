rm(list=ls())
setwd("C:/Users/14388/Desktop/reef_florida_keys_data")
library(rstan);library(loo);library(bayesplot);library(dplyr);library(tidyverse);library(stringr);library(lubridate);library(rlist)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

####Functions####
rvc_occurrence = function(x,GZ,sp){
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
    rvc_occs[[i]]=x5
  }
  return(rvc_occs)
}


occ_ts_rvc = function(x){ #Takes the output from the previous function
  ts<- list()
  for(i in 1:length(x)){
    x1 = x[[i]]
    x2= x1 %>% group_by(YEAR) %>% summarize(n.occ=sum(occ),n.surv=n(),p.occ=n.occ/n.surv,sp=unique(SPECIES_CD))
    for(z in 1:nrow(x2)){
      if(x2$p.occ[z]==0){
        x2$p.occ[z]=NA
      }
    }
    ts[[i]]=x2
  }
  return(ts)
}


reef_occurrence = function(R,GZ,sp,geog){ #function to trim dataframe and add in occurrence data by species
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

occ_ts_reef = function(R,GZ,sp,geog){
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
  
  surveyors<- TempDat2 %>% group_by(fish_memberid) %>% summarize(n=length(unique(formid))) #Determine survey effort of each diver in the region
  surveyors_trim<- subset(surveyors,n>=5) #only keep surveys by members with 5 or more dives
  
  TempDat3<- subset(TempDat2, fish_memberid %in% surveyors_trim$fish_memberid) #Trim out surveys from surveys with less than 5
  site_subsets<-TempDat3 %>% group_by(geogr) %>% summarize(n=n_distinct(formid)) %>% subset(n>=5) #Calculate surveys per site
  
  TempDat4<- TempDat3 %>% subset(fish_memberid %in% surveyors_trim$fish_memberid) %>% subset(geogr %in% site_subsets$geogr) #Only keep well surveyed sites (n>=5)
  
  for(i in 1:nrow(sp)){
    occ<- subset(TempDat4,speciesid==sp$speciesid[i]) #Subset out each species in the provided dataframe
    occ_by_year<- occ %>% group_by(year) %>% summarize(n.occ=sum(occ),n.surv=n(),p.occ=n.occ/n.surv,sp=unique(sp$commonname[i]))
    occ_list[[i]]<- occ_by_year
  }
  return(occ_list)
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


rvc_occs<- rvc_occurrence(fk_93_18,GZ='3403',sp=fish_reef)
rvc_ts<- occ_ts_rvc(rvc_occs)
rvc_ts_filter<- rlist::list.filter(rvc_ts,length(na.omit(p.occ))>16)
gc()
reef_occs<- reef_occurrence(R,GZ='3403',sp=fish_reef,geog=reef_geog_3403)
reef_ts<- occ_ts_reef(R,GZ='3403',sp=fish_reef,geog=reef_geog_3403)


rep_survs<- R %>% group_by(geogr,month,day,year) %>% summarize(n_survs=n_distinct(formid))

####Stan model
logit_test_1<-"data{
  int<lower=1> N;//number of observations (SSU surveys)
  int y[N]; //presence or absence on each survey
  int<lower=0> N_hab; //number of habitat classes
  int<lower=1,upper=N_hab> hab_class[N]; // vector of habitat class identities
}
parameters {
  //global intercept
  real alpha;
  
  //
  
  //deviations from intercept
  real dev_hab[N_hab]; //deviation between psus
  
  //st dev on the deviations
  real<lower = 0> sigma_hab;   
}

transformed parameters{
  real alpha_h[N_hab];

  for(h in 1:N_hab){
    alpha_h[h] = alpha + dev_hab[h];
  }
}

model{
  //priors
  alpha ~ normal(0,10);
  
  //standard deviations
  sigma_hab ~ cauchy(0, 5);
  
  //varying intercepts
  dev_hab ~ normal(0, sigma_hab);
 
  for(i in 1:N){
    y[i] ~ bernoulli_logit(alpha_h[hab_class[i]]);
  }
}
"

logit_test_yr<-"data{
  int<lower=1> N;//number of observations (SSU surveys)
  int y[N]; //presence or absence on each survey
  int<lower=0> N_hab; //number of habitat classes
  int<lower=1,upper=N_hab> hab_class[N]; // vector of habitat class identities
  int<lower=0> N_mth; //number of months
  int<lower=1,upper=N_mth> month[N]; // vector of months
  int<lower=0> N_yr; //number of years
  int<lower=1,upper=N_yr> year[N]; // vector of year
}
parameters {
  //global intercept
  real alpha;
  
  //
  
  //deviations from intercept
  real a_hab[N_hab]; //deviation between habitats
  real a_mth[N_mth]; //deviation between months
  real a_yr[N_yr]; //deviation between months
  
  //st dev on the deviations
  real<lower = 0> sigma_hab;
  real<lower = 0> sigma_mth;
  real<lower = 0> sigma_yr;
}

transformed parameters{
  vector[N] eta;
  
  for(n in 1:N){
    eta[n] = alpha + a_hab[hab_class[n]] + a_mth[month[n]] + a_yr[year[n]];
  }
}

model{
  //priors
  alpha ~ normal(0,10);
  
  //standard deviations
  sigma_hab ~ cauchy(0, 5);
  sigma_mth ~ cauchy(0, 5);
  sigma_yr ~ cauchy(0, 5);
  
  //varying intercepts
  a_hab ~ normal(0, sigma_hab);
  a_mth ~ normal(0, sigma_hab);
  a_yr ~ student_t(5, 0, sigma_yr);
 
  for(i in 1:N){
    y[i] ~ bernoulli_logit(eta[i]);
  }
}
"

logit_test_rvc<-"data{
  int<lower=1> N;//number of observations (SSU surveys)
  int y[N]; //presence or absence on each survey
  int<lower=0> N_hab; //number of habitat classes
  int<lower=1,upper=N_hab> hab_class[N]; // vector of habitat class identities
  int<lower=0> N_yr; //number of years
  int<lower=1,upper=N_yr> year[N]; // vector of year
  int K; // columns in the covariate matrix
  matrix[N,K] X; // design matrix X
}
parameters {
  real alpha; 
  //global intercept
  
  //deviations from intercept
  real a_hab[N_hab]; //deviation between habitats
  real a_yr[N_yr]; //deviation between years
 
  
  //st dev on the deviations
  real<lower = 0> sigma_hab; //sigma on habitat
  real<lower = 0> sigma_yr; //sigma on year
}

transformed parameters{
 
  vector[N] eta;
  
  for(n in 1:N){
    eta[n] = alpha + a_yr[year[n]] + a_hab[hab_class[n]]+ X[n,]*beta;
  }
}

model{
  //priors
  alpha ~ normal(0,10);
  beta ~ normal(0,2)

  //standard deviations
  sigma_hab ~ cauchy(0, 5);
  sigma_yr ~ cauchy(0, 5);
  
  //varying intercepts
  a_hab ~ normal(0, sigma_hab);
  a_yr ~student_t(5, 0, sigma_yr);
  
  for(i in 1:N){
    y[i] ~ bernoulli_logit(eta[i]);
  }
}
"

logit_test_SS_rvc<-"data{
  int TT; //timespan
  int N;//number of observations (SSU surveys)
  int y[N]; //presence or absence on each survey
  int N_hab; //number of habitat classes
  int<lower=1,upper=N_hab> hab_class[N]; // vector of habitat class identities
  int N_yr; //number of years
  int yr_index[N_yr]; //index of years
  int<lower=1,upper=N_yr> year_id[N]; // vector of year ids
}
parameters {
  //global intercept
  real x0;
  real x[TT];
  
  //deviations from intercept
  real a_hab[N_hab]; //deviation between habitats
  real a_yr[N_yr]; //deviation between years
  
  //st dev on the deviations
  real<lower = 0> sd_hab; //sigma on habitat
  real<lower = 0> sd_r; //sigma for observation error
  real<lower = 0> sd_q; //sigma on process error
}

transformed parameters{
   vector[N] eta;
  
  for(n in 1:N){
    eta[n] = a_yr[year_id[n]] + a_hab[hab_class[n]];
  }
}  
  
model{
  //priors
  x0 ~ normal(0,10);
  
  //standard deviations
  sd_hab ~ cauchy(0,5);
  sd_r ~ cauchy(0,5);
  sd_q ~ cauchy(0,5);
  
  //varying intercepts
  a_hab ~ normal(0, sd_hab);
  
  //state process 
  x[1] ~ normal(x0,sd_q);
  
  for(t in 2:TT){
    x[t] ~ normal(x[t-1],sd_q);
  }
  
  //observation process for year
  for(i in 1:N_yr){
    a_yr[i]~ normal(x[yr_index[i]],sd_r); 
  }
 
  //draw presences
  for(i in 1:N){
    y[i] ~ bernoulli_logit(eta[i]);
  }
}
"

logit_test_reef<-"data{
  int<lower=1> N;//number of observations (SSU surveys)
  int y[N]; //presence or absence on each survey
  int<lower=0> N_hab; //number of habitat classes
  int<lower=1,upper=N_hab> hab_class[N]; // vector of habitat class identities
  int<lower=0> N_site; //number of sites
  int<lower=1,upper=N_site> site[N]; // vector of site identities
  int<lower=0> N_dv; //number of divers
  int<lower=1,upper=N_dv> diver[N]; // vector of diver identities
  int<lower=0> N_yr; //number of years
  int<lower=1,upper=N_yr> year[N]; // vector of year
  int K; // columns in the covariate matrix
  matrix[N,K] X; // design matrix X
  
}
parameters {
  //global intercept
  real alpha;
  
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
 
  for(i in 1:N){
    y[i] ~ bernoulli_logit(eta[i]);
  }
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

blue_angel<- rvc_occs[[1]]
blue_angel_ts<-rvc_ts[[1]]
blue_angel_reef<- reef_occs[[1]]

year_index<- data.frame(yr=seq(1993,2018),y.ind=seq(1,26))
blue_angel$year_index=year_index$y.ind[match(blue_angel$YEAR,year_index$yr)]



blue_angel_reef<- reef_occs[[1]]
blue_angel_reef<- subset(blue_angel_reef,year<=2018)
blue_angel_reef_ts<- reef_ts[[1]]

X<- matrix(data=c(scale(as.numeric(blue_angel_reef$btime)),scale(as.numeric(blue_angel_reef$averagedepth)),scale(as.numeric(blue_angel_reef$visibility)),scale(as.numeric(blue_angel_reef$current))),ncol=4,nrow=nrow(blue_angel_reef))

test_1_comb<- rstan::stan(model_code = logit_test_SS_comb, data = list(y2 = blue_angel_reef$occ, 
                                                                    N_reef = nrow(blue_angel_reef),
                                                                    N_hab2 = length(unique(blue_angel_reef$hab_class)),
                                                                    hab_class2 =as.numeric(factor(blue_angel_reef$hab_class)),
                                                                    N_yr2 =length(unique(blue_angel_reef$year)),
                                                                    year_id2=as.numeric(factor(blue_angel_reef$year)),
                                                                    yr_index2=sort(unique(as.numeric(factor(blue_angel_reef$year)))),
                                                                    site=as.numeric(factor(blue_angel_reef$geogr)),
                                                                    N_site=length(unique(blue_angel_reef$geogr)),
                                                                    diver=as.numeric(factor(blue_angel_reef$fish_memberid)),
                                                                    N_dv=length(unique(blue_angel_reef$fish_memberid)),
                                                                    K=ncol(X),
                                                                    X=X,
                                                                    y1 = blue_angel$occ, 
                                                                    N_rvc = nrow(blue_angel),
                                                                    N_hab1 = length(unique(blue_angel$HABITAT_CD)),
                                                                    hab_class1=as.numeric(factor(blue_angel$HABITAT_CD)),
                                                                    N_yr1 =length(unique(blue_angel$YEAR)),
                                                                    TT=26,
                                                                    year_id1=as.numeric(factor(blue_angel$YEAR)),
                                                                    yr_index1=sort(unique(blue_angel$year_index))),
                          pars = c("x",'s','a_hab1','a_hab2','a_yr1','a_yr2','sd_r1','sd_r2','sd_q'),
                          control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 500, chains = 4, iter = 2000, thin = 10)

posterior <- as.array(test_1_comb)
dimnames(posterior)
color_scheme_set("viridis")
mcmc_areas(
  posterior,
  pars = c(paste('x[',seq(1:26),']',sep='')),
  prob = 0.95, # 80% intervals
  prob_outer = 1, # 99%
  point_est = "median"
)
mcmc_dens_overlay(posterior, c("alpha", paste('a_yr[',seq(1:23),']',sep=''),'sigma_hab','sigma_yr'))
mcmc_pairs(posterior, pars = c(paste('x[',seq(1:26),']',sep=''),'sd_q','sd_r'),
           off_diag_args = list(size = 1.5))
mcmc_trace(posterior, pars = c("sd_q"))
mcmc_trace(posterior, pars = c("sd_r"))
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

#gray angelfish
gray_angel<- rvc_occs[[4]]
gray_angel_ts<-rvc_ts[[4]]
gray_angel_reef<- reef_occs[[4]]

year_index<- data.frame(yr=seq(1993,2018),y.ind=seq(1,26))
gray_angel$year_index=year_index$y.ind[match(gray_angel$YEAR,year_index$yr)]
gray_angel_reef<- reef_occs[[4]]
gray_angel_reef_ts<- reef_ts[[4]]

X<- matrix(data=c(scale(as.numeric(gray_angel_reef$btime)),scale(as.numeric(gray_angel_reef$visibility)),scale(as.numeric(gray_angel_reef$current))),ncol=3,nrow=nrow(gray_angel_reef))

test_1_comb_ga<- rstan::stan(model_code = logit_test_SS_comb, data = list(y2 = gray_angel_reef$occ, 
                                                                       N_reef = nrow(gray_angel_reef),
                                                                       N_hab2 = length(unique(gray_angel_reef$hab_class)),
                                                                       hab_class2 =as.numeric(factor(gray_angel_reef$hab_class)),
                                                                       N_yr2 =length(unique(gray_angel_reef$year)),
                                                                       year_id2=as.numeric(factor(gray_angel_reef$year)),
                                                                       yr_index2=sort(unique(as.numeric(factor(gray_angel_reef$year)))),
                                                                       site=as.numeric(factor(gray_angel_reef$geogr)),
                                                                       N_site=length(unique(gray_angel_reef$geogr)),
                                                                       diver=as.numeric(factor(gray_angel_reef$fish_memberid)),
                                                                       N_dv=length(unique(gray_angel_reef$fish_memberid)),
                                                                       K=ncol(X),
                                                                       X=X,
                                                                       y1 = gray_angel$occ, 
                                                                       N_rvc = nrow(gray_angel),
                                                                       N_hab1 = length(unique(gray_angel$HABITAT_CD)),
                                                                       hab_class1=as.numeric(factor(gray_angel$HABITAT_CD)),
                                                                       N_yr1 =length(unique(gray_angel$YEAR)),
                                                                       TT=26,
                                                                       year_id1=as.numeric(factor(gray_angel$YEAR)),
                                                                       yr_index1=sort(unique(gray_angel$year_index))),
                          pars = c("x",'s','a_hab1','a_hab2','a_yr1','a_yr2','sd_r1','sd_r2','sd_q','sd_dv','sd_hab1','sd_hab2','sd_site','beta'),
                          control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 500, chains = 4, iter = 2000, thin = 1)

posterior <- as.array(test_1_comb_ga)
dimnames(posterior)
color_scheme_set("viridis")
mcmc_areas(
  posterior,
  pars = c(paste('lp__',sep='')),
  prob = 0.95, # 80% intervals
  prob_outer = 1, # 99%
  point_est = "median"
)
mcmc_dens_overlay(posterior, c("alpha", paste('a_yr[',seq(1:23),']',sep=''),'sigma_hab','sigma_yr'))
mcmc_pairs(posterior, pars = c(paste('x[',seq(1:26),']',sep=''),'sd_q','sd_r'),
           off_diag_args = list(size = 1.5))
mcmc_trace(posterior, pars = c("sd_q"))
mcmc_trace(posterior, pars = c("sd_r"))
mcmc_trace(posterior, pars = c("sd_r"))


ga_params<- rstan::extract(test_1_comb_ga)
x_mat<- data.frame(median=NA,l.95=NA,u.95=NA)
for(i in 1:26){
  x_mat[i,1]=median(plogis(ga_params$x[,i]))
  x_mat[i,2]=quantile(plogis(ga_params$x[,i]),0.975)
  x_mat[i,3]=quantile(plogis(ga_params$x[,i]),0.025)
}
y_mat_1<- data.frame(median=NA,l.95=NA,u.95=NA)
for(i in 1:23){
  y_mat_1[i,1]=median(plogis(ga_params$a_yr1[,i]))
  y_mat_1[i,2]=quantile(plogis(ga_params$a_yr1[,i]),0.975)
  y_mat_1[i,3]=quantile(plogis(ga_params$a_yr1[,i]),0.025)
}

y_mat_2<- data.frame(median=NA,l.95=NA,u.95=NA)
for(i in 1:26){
  y_mat_2[i,1]=median(plogis(ga_params$a_yr2[,i]))
  y_mat_2[i,2]=quantile(plogis(ga_params$a_yr2[,i]+ga_params$s),0.975)
  y_mat_2[i,3]=quantile(plogis(ga_params$a_yr2[,i]+ga_params$s),0.025)
}

comp_plot_SS_comb(ts1=rvc_ts[[4]],ts2=reef_ts[[4]],x_mat = x_mat,y_mat1=y_mat_1,y_mat2=y_mat_2,sp='Gray Angelfish',GZ='Key Largo')


####Separate state - gray angelfish
test_1_sep_ga<- rstan::stan(model_code = logit_test_SS_sep, data = list(y2 = gray_angel_reef$occ, 
                                                                          N_reef = nrow(gray_angel_reef),
                                                                          N_hab2 = length(unique(gray_angel_reef$hab_class)),
                                                                          hab_class2 =as.numeric(factor(gray_angel_reef$hab_class)),
                                                                          N_yr2 =length(unique(gray_angel_reef$year)),
                                                                          year_id2=as.numeric(factor(gray_angel_reef$year)),
                                                                          yr_index2=sort(unique(as.numeric(factor(gray_angel_reef$year)))),
                                                                          site=as.numeric(factor(gray_angel_reef$geogr)),
                                                                          N_site=length(unique(gray_angel_reef$geogr)),
                                                                          diver=as.numeric(factor(gray_angel_reef$fish_memberid)),
                                                                          N_dv=length(unique(gray_angel_reef$fish_memberid)),
                                                                          K=ncol(X),
                                                                          X=X,
                                                                          y1 = gray_angel$occ, 
                                                                          N_rvc = nrow(gray_angel),
                                                                          N_hab1 = length(unique(gray_angel$HABITAT_CD)),
                                                                          hab_class1=as.numeric(factor(gray_angel$HABITAT_CD)),
                                                                          N_yr1 =length(unique(gray_angel$YEAR)),
                                                                          TT=26,
                                                                          year_id1=as.numeric(factor(gray_angel$YEAR)),
                                                                          yr_index1=sort(unique(gray_angel$year_index))),
                             pars = c("x1","x2",'a_hab1','a_hab2','a_yr1','a_yr2','sd_r1','sd_r2','sd_q1','sd_q2','sd_dv','sd_hab1','sd_hab2','sd_site','beta'),
                             control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 300, chains = 2, iter = 800, thin = 1)

posterior <- as.array(test_1_sep_ga)
dimnames(posterior)
color_scheme_set("viridis")
mcmc_areas(
  posterior,
  pars = c('lp__'),
  prob = 0.95, # 80% intervals
  prob_outer = 1, # 99%
  point_est = "median"
)
mcmc_dens_overlay(posterior, c("alpha", paste('a_yr[',seq(1:23),']',sep=''),'sigma_hab','sigma_yr'))
mcmc_pairs(posterior, pars = c(paste('x[',seq(1:26),']',sep=''),'sd_q','sd_r'),
           off_diag_args = list(size = 1.5))
mcmc_trace(posterior, pars = c("beta[1]"))
mcmc_trace(posterior, pars = c("sd_r"))
mcmc_trace(posterior, pars = c("sd_r"))


ga_params<- rstan::extract(test_1_sep_ga)
x_mat1<- data.frame(median=NA,l.95=NA,u.95=NA)
for(i in 1:26){
  x_mat1[i,1]=median(plogis(ga_params$x1[,i]))
  x_mat1[i,2]=quantile(plogis(ga_params$x1[,i]),0.975)
  x_mat1[i,3]=quantile(plogis(ga_params$x1[,i]),0.025)
}
x_mat2<- data.frame(median=NA,l.95=NA,u.95=NA)
for(i in 1:26){
  x_mat2[i,1]=median(plogis(ga_params$x2[,i]))
  x_mat2[i,2]=quantile(plogis(ga_params$x2[,i]),0.975)
  x_mat2[i,3]=quantile(plogis(ga_params$x2[,i]),0.025)
}
y_mat_1<- data.frame(median=NA,l.95=NA,u.95=NA)
for(i in 1:23){
  y_mat_1[i,1]=median(plogis(ga_params$a_yr1[,i]))
  y_mat_1[i,2]=quantile(plogis(ga_params$a_yr1[,i]),0.975)
  y_mat_1[i,3]=quantile(plogis(ga_params$a_yr1[,i]),0.025)
}

y_mat_2<- data.frame(median=NA,l.95=NA,u.95=NA)
for(i in 1:26){
  y_mat_2[i,1]=median(plogis(ga_params$a_yr2[,i]))
  y_mat_2[i,2]=quantile(plogis(ga_params$a_yr2[,i]+ga_params$s),0.975)
  y_mat_2[i,3]=quantile(plogis(ga_params$a_yr2[,i]+ga_params$s),0.025)
}

comp_plot_SS_sep(ts1=rvc_ts[[4]],ts2=reef_ts[[4]],x_mat1 = x_mat1,x_mat2 = x_mat2,y_mat1=y_mat_1,y_mat2=y_mat_2,sp='Gray Angelfish',GZ='Key Largo')


###REEF state space
gray_angel_reef<- reef_occs[[4]]
X<- matrix(data=c(scale(as.numeric(gray_angel_reef$btime)),scale(as.numeric(gray_angel_reef$averagedepth)),scale(as.numeric(gray_angel_reef$visibility)),scale(as.numeric(gray_angel_reef$current))),ncol=4,nrow=nrow(gray_angel_reef))

test_ga_reef<- rstan::stan(model_code = logit_test_SS_reef, data = list(y = gray_angel_reef$occ, 
                                                                    N = nrow(gray_angel_reef),
                                                                    N_hab = length(unique(gray_angel_reef$hab_class)),
                                                                    hab_class=as.numeric(factor(gray_angel_reef$hab_class)),
                                                                    N_yr =length(unique(gray_angel_reef$year)),
                                                                    year_id=as.numeric(factor(gray_angel_reef$year)),
                                                                    site=as.numeric(factor(gray_angel_reef$geogr)),
                                                                    N_site=length(unique(gray_angel_reef$geogr)),
                                                                    diver=as.numeric(factor(gray_angel_reef$fish_memberid)),
                                                                    N_dv=length(unique(gray_angel_reef$fish_memberid)),
                                                                    yr_index=sort(as.numeric(unique(factor(gray_angel_reef$year)))),
                                                                    TT=26,
                                                                    K=ncol(X),
                                                                    X=X),
                          pars = c("x", "a_hab",'a_yr','sigma_hab','sd_r','sd_q','sigma_dv','sigma_site','beta'),
                          control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 500, chains = 4, iter = 1000, thin = 1)

posterior <- as.array(test_ga_reef)
dimnames(posterior)
color_scheme_set("viridis")
mcmc_areas(
  posterior,
  pars = c(paste('beta[',seq(1:3),']',sep='')),
  prob = 0.95, # 80% intervals
  prob_outer = 1, # 99%
  point_est = "median"
)
mcmc_dens_overlay(posterior, c("alpha", paste('a_yr[',seq(1:23),']',sep=''),'sigma_hab','sigma_yr'))
mcmc_pairs(posterior, pars = c(paste('x[',seq(1:26),']',sep=''),'sd_q','sd_r'),
           off_diag_args = list(size = 1.5))
mcmc_trace(posterior, pars = c("beta[1]"))
mcmc_trace(posterior, pars = c("sd_r"))
mcmc_trace(posterior, pars = c("sd_r"))

ga_params<- rstan::extract(test_ga_reef)
x_mat<- data.frame(median=NA,l.95=NA,u.95=NA)
for(i in 1:26){
  x_mat[i,1]=median(plogis(ga_params$x[,i]))
  x_mat[i,2]=quantile(plogis(ga_params$x[,i]),0.975)
  x_mat[i,3]=quantile(plogis(ga_params$x[,i]),0.025)
}
y_mat<- data.frame(median=NA,l.95=NA,u.95=NA)
for(i in 1:23){
  y_mat[i,1]=median(plogis(ga_params$a_yr[,i]))
  y_mat[i,2]=quantile(plogis(ga_params$a_yr[,i]),0.975)
  y_mat[i,3]=quantile(plogis(ga_params$a_yr[,i]),0.025)
}

comp_plot_SS_reef(ts=reef_ts[[4]],x_mat = x_mat,y_mat=y_mat,sp='Gray Angelfish',GZ='Key Largo')


###RVC state space
test_ga_rvc<- rstan::stan(model_code = logit_test_SS_rvc, data = list(y = gray_angel$occ, 
                                                                      N = nrow(gray_angel),
                                                                      N_hab = length(unique(gray_angel$HABITAT_CD)),
                                                                      hab_class=as.numeric(factor(gray_angel$HABITAT_CD)),
                                                                      N_yr =length(unique(gray_angel$YEAR)),
                                                                      TT=26,
                                                                      year_id=as.numeric(factor(gray_angel$YEAR)),
                                                                      yr_index=sort(unique(gray_angel$year_index))),
                          pars = c("x", "a_hab",'a_yr','sd_q','sd_r'),
                          control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 400, chains = 4, iter = 2000, thin = 1)

ga_params<- rstan::extract(test_ga_rvc)
x_mat<- data.frame(median=NA,l.95=NA,u.95=NA)
for(i in 1:26){
  x_mat[i,1]=median(plogis(ga_params$x[,i]))
  x_mat[i,2]=quantile(plogis(ga_params$x[,i]),0.975)
  x_mat[i,3]=quantile(plogis(ga_params$x[,i]),0.025)
}
y_mat<- data.frame(median=NA,l.95=NA,u.95=NA)
for(i in 1:23){
  y_mat[i,1]=median(plogis(ga_params$a_yr[,i]))
  y_mat[i,2]=quantile(plogis(ga_params$a_yr[,i]),0.975)
  y_mat[i,3]=quantile(plogis(ga_params$a_yr[,i]),0.025)
}

comp_plot_SS_rvc(ts=rvc_ts[[4]],x_mat = x_mat,y_mat=y_mat,sp='Gray Angelfish',GZ='Key Largo')

posterior <- as.array(test_1_y_g)
dimnames(posterior)
color_scheme_set("viridis")
mcmc_areas(
  posterior,
  pars = c(paste('a_hab[',seq(1:9),']',sep=''),'sd_q','sd_r'),
  prob = 0.95, # 80% intervals
  prob_outer = 1, # 99%
  point_est = "median"
)
mcmc_dens_overlay(posterior, c(paste('a_yr[',seq(1:23),']',sep=''),'sd_r','sd_q'))
mcmc_pairs(posterior, pars = c(paste('x[',seq(1:5),']',sep=''),'sd_q','sd_r'),
           off_diag_args = list(size = 1.5))
mcmc_trace(posterior, pars = c("sd_q"))

comp_plot(ts=rvc_ts[[4]],mod_mat = a_mat,sp='Gray Angelfish',GZ='Key Largo')


ypos <- Y[!is.na(Y)]
n<- n[!is.na(n)]
n_pos <- length(ypos)  #number on non-NA ys
indx_pos <- which(!is.na(Y), arr.ind = TRUE)  #index on the non-NAs
col_indx_pos <- as.vector(indx_pos[, "col"])
row_indx_pos <- as.vector(indx_pos[, "row"])
ts_pos<- ifelse(row_indx_pos==1,0,1)
ts_plot(ts=Y,sp='')

test_1_hab<- rstan::stan(model_code = logit_test_hab, data = list(y = blue_angel$occ, 
                                                          N = nrow(blue_angel),N_hab = length(unique(blue_angel$HABITAT_CD)),hab_class=as.numeric(factor(blue_angel$HABITAT_CD)),N_psu = length(unique(blue_angel$psu_id)),psu_id=blue_angel$psu_id), 
                     pars = c("alpha", "a_hab",'sigma_hab'),
                     control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 400, chains = 4, iter = 2000, thin = 1)

test_1_psu<- rstan::stan(model_code = logit_test_psu, data = list(y = blue_angel$occ, 
                                                          N = nrow(blue_angel),N_hab = length(unique(blue_angel$HABITAT_CD)),hab_class=as.numeric(factor(blue_angel$HABITAT_CD)),N_psu = length(unique(blue_angel$psu_id)),psu_id=blue_angel$psu_id), 
                     pars = c("alpha", "a_hab",'sigma_hab','sigma_psu'),
                     control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 400, chains = 4, iter = 2000, thin = 1)


posterior <- as.array(test_1)
dimnames(posterior)
color_scheme_set("viridis")
mcmc_areas(
  posterior,
  pars = c("alpha", "sigma_hab"),
  prob = 0.95, # 80% intervals
  prob_outer = 1, # 99%
  point_est = "median"
)
mcmc_dens_overlay(posterior, pars = c("alpha", "sigma_hab",'alpha_h[1]','alpha_h[2]','alpha_h[3]','alpha_h[4]'))
mcmc_pairs(posterior, pars = c("alpha", "sigma_hab",'alpha_h[1]','alpha_h[2]','alpha_h[3]','alpha_h[4]'),
           off_diag_args = list(size = 1.5))
mcmc_trace(posterior, pars = c("sd_q"))







test_1_y<- rstan::stan(model_code = logit_test, data = list(y = blue_angel$occ, 
                                                            N = nrow(blue_angel),N_psu = length(unique(blue_angel$psu_id)),N_hab = length(unique(blue_angel$HABITAT_CD)),psu_in_hab=psu_in_hab,hab_class=as.numeric(factor(blue_angel$HABITAT_CD)),psu_id=blue_angel$psu_id), 
                       pars = c("alpha", "alpha_h","alpha_p",'sigma_hab','sigma_psu'),
                       control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 400, chains = 4, iter = 2000, thin = 1)

test_1_y<- rstan::stan(model_code = logit_test_SS, data = list(y = blue_angel$occ, 
                                                     N = nrow(blue_angel),
                                                     N_hab = length(unique(blue_angel$HABITAT_CD)),
                                                     hab_class=as.numeric(factor(blue_angel$HABITAT_CD)),
                                                     N_yr =length(unique(blue_angel$YEAR)),
                                                     TT=26,
                                                     year_id=as.numeric(factor(blue_angel$YEAR)),
                                                     yr_index=sort(unique(blue_angel$year_index))),
                                                    pars = c("x", "a_hab",'a_yr','sd_q','sd_r'),
                                                  control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 400, chains = 4, iter = 3000, thin = 10)

ba_params<- rstan::extract(test_1_y)
x_mat<- data.frame(median=NA,l.95=NA,u.95=NA)
for(i in 1:26){
  x_mat[i,1]=median(plogis(ba_params$x[,i]))
  x_mat[i,2]=quantile(plogis(ba_params$x[,i]),0.975)
  x_mat[i,3]=quantile(plogis(ba_params$x[,i]),0.025)
}
y_mat<- data.frame(median=NA,l.95=NA,u.95=NA)
for(i in 1:23){
  y_mat[i,1]=median(plogis(ba_params$a_yr[,i]))
  y_mat[i,2]=quantile(plogis(ba_params$a_yr[,i]),0.975)
  y_mat[i,3]=quantile(plogis(ba_params$a_yr[,i]),0.025)
}

comp_plot_SS_rvc(ts=rvc_ts[[1]],x_mat = x_mat,y_mat=y_mat,sp='Blue Angelfish',GZ='Key Largo')


posterior <- as.array(test_1_y)
dimnames(posterior)
color_scheme_set("viridis")
mcmc_areas(
  posterior,
  pars = c(paste('x[',seq(1:26),']',sep=''),'sd_q','sd_r'),
  prob = 0.95, # 80% intervals
  prob_outer = 1, # 99%
  point_est = "median"
)
mcmc_dens_overlay(posterior, c("alpha", paste('a_yr[',seq(1:23),']',sep=''),'sigma_hab','sigma_yr'))
mcmc_pairs(posterior, pars = c(paste('x[',seq(1:26),']',sep=''),'sd_q','sd_r'),
           off_diag_args = list(size = 1.5))
mcmc_trace(posterior, pars = c("sd_q"))
mcmc_trace(posterior, pars = c("sd_r"))
mcmc_trace(posterior, pars = c("sd_r"))



test_1_y2<- rstan::stan(model_code = logit_test_yr2, data = list(y = blue_angel$occ, 
                                                               N = nrow(blue_angel),
                                                               N_hab = length(unique(blue_angel$HABITAT_CD)),
                                                               hab_class=as.numeric(factor(blue_angel$HABITAT_CD)),
                                                               N_yr =length(unique(blue_angel$YEAR)),
                                                               year=as.numeric(factor(blue_angel$YEAR))),
                       pars = c("alpha", "a_hab",'sigma_hab','sigma_yr'),
                       control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 400, chains = 4, iter = 2000, thin = 1)

test_1_y2<- rstan::stan(model_code = logit_test_yr2, data = list(y = blue_angel$occ, 
                                                                 N = nrow(blue_angel),
                                                                 N_hab = length(unique(blue_angel$HABITAT_CD)),
                                                                 hab_class=as.numeric(factor(blue_angel$HABITAT_CD)),
                                                                 N_yr =length(unique(blue_angel$YEAR)),
                                                                 year=as.numeric(factor(blue_angel$YEAR))),
                        pars = c("alpha", "a_hab",'a_yr','sigma_hab','sigma_yr'),
                        control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 400, chains = 4, iter = 2000, thin = 1)

ba_params<- rstan::extract(test_1_y2)
a_mat<- data.frame(median=NA,l.95=NA,u.95=NA)
for(i in 1:23){
  a_mat[i,1]=median(plogis(ba_params$alpha+ba_params$a_yr[,i]))
  a_mat[i,2]=quantile(plogis(ba_params$alpha+ba_params$a_yr[,i]),0.975)
  a_mat[i,3]=quantile(plogis(ba_params$alpha+ba_params$a_yr[,i]),0.025)
}

comp_plot_rvc(ts=rvc_ts[[1]],mod_mat = a_mat,sp='Blue Angelfish',GZ='Key Largo')

gray_angel<- rvc_occs[[4]] #gray angelfish in key largo
year_index<- data.frame(yr=seq(1993,2018),y.ind=seq(1,26))
gray_angel$year_index=year_index$y.ind[match(gray_angel$YEAR,year_index$yr)]

test_1_y2_g<- rstan::stan(model_code = logit_test_yr2, data = list(y = gray_angel$occ, 
                                                                 N = nrow(gray_angel),
                                                                 N_hab = length(unique(gray_angel$HABITAT_CD)),
                                                                 hab_class=as.numeric(factor(gray_angel$HABITAT_CD)),
                                                                 N_yr =length(unique(gray_angel$YEAR)),
                                                                 year=as.numeric(factor(gray_angel$YEAR))),
                        pars = c("alpha", "a_hab",'a_yr','sigma_hab','sigma_yr'),
                        control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 400, chains = 4, iter = 2000, thin = 1)

ga_params<- rstan::extract(test_1_y2_g)
a_mat<- data.frame(median=NA,l.95=NA,u.95=NA)
for(i in 1:23){
  a_mat[i,1]=median(plogis(ga_params$alpha+ga_params$a_yr[,i]))
  a_mat[i,2]=quantile(plogis(ga_params$alpha+ga_params$a_yr[,i]),0.975)
  a_mat[i,3]=quantile(plogis(ga_params$alpha+ga_params$a_yr[,i]),0.025)
}

comp_plot_rvc(ts=rvc_ts[[4]],mod_mat = a_mat,sp='Gray Angelfish',GZ='Key Largo')


test_ga_rvc<- rstan::stan(model_code = logit_test_SS_rvc, data = list(y = gray_angel$occ, 
                                                               N = nrow(gray_angel),
                                                               N_hab = length(unique(gray_angel$HABITAT_CD)),
                                                               hab_class=as.numeric(factor(gray_angel$HABITAT_CD)),
                                                               N_yr =length(unique(gray_angel$YEAR)),
                                                               TT=26,
                                                               year_id=as.numeric(factor(gray_angel$YEAR)),
                                                               yr_index=sort(unique(gray_angel$year_index))),
                       pars = c("x", "a_hab",'a_yr','sd_q','sd_r'),
                       control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 400, chains = 4, iter = 2000, thin = 1)

ga_params<- rstan::extract(test_1_y_g)
x_mat<- data.frame(median=NA,l.95=NA,u.95=NA)
for(i in 1:26){
  x_mat[i,1]=median(plogis(ga_params$x[,i]))
  x_mat[i,2]=quantile(plogis(ga_params$x[,i]),0.975)
  x_mat[i,3]=quantile(plogis(ga_params$x[,i]),0.025)
}
y_mat<- data.frame(median=NA,l.95=NA,u.95=NA)
for(i in 1:23){
  y_mat[i,1]=median(plogis(ga_params$a_yr[,i]))
  y_mat[i,2]=quantile(plogis(ga_params$a_yr[,i]),0.975)
  y_mat[i,3]=quantile(plogis(ga_params$a_yr[,i]),0.025)
}

comp_plot_SS_rvc(ts=rvc_ts[[4]],x_mat = x_mat,y_mat=y_mat,sp='Gray Angelfish',GZ='Key Largo')

posterior <- as.array(test_1_y_g)
dimnames(posterior)
color_scheme_set("viridis")
mcmc_areas(
  posterior,
  pars = c(paste('a_hab[',seq(1:9),']',sep=''),'sd_q','sd_r'),
  prob = 0.95, # 80% intervals
  prob_outer = 1, # 99%
  point_est = "median"
)
mcmc_dens_overlay(posterior, c(paste('a_yr[',seq(1:23),']',sep=''),'sd_r','sd_q'))
mcmc_pairs(posterior, pars = c(paste('x[',seq(1:5),']',sep=''),'sd_q','sd_r'),
           off_diag_args = list(size = 1.5))
mcmc_trace(posterior, pars = c("sd_q"))




comp_plot(ts=rvc_ts[[4]],mod_mat = a_mat,sp='Gray Angelfish',GZ='Key Largo')

queen_angel<- rvc_occs[[5]] #queen angelfish in key largo
test_1_y2_g2<- rstan::stan(model_code = logit_test_yr2, data = list(y = queen_angel$occ, 
                                                                   N = nrow(queen_angel),
                                                                   N_hab = length(unique(queen_angel$HABITAT_CD)),
                                                                   hab_class=as.numeric(factor(queen_angel$HABITAT_CD)),
                                                                   N_yr =length(unique(queen_angel$YEAR)),
                                                                   year=as.numeric(factor(queen_angel$YEAR))),
                          pars = c("alpha", "a_hab",'a_yr','sigma_hab','sigma_yr'),
                          control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 400, chains = 4, iter = 2000, thin = 1)

ga_params<- rstan::extract(test_1_y2_g2)
a_mat<- data.frame(median=NA,l.95=NA,u.95=NA)
for(i in 1:23){
  a_mat[i,1]=median(plogis(ga_params$alpha+ga_params$a_yr[,i]))
  a_mat[i,2]=quantile(plogis(ga_params$alpha+ga_params$a_yr[,i]),0.975)
  a_mat[i,3]=quantile(plogis(ga_params$alpha+ga_params$a_yr[,i]),0.025)
}

comp_plot(ts=rvc_ts[[5]],mod_mat = a_mat,sp='Gray Angelfish',GZ='Key West')

posterior <- as.array(test_1_y2_g)
dimnames(posterior)
color_scheme_set("viridis")
mcmc_areas(
  posterior,
  pars = c("alpha", paste('a_hab[',seq(1:9),']',sep=''),'sigma_hab','sigma_yr'),
  prob = 0.95, # 80% intervals
  prob_outer = 1, # 99%
  point_est = "median"
)
mcmc_dens_overlay(posterior, c("alpha", paste('a_yr[',seq(1:23),']',sep=''),'sigma_hab','sigma_yr'))
mcmc_pairs(posterior, pars = c("alpha", paste('a_yr[',seq(1:23),']',sep=''),'sigma_hab','sigma_yr'),
           off_diag_args = list(size = 1.5))
mcmc_trace(posterior, pars = c("sd_q"))


###test for reef model
gray_angel_reef<- reef_occs[[4]]
X<- matrix(data=c(scale(as.numeric(gray_angel_reef$btime)),scale(as.numeric(gray_angel_reef$visibility)),scale(as.numeric(gray_angel_reef$current))),ncol=3,nrow=nrow(gray_angel_reef))

test_1_reef<- rstan::stan(model_code = logit_test_reef, data = list(y = gray_angel_reef$occ, 
                                                                    N = nrow(gray_angel_reef),
                                                                    N_hab = length(unique(gray_angel_reef$hab_class)),
                                                                    hab_class=as.numeric(factor(gray_angel_reef$hab_class)),
                                                                    N_yr =length(unique(gray_angel_reef$year)),
                                                                    year=as.numeric(factor(gray_angel_reef$year)),
                                                                    site=as.numeric(factor(gray_angel_reef$geogr)),
                                                                    N_site=length(unique(gray_angel_reef$geogr)),
                                                                    diver=as.numeric(factor(gray_angel_reef$fish_memberid)),
                                                                    N_dv=length(unique(gray_angel_reef$fish_memberid)),
                                                                    K=ncol(X),
                                                                    X=X),
                           pars = c("alpha", "a_hab",'a_yr','sigma_hab','sigma_yr','sigma_dv','sigma_site','beta'),
                           control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 100, chains = 4, iter = 500, thin = 1)

ba_params<- rstan::extract(test_1_reef)
a_mat<- data.frame(median=NA,l.95=NA,u.95=NA)
for(i in 1:23){
  a_mat[i,1]=median(plogis(ba_params$alpha+ba_params$a_yr[,i]))
  a_mat[i,2]=quantile(plogis(ba_params$alpha+ba_params$a_yr[,i]),0.975)
  a_mat[i,3]=quantile(plogis(ba_params$alpha+ba_params$a_yr[,i]),0.025)
}

comp_plot_reef(ts=occs_ts_reef[[4]][1:23,],mod_mat = a_mat,sp='Blue Angelfish',GZ='Key Largo')

posterior <- as.array(test_1_reef)
dimnames(posterior)
color_scheme_set("viridis")
mcmc_areas(
  posterior,
  pars = c("alpha", paste('a_hab[',seq(1:5),']',sep=''),'sigma_hab','sigma_yr'),
  prob = 0.95, # 80% intervals
  prob_outer = 1, # 99%
  point_est = "median"
)
mcmc_dens_overlay(posterior, c("alpha", paste('a_yr[',seq(1:23),']',sep=''),'sigma_hab','sigma_yr'))
mcmc_pairs(posterior, pars = c("alpha", paste('a_yr[',seq(1:23),']',sep=''),'sigma_hab','sigma_yr'),
           off_diag_args = list(size = 1.5))
mcmc_trace(posterior, pars = c("sd_q"))

mod_1<-"data{
  int<lower=1> N;//number of observations (SSU surveys)
  vector y[N]; //presence or absence on each survey
  int<lower=0> N_psu; //number of primary sample units
  int<lower=0> N_hab; //number of habitat classes
  int<lower=0> N_month; //number of months
  int<lower=0> N_years;//number of years of survey
  int<lower=1,upper=N_psu> psu_id[N] // vector of psu identities
  int<lower=1,upper=h>    psu_in_hab[N_psu] // index of psu in each habitat class
  int<lower=1,upper=N_hab> hab_class[N] // vector of habitat class identities
  int<lower=1,upper=N_years> year[N] // vector of year identities

}
parameters {
  real alpha; // mean log-odds (logit)
  real dev_psu[N_psu]; //deviation between psus
  real dev_hab[N_hab]; //deviation between psus
  real dev_mth[N_psu]; //deviation between psus
  real dev_year[N_hab]; //deviation between psus
  real<lower = 0> sigma_psu;
  real<lower = 0> sigma_hab;   
  real alpha;
  real beta;
  vector[N_subj] u;
}
transformed parameters{
  for(h in 1:N_hab){
    alpha_h[h] = alpha + dev_hab[h];
  }
  for(p in 1:N_psu){
    alpha_p[p] = mu_hab[psu_in_hab[p]] + dev_psu[p];
  }
}
model{
 
  for(i in 1:N){
    y[i] ~ bernoulli_logit(alpha_p[psu_id[i]]);
  }
}
"


blue_angelfish<- subset(fk_93_18,SPECIES_CD=='HOL BERM')

ts_comp[[1]]

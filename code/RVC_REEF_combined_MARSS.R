#RVC data & MARSS analysis
rm(list=ls())
setwd("C:/Users/14388/Desktop/reef_florida_keys_data")
library(dplyr);library(magrittr);library(tidyverse); library(MARSS)
library(stringr);library(lubridate)

####1. Functions ####
#not in filter
`%notin%` <- Negate(`%in%`)
#geometric mean
gm_mean<- function(x){
  prod(x)^(1/length(x))
}

#Aggregate by SSU (15m, 15min diameter survey;some samples are broken down to multiple rows per SSU survey)
ssu_density = function(x){
  x %>% dplyr::group_by(STATION_NR,PRIMARY_SAMPLE_UNIT,YEAR) %>%
    dplyr::summarise(density=sum(na.omit(NUM)),occ=NA) %>% #Sums up the number of counts
    mutate(occ=ifelse(density>0,1,0)) %>% arrange(PRIMARY_SAMPLE_UNIT,YEAR,STATION_NR) #Also scores presence/absence at the SSU level
}


#Aggregate by PSU (200x200m primary sample location; 1-4 SSUs per PSU)
#Calculates number of SSUs, the total counts across the PSU, average count per SSU, and occupancy at the PSU level
psu_density = function(x) {
  ssu_density(x) %>% dplyr::group_by(PRIMARY_SAMPLE_UNIT,YEAR) %>%
    summarise(m=length(STATION_NR),ssu.var=var(density),mean.ssu.abundance=mean(density),occ=NA) %>% 
    mutate(occ=ifelse(mean.ssu.abundance>0,1,0)) %>% arrange(PRIMARY_SAMPLE_UNIT,YEAR)
}

#Transform PSU data to relative abundance and back-transform to the min. threshold
psu_relative_abund = function(x) {
  x %>% ssu_density() %>% dplyr::group_by(PRIMARY_SAMPLE_UNIT,YEAR) %>%
    summarise(m=length(STATION_NR),var=var(na.omit(density)),abundance=sum(na.omit(density)),RA=NA,abund_trans=NA) %>%  
    mutate(
      RA=ifelse(abundance==0,0,RA),
      RA=ifelse(abundance>0&abundance<=1,1,RA),
      RA=ifelse(abundance>1&abundance<=10,2,RA),
      RA=ifelse(abundance>10&abundance<=100,3,RA),
      RA=ifelse(abundance>100,4,RA)) %>%
    mutate(
      abund_trans=ifelse(RA==0,0,abund_trans),
      abund_trans=ifelse(RA==1,1,abund_trans),
      abund_trans=ifelse(RA==2,2,abund_trans),
      abund_trans=ifelse(RA==3,11,abund_trans),
      abund_trans=ifelse(RA==4,101,abund_trans))
}

#Create a yearly aggregate of density/total abundance and fill in missing years
psu_density_year = function(x) {
  x %>% psu_density() %>% 
    group_by(YEAR) %>% summarise(n=length(PRIMARY_SAMPLE_UNIT), nm=sum(m), mean_ssu_abundance=mean(mean.ssu.abundance),var_abundance=var(mean.ssu.abundance),n.occ=sum(occ),p.occ=sum(occ)/n()) %>%
    complete(YEAR=seq(min(na.omit(x$YEAR)),max(na.omit(x$YEAR)),by=1))
}

#Create a yearly aggregate of relative abundance counts and fill in missing years
psu_relative_abund_year = function(x) {x %>% psu_relative_abund() %>%
    group_by(YEAR,RA) %>%  summarise(n=n()) %>% spread(key=RA, n, fill = 0, sep = ".") %>% as.data.frame() %>%
    complete(YEAR=seq(min(x$YEAR),max(x$YEAR),by=1))
}

#Create a yearly aggregate of relative and transformed abundance and fill in missing years
psu_transform_year = function(x) {x %>% psu_relative_abund() %>%
    group_by(YEAR) %>% summarise(mean.RA.cat=mean(RA),mean.abund.trans=mean(abund_trans),var.abund.trans=var(abund_trans)) %>%
    complete(YEAR=seq(min(x$YEAR),max(x$YEAR),by=1))
}

rvc_batch = function(x,GZ,sp,geog){
  x$SSU_YEAR<- paste(x$PRIMARY_SAMPLE_UNIT,x$STATION_NR,x$YEAR,sep='_')
  x1= x %>% subset(region.id==GZ) %>% subset(LAT_LON %in% geog$lat_lon) #Trim down full dataset to those in the selected plots (eg. only spur and groove) 
  x1= complete(x1,SSU_YEAR,nesting(SPECIES_CD),fill=list(NUM=0)) #This makes sure the zeros in each sample unit are recorded for each species
  
  RVC_TS=list()
  for(i in 1:nrow(sp)){
    x2= x1 %>% subset(SPECIES_CD==sp$rvc_code[i]) %>% psu_density_year() #For each species calculate annual densities in the SSUs
    
    for(z in 1:nrow(x2)){ #Filter out years that lack enough data
      if(is.na(x2$nm[z])==T){next}
      if(x2$nm[z]<20){ #less than 20 sample units into NAs
        x2[z,4:7]<- NA 
      }
      if(is.na(x2$mean_ssu_abundance[z])==T){next}
      if(x2$mean_ssu_abundance[z]==0){ #zeros turned into NAs
        x2[z,4:7]<- NA 
      }
     
    }
    x_full<- cbind(x2,Species=rep(unique(sp$commonname)[i],nrow(x2)))
    RVC_TS[[i]]=x_full
  }
  return(RVC_TS)
}

# Create function to calculate spp specific abundance scores by year for a given region
REEF_pull <- function(R,GZ,sp,geog){
  # R is the REEF data in flat format (raw format provided by REEF)
  # GZ is the geozone you are interested in summarizing (can be any length)
  # sp is the matrix of sp you have previously used to subset the R raw data
  # geog are the dive sites (lat/lons) you want to include for the time-series
  Time_series<- list()
  #pull only those data (rows) from the geozone given to function
  TempDat<-R %>% subset(site4==GZ) %>% subset(site %in% geog$geogid) %>% select('formid','speciesid','abundance',everything())
  Zeros<- complete(TempDat,formid,nesting(speciesid),
                   fill=list(abundance=0)) %>%
    anti_join(TempDat) %>% 
    select('formid','speciesid','abundance',everything())
  
  m<- match(Zeros$formid,TempDat$formid) #Match the zero data frame to the rest of the data
  Zeros[,4:ncol(Zeros)]<- TempDat[m,4:ncol(TempDat)] #Replicate the survey-level data (except abundance)
  TempDat2<- rbind(TempDat,Zeros) %>% arrange(formid,speciesid) #Combine presence and absence data for all species
  TempDat2$abund_trans <-NA #Transform abundance categories into the minimum counts
  TempDat2<- TempDat2%>%
    mutate(
      abund_trans=ifelse(abundance==0,0,abund_trans),
      abund_trans=ifelse(abundance==1,1,abund_trans),
      abund_trans=ifelse(abundance==2,2,abund_trans),
      abund_trans=ifelse(abundance==3,11,abund_trans),
      abund_trans=ifelse(abundance==4,101,abund_trans)
    )
  
  for(i in 1:nrow(sp)){
    Temp_sp<- subset(TempDat2,speciesid==sp$speciesid[i] & abundance>0) #Extract out all occurrences of species
    if(nrow(Temp_sp)==0){next} 
    TempDat_sites<- subset(TempDat2, site %in% Temp_sp$site & speciesid==sp$speciesid[i]) #Now include all surveys at these sites
    
    Occs<- as.data.frame(Temp_sp %>% #Summarize for all surveys where the target species was found
                              group_by(site,year) %>%
                               group_by(site,year)%>% 
                               summarize(mean_abund=mean(abund_trans),n=n_distinct(formid))  %>%
                              group_by(year) %>%
                              summarize(mean_site_abund.occ=mean(mean_abund),n.sites.occ=n_distinct(site),n.survs.occ=sum(n)))
                  
  
    Survs<- as.data.frame(TempDat_sites %>% #Summarize for all surveys
                            group_by(site,year)%>% 
                            summarize(mean_abund=mean(abund_trans),n=n_distinct(formid)) %>% 
                            group_by(year) %>%
                            summarize(mean_site_abund=mean(mean_abund),sd_site_abund=sd(mean_abund),n.sites=n_distinct(site),n.survs=sum(n)))
    
    TS<- full_join(Survs,Occs,by='year')
    TS$comName<- rep(sp$commonname[i],nrow(TS))
    for(z in 1:nrow(TS)){
      if(TS$n.survs[z]<20){ #Make years with less than 20 surveys NAs
        TS$mean_site_abund[z]=NA
      }
      if(TS$n.sites[z]<5){ #Make years with less than 5 sites surveyed into NAs
        TS$mean_site_abund[z]=NA
      }
      if(TS$mean_site_abund[z]==0){ #turn zeros into NAs
        TS$mean_site_abund[z]=NA
      }
    }
    
    Time_series[[i]]<- TS
  }
  
  # calculate the sum of reported abundance codes for each spp by year (when reported for a year!!)
  Time_series<- plyr::compact(Time_series)
  return(Time_series)
} #end function

#Plot function


#####2. Data Loading and site matching####
###REEF data
REEF<- read.csv("./REEF Tropical Western Atlantic/TWA.csv") #Reef sightings 
REEF_survey<- read.csv("./REEF Tropical Western Atlantic/TWAsurveys.csv") #Survey metadata
m<- match(REEF$formid,REEF_survey$formid) #match sightings to surveys
REEF[,14:30]<- REEF_survey[m,4:20] #merge survey level data

R<- REEF %>% group_by(formid) %>% summarize(n=n(),sum.abundance=sum(abundance)) 
R_filter<- subset(R,n==sum.abundance) #Some surveys have all 1s for every species - may be recording occurrence only, removing these

#Remove these surveys
R<- subset(REEF, formid %notin% R_filter$formid)

# Get rid of  rows of data with no date reported
R<-R[-which(R$date=="0000-00-00"),]

# let's get the year & month from the Date field (requires lubridate)
R$date<-ymd(R$date)#put into proper date format
R<-cbind(R,year=year(R$date))
R<-cbind(R,month=month(R$date))

###REEF geog
library(sp)
reef_geog<- read.csv("./REEF Tropical Western Atlantic/TWAgeog.csv")
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

# Remove sites without lat-long (b/c those sites are probably janky anyway)
R<-filter(R, R$geogr %in% reef_geog$geogid)

# Thin  raw data to exclude dives shorter than 20min, longer than 120min
R<-R[as.numeric(R$btime)>20,]
R<-R[as.numeric(R$btime)<120,]

# Thin the data to remove night dives (start before 5am or after 8pm)
R<-R[as.numeric(R$start)>5,]
R<-R[as.numeric(R$start)<20,]


####Florida keys RVC data
fk_99_18<- read.csv("./RVC/Florida_Keys_RVC.csv") #RVC survey data from 1999 to 2018
fk_79_98<- read.csv("./RVC/Florida_keys_pre_1999.txt") #RVC survey dating from 1979 to 1998

fk_79_98$len<- as.numeric(ifelse(fk_79_98$len==-9,0,fk_79_98$len)) #Replace -9, these represent 0s

t<- fk_79_98 %>% group_by(station_nr,psu,YEAR,spcode) %>%  #For each SSU and species
  summarise(NUM=n(),sum=sum(len),id=ID[1]) %>% arrange(psu,YEAR,station_nr,spcode) #Summarize the number of individual entries (non-zero lengths indicate sightings)
t$NUM<- ifelse(t$sum==0,0,t$NUM) #NUM also counts zeros, now only presences

#Put back NUM into original data-frame, now each row represents a SSU survey
m<- match(fk_79_98$ID,t$id) #match back the sums to the survey id 
fk_79_98$NUM<- t$NUM[m] #Create the number column in this survey set
fk_79_98<- fk_79_98[complete.cases(fk_79_98$NUM),] #Reduce down to the summed count only

#Synonymize column names
colnames(fk_79_98)<- c(colnames(fk_99_18)[1],colnames(fk_99_18)[3:5],colnames(fk_99_18)[2],colnames(fk_99_18)[6],colnames(fk_99_18)[7:18],colnames(fk_99_18)[20],colnames(fk_99_18)[19]) 
#Join together datasets
fk_79_18<- full_join(fk_99_18,fk_79_98)

#Extract out florida key subregions for select zones
reef_3403_sites<- subset(reef_geog,region.id=='3403')
reef_3404_sites<- subset(reef_geog,region.id=='3404')
reef_3408_sites<- subset(reef_geog,region.id=='3408')

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
fk_93_18<- subset(fk_79_18,YEAR>=1993) #Subset for the dataset from 1993 to match the first year of REEF surveys

####3. Geographic filtering ####
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
reef_geog$hab_class<- rvc_grid$habclass[reef_pts$grid_match]

grid_match_rvc<- st_intersects(rvc_pts,rvc_grid_84,sparse=T)
rvc_sites$grid_match<- NA
for(i in 1:nrow(rvc_sites)){
  if(length(grid_match_rvc[[i]])==1){
    rvc_sites$grid_match[i]=grid_match_rvc[[i]]  
  }else{
    rvc_sites$grid_match[i]=NA
  }
}
rvc_sites$hab_class<- rvc_grid$habclass[rvc_sites$grid_match]

reef_sites<- R %>% group_by(geogr) %>% summarize(n=n_distinct(formid)) #Number of surveys per site
reef_sites_filter<- subset(reef_sites,n>5) #only keep sites surveyed at least 5 times
R<- subset(R, geogr %in% reef_sites_filter$geogr)

#### 4. Creating RVC time-series ####
#Fish data from REEF - remove ultra rare and basket species designations
fish_reef<- read.csv("Caribbean_fish_trait_matrix.csv") #fish species found in the Tropical Western Atlantic
fish_reef<- subset(fish_dat,expert_sighting_freq>1) #Take out the very rare species
fish_rvc<- read.csv("Florida_keys_taxonomic_data.csv")
fish_rvc<- subset(fish_rvc,gsub('.*\\ ', '', fish_rvc$SCINAME)!='sp.') #remove unknown species
fish_rvc<- subset(fish_rvc, SCINAME %in% fish_reef$sciname2)
m<- match(fish_reef$sciname2,fish_rvc$SCINAME)
fish_reef$rvc_code<- fish_rvc$SPECIES_CD[m]

#Separate out RVC data into the REEF sub-regions
rvc_3403_geog<- subset(rvc_sites,region.id=='3403') #Sites in the Key Largo sub-region
write.csv(rvc_3403_geog,'rvc_sites_3403.csv')
rvc_3403_geog_SG<- subset(rvc_3403_geog,hab_class=='SPGR_LR'|hab_class=='SPGR_HR') #Subset down to the spur and groove habitat

rvc_3404_geog<- subset(rvc_sites,region.id=='3404') #Sites in the Islamorada sub-region
write.csv(rvc_3404_geog,'rvc_sites_3404.csv')
rvc_3404_geog_SG<- subset(rvc_3404_geog,hab_class=='SPGR_LR'|hab_class=='SPGR_HR') #Subset down to the spur and groove habitat

rvc_3408_geog<- subset(rvc_sites,region.id=='3408') #Sites in the Key West sub-region
write.csv(rvc_3408_geog,'rvc_sites_3408.csv')
rvc_3408_geog_SG<- subset(rvc_3408_geog,hab_class=='SPGR_LR'|hab_class=='SPGR_HR') #Subset down to the spur and groove habitat

#Determine abundance trends for selected species in spur & groove surveys from each sub-region
rvc_trends_3403<- rvc_batch(fk_93_18,GZ='3403',sp=fish_reef,geog=rvc_3403_geog_SG)
rvc_3403_green<- rlist::list.filter(rvc_trends_3403,length(na.omit(mean_ssu_abundance))>=18) #Green list filtering down to well sampled species - 175 species

rvc_trends_3404<- rvc_batch(fk_93_18,GZ='3404',sp=fish_reef,geog=rvc_3404_geog_SG)
rvc_3404_green<- rlist::list.filter(rvc_trends_3404,length(na.omit(mean_ssu_abundance))>=18) #Green list filtering down to well sampled species - 175 species

rvc_trends_3408<- rvc_batch(fk_93_18,GZ='3408',sp=fish_reef,geog=rvc_3408_geog_SG)
rvc_3408_green<- rlist::list.filter(rvc_trends_3408,length(na.omit(mean_ssu_abundance))>=18) #Green list filtering down to well sampled species - 175 species


####5. Creating REEF time-series
#
#
reef_geog_3403<- reef_3403_pts %>% subset(is.na(grid_match)==F& hab_class=='SPGR_HR'|hab_class=='SPGR_LR'|hab_class=='ISOL_MR') %>% subset(n >= 5)
reef_geog_3404<- subset(reef_3404_pts,is.na(grid_match)==F& n>=5)
reef_geog_3408<- subset(reef_3408_pts,is.na(grid_match)==F& n>=5)

####Reef & RVC trends####
REEF_3403<- REEF_pull(R,GZ='3403',sp=spp[1],geog=reef_geog_3403)
REEF_3403_exp<- REEF_pull(R_exp,GZ='3403',sp=spp,geog=reef_geog_3403)

REEF_3404<- REEF_pull(R,GZ='3404',sp=spp)
REEF_3404_exp<- REEF_pull(R_exp,GZ='3404',sp=spp)
#REEF_3405<- REEF_pull(R,GZ='3405',sp=spp)
#REEF_3406<- REEF_pull(R,GZ='3406',sp=spp)
REEF_3408<- REEF_pull(R,GZ='3408',sp=spp)
REEF_3408_exp<- REEF_pull(R_exp,GZ='3408',sp=spp)

#Filter out low data (At least 20 years of sighting data)
#REEF_3302_trim<- list.filter(REEF_3302,min(tot.surveys)>=10) #Basically all of 3302
REEF_3403_green<- list.filter(REEF_3403,length(na.omit(meanAbund))>=18) #Green list - 175 species
REEF_3403_exp_green<- list.filter(REEF_3403_exp,length(na.omit(meanAbund))>=18) #Green list - 175 species

REEF_3404_green<- list.filter(REEF_3404,length(na.omit(meanAbund))>=18)#151 left
REEF_3404_exp_green<- list.filter(REEF_3404_exp,length(na.omit(meanAbund))>=18) #Green list - 175 species
#REEF_3405_trim<- list.filter(REEF_3405,min(tot.surveys)>=10)
#REEF_3406_trim<- list.filter(REEF_3405,min(tot.surveys)>=10)
#REEF_3406_green<- list.filter(REEF_3406,length(na.omit(meanAbund))>=18)# 129 left
REEF_3408_green<- list.filter(REEF_3408,length(na.omit(meanAbund))>=18)# 108 left
REEF_3408_exp_green<- list.filter(REEF_3408_exp,length(na.omit(meanAbund))>=18) #Green list - 175 species

#Trim down RVC dataset to these species
fk_93_18_trim<- subset(fk_93_18,SPECIES_CD %in% fish_reef$rvc_code)

#Separate out RVC data into the REEF sub-regions
rvc_3403<- subset(fk_93_18_trim,region.id=='3403')
rvc_3404<- subset(fk_93_18_trim,region.id=='3404')
rvc_3408<- subset(fk_93_18_trim,region.id=='3408')

#Only include sites in the sampling grid
rvc_geog_3403<- subset(rvc_3403_pts,hab_class=='SPGR_LR'|hab_class=='SPGR_HR'|hab_class=='ISOL_MR')
rvc_geog_3404<- subset(rvc_3404_pts,hab_class=='SPGR_LR'|hab_class=='SPGR_HR')
rvc_geog_3408<- subset(rvc_3408_pts,hab_class=='SPGR_LR'|hab_class=='SPGR_HR')

#Only include spur and groove
rvc_3403_sg_plus<- filter(rvc_3403, LAT_LON %in% rvc_geog_3403$LAT_LON)
rvc_3404_sg<- filter(rvc_3404, LAT_LON %in% rvc_geog_3404$LAT_LON)
rvc_3408_sg<- filter(rvc_3408, LAT_LON %in% rvc_geog_3408$LAT_LON)

#
rvc_trends_3403<- rvc_batch(rvc_3403_sg_plus)
rvc_trends_3404<- rvc_batch(rvc_3404_sg)
rvc_trends_3408<- rvc_batch(rvc_3408_sg)

rvc_3403_green<-  list.filter(rvc_trends_3403,length(na.omit(mean_ssu_abundance))>=18) #Records from 1993 to 2018 (23 years of data, 5 missing years = 18), 109 spp
rvc_3404_green<-  list.filter(rvc_trends_3404,length(na.omit(mean_ssu_abundance))>=16) #Records from 1995 to 2018 (21 years of data, 5 missing = 16),75 spp
rvc_3408_green<-  list.filter(rvc_trends_3408,length(na.omit(mean_ssu_abundance))>=18) #Records from 1993 to 2018 (23 years of data, 5 missing years = 18), 115 spp 


####Mars batches - site weighted, expert data####
REEF_3403<- REEF_pull(R,GZ='3403',sp=spp,geog=reef_geog_3403)
REEF_3403_exp<- REEF_pull(R_exp,GZ='3403',sp=spp,geog=reef_geog_3403)

REEF_3403_green<- list.filter(REEF_3403,length(na.omit(meanAbund))>=18) #Green list - 175 species
REEF_3403_exp_green<- list.filter(REEF_3403_exp,length(na.omit(meanAbund))>=18) #Green list - 175 species

REEF.3403.comb.green<- do.call(rbind, lapply(REEF_3403_green, data.frame, stringsAsFactors=FALSE))
REEF.green.3403.sp<- unique(REEF.3403.comb.green$sciName)

rvc.3403.comb<- do.call(rbind, lapply(rvc_3403_green, data.frame, stringsAsFactors=FALSE))
rvc.green.3403.sp<- unique(rvc.3403.comb$Species)

ts_comp<- list()
years<- seq(1993,2018,by=1)
n<- NA
R_matrix<- list()

setwd("C:/Users/14388/Desktop/Scripps - Project 1 RVC and REEF/RVC/timeseries/1993 expert timeseries - site-weighted sg_isol")
mars_3403<- data.frame(SP=NA,conv.1=NA,conv.2=NA,conv.3=NA,tier.1.AICc=NA,tier2.AICc=NA,tier3.AICc=NA,Q1=NA,Q2=NA,Q3.rvc=NA,Q3.reef=NA,U1=NA,U2.rvc=NA,U2.reef=NA,U3.rvc=NA,U3.reef=NA,R1.rvc=NA,R1.reef=NA,R2.rvc=NA,R2.reef=NA,R3.rvc=NA,R3.reef=NA,cName=NA,obs.grp=NA,mAbund.rvc=NA,mAbund.reef=NA,sdAbund.rvc=NA,sdAbund.reef=NA,dAIC1=NA,dAIC2=NA,dAIC3=NA,mod=NA,years.rvc=NA,years.reef=NA)
for(i in 1:length(rvc.green.3403.sp)){
  spp<- filter(fish_reef,rvc_code==rvc.green.3403.sp[i])
  ts_comp[[i]]<- t(log10(rvc_3403_green[[i]]$mean_ssu_abundance)) #extract logged time-series
  m<- match(spp$scientificname,REEF.green.3403.sp)
  if(is.na(m)==T){next}
  ts_comp[[i]]<- rbind(ts_comp[[i]],t(log10(REEF_3403_green[[m]]$SiteAbund[1:26])))
  rownames(ts_comp[[i]])<- c(paste(spp$scientificname,'rvc',sep=":"),paste(spp$scientificname,'REEF',sep=":"))
  colnames(ts_comp[[i]])<- years
  
  Z_1=factor(c('rvc_reef','rvc_reef')) #Common state process for both time-series
  Z_2=factor(c('rvc','reef')) #Different state process for each time-series
  R_matrix[[i]]=matrix(list(0),nrow(ts_comp[[i]]),nrow(ts_comp[[i]])) #create empty Obs. error matrix
  diag(R_matrix[[i]])=c(rep("rvc",times=1),rep("reef",times=1))
  
  #tier 1 - common state process
  tier_1<- list(Z=Z_1,
                Q = 'diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 2 - Diff. U, same process error
  tier_2<- list(Z=Z_2,
                Q ='diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 3
  tier_3<- list(Z=Z_2,
                Q = 'diagonal and unequal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  fit_1<-  MARSS(ts_comp[[i]], model = tier_1,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  fit_2<-  MARSS(ts_comp[[i]], model = tier_2,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  fit_3<-  MARSS(ts_comp[[i]], model = tier_3,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  params.1<- MARSSparamCIs(fit_1)
  params.2<- MARSSparamCIs(fit_2)
  params.3<- MARSSparamCIs(fit_3)
  
  mars_3403[i,1]=spp$scientificname
  mars_3403[i,2]=ifelse(is.null(fit_1$errors)==T,1,0)
  mars_3403[i,3]=ifelse(is.null(fit_2$errors)==T,1,0)  
  mars_3403[i,4]=ifelse(is.null(fit_3$errors)==T,1,0)  
  mars_3403[i,5]=fit_1$AICc
  mars_3403[i,6]=fit_2$AICc  
  mars_3403[i,7]=fit_3$AICc  
  mars_3403[i,8]=params.1$parMean[5]
  mars_3403[i,9]=params.2$parMean[5]
  mars_3403[i,10]=params.3$parMean[5]
  mars_3403[i,11]=params.3$parMean[6]
  mars_3403[i,12]=params.1$parMean[4]
  mars_3403[i,13]=params.2$parMean[3] 
  mars_3403[i,14]=params.2$parMean[4]
  mars_3403[i,15]=params.3$parMean[3]
  mars_3403[i,16]=params.3$parMean[4]
  mars_3403[i,17]=params.1$parMean[2] 
  mars_3403[i,18]=params.1$parMean[3]
  mars_3403[i,19]=params.2$parMean[1] 
  mars_3403[i,20]=params.2$parMean[2]
  mars_3403[i,21]=params.3$parMean[1] 
  mars_3403[i,22]=params.3$parMean[2]
  mars_3403[i,23]=spp$commonname
  mars_3403[i,24]=spp$group_assignment
  mars_3403[i,25]=exp(mean(na.omit(ts_comp[[i]][1,])))
  mars_3403[i,26]=exp(mean(na.omit(ts_comp[[i]][2,])))
  mars_3403[i,27]=sd(exp(na.omit(ts_comp[[i]][1,])))
  mars_3403[i,28]=sd(exp(na.omit(ts_comp[[i]][2,])))
  mars_3403[i,29]=mars_3403[i,5]-min(mars_3403[i,5:7])
  mars_3403[i,30]=mars_3403[i,6]-min(mars_3403[i,5:7])
  mars_3403[i,31]=mars_3403[i,7]-min(mars_3403[i,5:7])
  if(mars_3403[i,29]==0){mars_3403[i,32]=1}
  if(mars_3403[i,30]==0){mars_3403[i,32]=2}
  if(mars_3403[i,31]==0){mars_3403[i,32]=3}
  mars_3403[i,33]=length(na.omit(ts_comp[[i]][1,]))
  mars_3403[i,34]=length(na.omit(ts_comp[[i]][2,]))
  
  
  TS.plot.pdf(x=ts_comp[[i]],m=mars_3403$mod[i],sp=mars_3403$cName[i],GZ='Key Largo')
  dev.off()
  print(i)
}

setwd("C:/Users/14388/Desktop/Scripps - Project 1 RVC and REEF/RVC")

write.csv(mars_3403,'MARSS_model_comparisons_3403_1993_expert_sg_isol_siteweight.csv')

#expert data
setwd("C:/Users/14388/Desktop/Scripps - Project 1 RVC and REEF/RVC/timeseries/1993 expert timeseries - site-weighted spur and groove")
mars_3403<- data.frame(SP=NA,conv.1=NA,conv.2=NA,conv.3=NA,tier.1.AICc=NA,tier2.AICc=NA,tier3.AICc=NA,Q1=NA,Q2=NA,Q3.rvc=NA,Q3.reef=NA,U1=NA,U2.rvc=NA,U2.reef=NA,U3.rvc=NA,U3.reef=NA,R1.rvc=NA,R1.reef=NA,R2.rvc=NA,R2.reef=NA,R3.rvc=NA,R3.reef=NA,cName=NA,obs.grp=NA,mAbund.rvc=NA,mAbund.reef=NA,sdAbund.rvc=NA,sdAbund.reef=NA,dAIC1=NA,dAIC2=NA,dAIC3=NA,mod=NA,years.rvc=NA,years.reef=NA)
for(i in 1:length(rvc.green.3403.sp)){
  spp<- filter(fish_reef,rvc_code==rvc.green.3403.sp[i])
  ts_comp[[i]]<- t(log10(rvc_3403_green[[i]]$mean_ssu_abundance)) #extract logged time-series
  m<- match(spp$scientificname,REEF.green.exp.3403.sp)
  if(is.na(m)==T){next}
  ts_comp[[i]]<- rbind(ts_comp[[i]],t(log10(REEF_3403_exp_green[[m]]$SiteAbund[1:26])))
  rownames(ts_comp[[i]])<- c(paste(spp$scientificname,'rvc',sep=":"),paste(spp$scientificname,'REEF',sep=":"))
  colnames(ts_comp[[i]])<- years
  
  Z_1=factor(c('rvc_reef','rvc_reef')) #Common state process for both time-series
  Z_2=factor(c('rvc','reef')) #Different state process for each time-series
  R_matrix[[i]]=matrix(list(0),nrow(ts_comp[[i]]),nrow(ts_comp[[i]])) #create empty Obs. error matrix
  diag(R_matrix[[i]])=c(rep("rvc",times=1),rep("reef",times=1))
  
  #tier 1 - common state process
  tier_1<- list(Z=Z_1,
                Q = 'diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 2 - Diff. U, same process error
  tier_2<- list(Z=Z_2,
                Q ='diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 3
  tier_3<- list(Z=Z_2,
                Q = 'diagonal and unequal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  fit_1<-  MARSS(ts_comp[[i]], model = tier_1,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  fit_2<-  MARSS(ts_comp[[i]], model = tier_2,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  fit_3<-  MARSS(ts_comp[[i]], model = tier_3,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  params.1<- MARSSparamCIs(fit_1)
  params.2<- MARSSparamCIs(fit_2)
  params.3<- MARSSparamCIs(fit_3)
  
  mars_3403[i,1]=spp$scientificname
  mars_3403[i,2]=ifelse(is.null(fit_1$errors)==T,1,0)
  mars_3403[i,3]=ifelse(is.null(fit_2$errors)==T,1,0)  
  mars_3403[i,4]=ifelse(is.null(fit_3$errors)==T,1,0)  
  mars_3403[i,5]=fit_1$AICc
  mars_3403[i,6]=fit_2$AICc  
  mars_3403[i,7]=fit_3$AICc  
  mars_3403[i,8]=params.1$parMean[5]
  mars_3403[i,9]=params.2$parMean[5]
  mars_3403[i,10]=params.3$parMean[5]
  mars_3403[i,11]=params.3$parMean[6]
  mars_3403[i,12]=params.1$parMean[4]
  mars_3403[i,13]=params.2$parMean[3] 
  mars_3403[i,14]=params.2$parMean[4]
  mars_3403[i,15]=params.3$parMean[3]
  mars_3403[i,16]=params.3$parMean[4]
  mars_3403[i,17]=params.1$parMean[2] 
  mars_3403[i,18]=params.1$parMean[3]
  mars_3403[i,19]=params.2$parMean[1] 
  mars_3403[i,20]=params.2$parMean[2]
  mars_3403[i,21]=params.3$parMean[1] 
  mars_3403[i,22]=params.3$parMean[2]
  mars_3403[i,23]=spp$commonname
  mars_3403[i,24]=spp$group_assignment
  mars_3403[i,25]=exp(mean(na.omit(ts_comp[[i]][1,])))
  mars_3403[i,26]=exp(mean(na.omit(ts_comp[[i]][2,])))
  mars_3403[i,27]=sd(exp(na.omit(ts_comp[[i]][1,])))
  mars_3403[i,28]=sd(exp(na.omit(ts_comp[[i]][2,])))
  mars_3403[i,29]=mars_3403[i,5]-min(mars_3403[i,5:7])
  mars_3403[i,30]=mars_3403[i,6]-min(mars_3403[i,5:7])
  mars_3403[i,31]=mars_3403[i,7]-min(mars_3403[i,5:7])
  if(mars_3403[i,29]==0){mars_3403[i,32]=1}
  if(mars_3403[i,30]==0){mars_3403[i,32]=2}
  if(mars_3403[i,31]==0){mars_3403[i,32]=3}
  mars_3403[i,33]=length(na.omit(ts_comp[[i]][1,]))
  mars_3403[i,34]=length(na.omit(ts_comp[[i]][2,]))
  
  
  TS.plot.pdf(x=ts_comp[[i]],m=mars_3403$mod[i],sp=mars_3403$cName[i],GZ='Key Largo')
  dev.off()
  print(i)
}

setwd("C:/Users/14388/Desktop/Scripps - Project 1 RVC and REEF/RVC")

write.csv(mars_3403,'MARSS_model_comparisons_3403_1993_expert_spurgroove_siteweight.csv')


#Spur & Groove; Contiguous low relief; Isolated medium relief
REEF_3403_exp<- REEF_pull(R_exp,GZ='3403',sp=spp,geog=reef_geog_3403)
REEF_3403_exp_green<- list.filter(REEF_3403_exp,length(na.omit(meanAbund))>=18) #Green list - 175 species

REEF.3403.exp.comb.green<- do.call(rbind, lapply(REEF_3403_exp_green, data.frame, stringsAsFactors=FALSE))
REEF.green.exp.3403.sp<- unique(REEF.3403.exp.comb.green$sciName)

rvc.3403.comb<- do.call(rbind, lapply(rvc_3403_green, data.frame, stringsAsFactors=FALSE))
rvc.green.3403.sp<- unique(rvc.3403.comb$Species)

ts_comp<- list()
years<- seq(1993,2018,by=1)
n<- NA
R_matrix<- list()

setwd("C:/Users/14388/Desktop/Scripps - Project 1 RVC and REEF/RVC/timeseries/1993 expert timeseries - site-weighted sg_cont_isol")
mars_3403<- data.frame(SP=NA,conv.1=NA,conv.2=NA,conv.3=NA,tier.1.AICc=NA,tier2.AICc=NA,tier3.AICc=NA,Q1=NA,Q2=NA,Q3.rvc=NA,Q3.reef=NA,U1=NA,U2.rvc=NA,U2.reef=NA,U3.rvc=NA,U3.reef=NA,R1.rvc=NA,R1.reef=NA,R2.rvc=NA,R2.reef=NA,R3.rvc=NA,R3.reef=NA,cName=NA,obs.grp=NA,mAbund.rvc=NA,mAbund.reef=NA,sdAbund.rvc=NA,sdAbund.reef=NA,dAIC1=NA,dAIC2=NA,dAIC3=NA,mod=NA,years.rvc=NA,years.reef=NA)
for(i in 1:length(rvc.green.3403.sp)){
  spp<- filter(fish_reef,rvc_code==rvc.green.3403.sp[i])
  ts_comp[[i]]<- t(log10(rvc_3403_green[[i]]$mean_ssu_abundance)) #extract logged time-series
  m<- match(spp$scientificname,REEF.green.exp.3403.sp)
  if(is.na(m)==T){next}
  ts_comp[[i]]<- rbind(ts_comp[[i]],t(log10(REEF_3403_exp_green[[m]]$SiteAbund[1:26])))
  rownames(ts_comp[[i]])<- c(paste(spp$scientificname,'rvc',sep=":"),paste(spp$scientificname,'REEF',sep=":"))
  colnames(ts_comp[[i]])<- years
  
  Z_1=factor(c('rvc_reef','rvc_reef')) #Common state process for both time-series
  Z_2=factor(c('rvc','reef')) #Different state process for each time-series
  R_matrix[[i]]=matrix(list(0),nrow(ts_comp[[i]]),nrow(ts_comp[[i]])) #create empty Obs. error matrix
  diag(R_matrix[[i]])=c(rep("rvc",times=1),rep("reef",times=1))
  
  #tier 1 - common state process
  tier_1<- list(Z=Z_1,
                Q = 'diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 2 - Diff. U, same process error
  tier_2<- list(Z=Z_2,
                Q ='diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 3
  tier_3<- list(Z=Z_2,
                Q = 'diagonal and unequal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  fit_1<-  MARSS(ts_comp[[i]], model = tier_1,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  fit_2<-  MARSS(ts_comp[[i]], model = tier_2,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  fit_3<-  MARSS(ts_comp[[i]], model = tier_3,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  params.1<- MARSSparamCIs(fit_1)
  params.2<- MARSSparamCIs(fit_2)
  params.3<- MARSSparamCIs(fit_3)
  
  mars_3403[i,1]=spp$scientificname
  mars_3403[i,2]=ifelse(is.null(fit_1$errors)==T,1,0)
  mars_3403[i,3]=ifelse(is.null(fit_2$errors)==T,1,0)  
  mars_3403[i,4]=ifelse(is.null(fit_3$errors)==T,1,0)  
  mars_3403[i,5]=fit_1$AICc
  mars_3403[i,6]=fit_2$AICc  
  mars_3403[i,7]=fit_3$AICc  
  mars_3403[i,8]=params.1$parMean[5]
  mars_3403[i,9]=params.2$parMean[5]
  mars_3403[i,10]=params.3$parMean[5]
  mars_3403[i,11]=params.3$parMean[6]
  mars_3403[i,12]=params.1$parMean[4]
  mars_3403[i,13]=params.2$parMean[3] 
  mars_3403[i,14]=params.2$parMean[4]
  mars_3403[i,15]=params.3$parMean[3]
  mars_3403[i,16]=params.3$parMean[4]
  mars_3403[i,17]=params.1$parMean[2] 
  mars_3403[i,18]=params.1$parMean[3]
  mars_3403[i,19]=params.2$parMean[1] 
  mars_3403[i,20]=params.2$parMean[2]
  mars_3403[i,21]=params.3$parMean[1] 
  mars_3403[i,22]=params.3$parMean[2]
  mars_3403[i,23]=spp$commonname
  mars_3403[i,24]=spp$group_assignment
  mars_3403[i,25]=exp(mean(na.omit(ts_comp[[i]][1,])))
  mars_3403[i,26]=exp(mean(na.omit(ts_comp[[i]][2,])))
  mars_3403[i,27]=sd(exp(na.omit(ts_comp[[i]][1,])))
  mars_3403[i,28]=sd(exp(na.omit(ts_comp[[i]][2,])))
  mars_3403[i,29]=mars_3403[i,5]-min(mars_3403[i,5:7])
  mars_3403[i,30]=mars_3403[i,6]-min(mars_3403[i,5:7])
  mars_3403[i,31]=mars_3403[i,7]-min(mars_3403[i,5:7])
  if(mars_3403[i,29]==0){mars_3403[i,32]=1}
  if(mars_3403[i,30]==0){mars_3403[i,32]=2}
  if(mars_3403[i,31]==0){mars_3403[i,32]=3}
  mars_3403[i,33]=length(na.omit(ts_comp[[i]][1,]))
  mars_3403[i,34]=length(na.omit(ts_comp[[i]][2,]))
  
  
  TS.plot.pdf(x=ts_comp[[i]],m=mars_3403$mod[i],sp=mars_3403$cName[i],GZ='Key Largo')
  dev.off()
  print(i)
}

setwd("C:/Users/14388/Desktop/Scripps - Project 1 RVC and REEF/RVC")

write.csv(mars_3403,'MARSS_model_comparisons_3403_1993_expert_sg_isol_cont_siteweight.csv')


####Temporal Matching####

#Expert
reef_matched_3403<- list()

for(i in 1:length(rvc_3403_green)){
  reef_matched_3403[[i]]<- REEF_pull_time_window(rvc_3403_green[[i]],R=R_exp,GZ='3403',geog=reef_geog_3403)
}


REEF_3403_green<- list.filter(reef_matched_3403,length(na.omit(meanAbund))>=18) #Green list - 77 species

REEF.3403.comb.green<- do.call(rbind, lapply(REEF_3403_green, data.frame, stringsAsFactors=FALSE))
REEF.green.3403.sp<- unique(REEF.3403.comb.green$sciName)

rvc.3403.comb<- do.call(rbind, lapply(rvc_3403_green, data.frame, stringsAsFactors=FALSE))
rvc.green.3403.sp<- unique(rvc.3403.comb$Species)

ts_comp<- list()
years<- seq(1993,2018,by=1)
n<- NA
R_matrix<- list()


setwd("C:/Users/14388/Desktop/Scripps - Project 1 RVC and REEF/RVC/timeseries/1993 expert timeseries - temporal match spur_groove 3403")
mars_3403<- data.frame(SP=NA,conv.1=NA,conv.2=NA,conv.3=NA,tier.1.AICc=NA,tier2.AICc=NA,tier3.AICc=NA,Q1=NA,Q2=NA,Q3.rvc=NA,Q3.reef=NA,U1=NA,U2.rvc=NA,U2.reef=NA,U3.rvc=NA,U3.reef=NA,R1.rvc=NA,R1.reef=NA,R2.rvc=NA,R2.reef=NA,R3.rvc=NA,R3.reef=NA,cName=NA,obs.grp=NA,mAbund.rvc=NA,mAbund.reef=NA,sdAbund.rvc=NA,sdAbund.reef=NA,dAIC1=NA,dAIC2=NA,dAIC3=NA,mod=NA,years.rvc=NA,years.reef=NA)
for(i in 1:length(rvc.green.3403.sp)){
  spp<- filter(fish_reef,rvc_code==rvc.green.3403.sp[i])
  ts_comp[[i]]<- t(log10(rvc_3403_green[[i]]$mean_ssu_abundance)) #extract logged time-series
  m<- match(spp$scientificname,REEF.green.3403.sp)
  if(is.na(m)==T){next}
  ts_comp[[i]]<- rbind(ts_comp[[i]],t(log10(REEF_3403_green[[m]]$SiteAbund[1:26])))
  rownames(ts_comp[[i]])<- c(paste(spp$scientificname,'rvc',sep=":"),paste(spp$scientificname,'REEF',sep=":"))
  colnames(ts_comp[[i]])<- years
  
  Z_1=factor(c('rvc_reef','rvc_reef')) #Common state process for both time-series
  Z_2=factor(c('rvc','reef')) #Different state process for each time-series
  R_matrix[[i]]=matrix(list(0),nrow(ts_comp[[i]]),nrow(ts_comp[[i]])) #create empty Obs. error matrix
  diag(R_matrix[[i]])=c(rep("rvc",times=1),rep("reef",times=1))
  
  #tier 1 - common state process
  tier_1<- list(Z=Z_1,
                Q = 'diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 2 - Diff. U, same process error
  tier_2<- list(Z=Z_2,
                Q ='diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 3
  tier_3<- list(Z=Z_2,
                Q = 'diagonal and unequal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  fit_1<-  MARSS(ts_comp[[i]], model = tier_1,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  fit_2<-  MARSS(ts_comp[[i]], model = tier_2,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  fit_3<-  MARSS(ts_comp[[i]], model = tier_3,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  params.1<- MARSSparamCIs(fit_1)
  params.2<- MARSSparamCIs(fit_2)
  params.3<- MARSSparamCIs(fit_3)
  
  mars_3403[i,1]=spp$scientificname
  mars_3403[i,2]=ifelse(is.null(fit_1$errors)==T,1,0)
  mars_3403[i,3]=ifelse(is.null(fit_2$errors)==T,1,0)  
  mars_3403[i,4]=ifelse(is.null(fit_3$errors)==T,1,0)  
  mars_3403[i,5]=fit_1$AICc
  mars_3403[i,6]=fit_2$AICc  
  mars_3403[i,7]=fit_3$AICc  
  mars_3403[i,8]=params.1$parMean[5]
  mars_3403[i,9]=params.2$parMean[5]
  mars_3403[i,10]=params.3$parMean[5]
  mars_3403[i,11]=params.3$parMean[6]
  mars_3403[i,12]=params.1$parMean[4]
  mars_3403[i,13]=params.2$parMean[3] 
  mars_3403[i,14]=params.2$parMean[4]
  mars_3403[i,15]=params.3$parMean[3]
  mars_3403[i,16]=params.3$parMean[4]
  mars_3403[i,17]=params.1$parMean[2] 
  mars_3403[i,18]=params.1$parMean[3]
  mars_3403[i,19]=params.2$parMean[1] 
  mars_3403[i,20]=params.2$parMean[2]
  mars_3403[i,21]=params.3$parMean[1] 
  mars_3403[i,22]=params.3$parMean[2]
  mars_3403[i,23]=spp$commonname
  mars_3403[i,24]=spp$group_assignment
  mars_3403[i,25]=exp(mean(na.omit(ts_comp[[i]][1,])))
  mars_3403[i,26]=exp(mean(na.omit(ts_comp[[i]][2,])))
  mars_3403[i,27]=sd(exp(na.omit(ts_comp[[i]][1,])))
  mars_3403[i,28]=sd(exp(na.omit(ts_comp[[i]][2,])))
  mars_3403[i,29]=mars_3403[i,5]-min(mars_3403[i,5:7])
  mars_3403[i,30]=mars_3403[i,6]-min(mars_3403[i,5:7])
  mars_3403[i,31]=mars_3403[i,7]-min(mars_3403[i,5:7])
  if(mars_3403[i,29]==0){mars_3403[i,32]=1}
  if(mars_3403[i,30]==0){mars_3403[i,32]=2}
  if(mars_3403[i,31]==0){mars_3403[i,32]=3}
  mars_3403[i,33]=length(na.omit(ts_comp[[i]][1,]))
  mars_3403[i,34]=length(na.omit(ts_comp[[i]][2,]))
  
  
  TS.plot.pdf(x=ts_comp[[i]],m=mars_3403$mod[i],sp=mars_3403$cName[i],GZ='Key Largo')
  dev.off()
  print(i)
}

setwd("C:/Users/14388/Desktop/Scripps - Project 1 RVC and REEF/RVC")

write.csv(mars_3403,'MARSS_model_comparisons_3403_time_matched_expert_spur_groove.csv')

#Nov + Expert
reef_matched_3403<- list()

for(i in 1:length(rvc_3403_green)){
  reef_matched_3403[[i]]<- REEF_pull_time_window(rvc_3403_green[[i]],R=R,GZ='3403',geog=reef_geog_3403)
}


REEF_3403_green<- list.filter(reef_matched_3403,length(na.omit(meanAbund))>=18) #Green list - 77 species

REEF.3403.comb.green<- do.call(rbind, lapply(REEF_3403_green, data.frame, stringsAsFactors=FALSE))
REEF.green.3403.sp<- unique(REEF.3403.comb.green$sciName)

rvc.3403.comb<- do.call(rbind, lapply(rvc_3403_green, data.frame, stringsAsFactors=FALSE))
rvc.green.3403.sp<- unique(rvc.3403.comb$Species)

ts_comp<- list()
years<- seq(1993,2018,by=1)
n<- NA
R_matrix<- list()


setwd("C:/Users/14388/Desktop/Scripps - Project 1 RVC and REEF/RVC/timeseries/1993 nov+expert timeseries - temporal match spur_groove 3403")
mars_3403<- data.frame(SP=NA,conv.1=NA,conv.2=NA,conv.3=NA,tier.1.AICc=NA,tier2.AICc=NA,tier3.AICc=NA,Q1=NA,Q2=NA,Q3.rvc=NA,Q3.reef=NA,U1=NA,U2.rvc=NA,U2.reef=NA,U3.rvc=NA,U3.reef=NA,R1.rvc=NA,R1.reef=NA,R2.rvc=NA,R2.reef=NA,R3.rvc=NA,R3.reef=NA,cName=NA,obs.grp=NA,mAbund.rvc=NA,mAbund.reef=NA,sdAbund.rvc=NA,sdAbund.reef=NA,dAIC1=NA,dAIC2=NA,dAIC3=NA,mod=NA,years.rvc=NA,years.reef=NA)
for(i in 1:length(rvc.green.3403.sp)){
  spp<- filter(fish_reef,rvc_code==rvc.green.3403.sp[i])
  ts_comp[[i]]<- t(log10(rvc_3403_green[[i]]$mean_ssu_abundance)) #extract logged time-series
  m<- match(spp$scientificname,REEF.green.3403.sp)
  if(is.na(m)==T){next}
  ts_comp[[i]]<- rbind(ts_comp[[i]],t(log10(REEF_3403_green[[m]]$SiteAbund[1:26])))
  rownames(ts_comp[[i]])<- c(paste(spp$scientificname,'rvc',sep=":"),paste(spp$scientificname,'REEF',sep=":"))
  colnames(ts_comp[[i]])<- years
  
  Z_1=factor(c('rvc_reef','rvc_reef')) #Common state process for both time-series
  Z_2=factor(c('rvc','reef')) #Different state process for each time-series
  R_matrix[[i]]=matrix(list(0),nrow(ts_comp[[i]]),nrow(ts_comp[[i]])) #create empty Obs. error matrix
  diag(R_matrix[[i]])=c(rep("rvc",times=1),rep("reef",times=1))
  
  #tier 1 - common state process
  tier_1<- list(Z=Z_1,
                Q = 'diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 2 - Diff. U, same process error
  tier_2<- list(Z=Z_2,
                Q ='diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 3
  tier_3<- list(Z=Z_2,
                Q = 'diagonal and unequal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  fit_1<-  MARSS(ts_comp[[i]], model = tier_1,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  fit_2<-  MARSS(ts_comp[[i]], model = tier_2,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  fit_3<-  MARSS(ts_comp[[i]], model = tier_3,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  params.1<- MARSSparamCIs(fit_1)
  params.2<- MARSSparamCIs(fit_2)
  params.3<- MARSSparamCIs(fit_3)
  
  mars_3403[i,1]=spp$scientificname
  mars_3403[i,2]=ifelse(is.null(fit_1$errors)==T,1,0)
  mars_3403[i,3]=ifelse(is.null(fit_2$errors)==T,1,0)  
  mars_3403[i,4]=ifelse(is.null(fit_3$errors)==T,1,0)  
  mars_3403[i,5]=fit_1$AICc
  mars_3403[i,6]=fit_2$AICc  
  mars_3403[i,7]=fit_3$AICc  
  mars_3403[i,8]=params.1$parMean[5]
  mars_3403[i,9]=params.2$parMean[5]
  mars_3403[i,10]=params.3$parMean[5]
  mars_3403[i,11]=params.3$parMean[6]
  mars_3403[i,12]=params.1$parMean[4]
  mars_3403[i,13]=params.2$parMean[3] 
  mars_3403[i,14]=params.2$parMean[4]
  mars_3403[i,15]=params.3$parMean[3]
  mars_3403[i,16]=params.3$parMean[4]
  mars_3403[i,17]=params.1$parMean[2] 
  mars_3403[i,18]=params.1$parMean[3]
  mars_3403[i,19]=params.2$parMean[1] 
  mars_3403[i,20]=params.2$parMean[2]
  mars_3403[i,21]=params.3$parMean[1] 
  mars_3403[i,22]=params.3$parMean[2]
  mars_3403[i,23]=spp$commonname
  mars_3403[i,24]=spp$group_assignment
  mars_3403[i,25]=exp(mean(na.omit(ts_comp[[i]][1,])))
  mars_3403[i,26]=exp(mean(na.omit(ts_comp[[i]][2,])))
  mars_3403[i,27]=sd(exp(na.omit(ts_comp[[i]][1,])))
  mars_3403[i,28]=sd(exp(na.omit(ts_comp[[i]][2,])))
  mars_3403[i,29]=mars_3403[i,5]-min(mars_3403[i,5:7])
  mars_3403[i,30]=mars_3403[i,6]-min(mars_3403[i,5:7])
  mars_3403[i,31]=mars_3403[i,7]-min(mars_3403[i,5:7])
  if(mars_3403[i,29]==0){mars_3403[i,32]=1}
  if(mars_3403[i,30]==0){mars_3403[i,32]=2}
  if(mars_3403[i,31]==0){mars_3403[i,32]=3}
  mars_3403[i,33]=length(na.omit(ts_comp[[i]][1,]))
  mars_3403[i,34]=length(na.omit(ts_comp[[i]][2,]))
  
  
  TS.plot.pdf(x=ts_comp[[i]],m=mars_3403$mod[i],sp=mars_3403$cName[i],GZ='Key Largo')
  dev.off()
  print(i)
}

setwd("C:/Users/14388/Desktop/Scripps - Project 1 RVC and REEF/RVC")

write.csv(mars_3403,'MARSS_model_comparisons_3403_time_matched_nov_expert_spur_groove.csv')

#
##Data from 1979 to 2018

#Trim down RVC dataset to these species
fk_79_18_trim<- subset(fk_79_18,SPECIES_CD %in% fish_reef$rvc_code)

#Separate out RVC data into the REEF sub-regions
rvc_3403<- subset(fk_79_18_trim,region.id=='3403')
rvc_3404<- subset(fk_79_18_trim,region.id=='3404')
rvc_3408<- subset(fk_79_18_trim,region.id=='3408')

#Only include sites in the sampling grid
rvc_geog_3403<- subset(rvc_3403_pts,hab_class=='SPGR_LR'|hab_class=='SPGR_HR')
rvc_geog_3404<- subset(rvc_3404_pts,hab_class=='SPGR_LR'|hab_class=='SPGR_HR')
rvc_geog_3408<- subset(rvc_3408_pts,hab_class=='SPGR_LR'|hab_class=='SPGR_HR')

#Only include spur and groove
rvc_3403_sg<- filter(rvc_3403, LAT_LON %in% rvc_geog_3403$LAT_LON)
rvc_3404_sg<- filter(rvc_3404, LAT_LON %in% rvc_geog_3404$LAT_LON)
rvc_3408_sg<- filter(rvc_3408, LAT_LON %in% rvc_geog_3408$LAT_LON)

#
rvc_trends_3403<- rvc_batch(rvc_3403_sg)
rvc_trends_3404<- rvc_batch(rvc_3404_sg)
rvc_trends_3408<- rvc_batch(rvc_3408_sg)

rvc_3403_green<-  list.filter(rvc_trends_3403,length(na.omit(mean_ssu_abundance))>=18) #Records from 1993 to 2018 (23 years of data, 5 missing years = 18), 109 spp
rvc_3404_green<-  list.filter(rvc_trends_3404,length(na.omit(mean_ssu_abundance))>=16) #Records from 1995 to 2018 (21 years of data, 5 missing = 16),75 spp
rvc_3408_green<-  list.filter(rvc_trends_3408,length(na.omit(mean_ssu_abundance))>=18) #Records from 1993 to 2018 (23 years of data, 5 missing years = 18), 115 spp 

reef_matched_3403<- list()

for(i in 1:length(rvc_3403_green)){
  reef_matched_3403[[i]]<- REEF_pull_time_window(rvc_3403_green[[i]],R=R_exp,GZ='3403',geog=reef_geog_3403)
}


REEF_3403_green<- list.filter(reef_matched_3403,length(na.omit(meanAbund))>=18) #Green list - 77 species

REEF.3403.comb.green<- do.call(rbind, lapply(REEF_3403_green, data.frame, stringsAsFactors=FALSE))
REEF.green.3403.sp<- unique(REEF.3403.comb.green$sciName)

rvc.3403.comb<- do.call(rbind, lapply(rvc_3403_green, data.frame, stringsAsFactors=FALSE))
rvc.green.3403.sp<- unique(rvc.3403.comb$Species)

ts_comp<- list()
years<- seq(1980,2018,by=1)
n<- NA
R_matrix<- list()

setwd("C:/Users/14388/Desktop/Scripps - Project 1 RVC and REEF/RVC/timeseries/1980 expert timeseries - temporal match spur_groove 3403")
mars_3403<- data.frame(SP=NA,conv.1=NA,conv.2=NA,conv.3=NA,tier.1.AICc=NA,tier2.AICc=NA,tier3.AICc=NA,Q1=NA,Q2=NA,Q3.rvc=NA,Q3.reef=NA,U1=NA,U2.rvc=NA,U2.reef=NA,U3.rvc=NA,U3.reef=NA,R1.rvc=NA,R1.reef=NA,R2.rvc=NA,R2.reef=NA,R3.rvc=NA,R3.reef=NA,cName=NA,obs.grp=NA,mAbund.rvc=NA,mAbund.reef=NA,sdAbund.rvc=NA,sdAbund.reef=NA,dAIC1=NA,dAIC2=NA,dAIC3=NA,mod=NA,years.rvc=NA,years.reef=NA)
for(i in 1:length(rvc.green.3403.sp)){
  spp<- filter(fish_reef,rvc_code==rvc.green.3403.sp[i])
  ts_comp[[i]]<- t(log10(rvc_3403_green[[i]]$mean_ssu_abundance[2:40])) #extract logged time-series
  m<- match(spp$scientificname,REEF.green.3403.sp)
  if(is.na(m)==T){next}
  ts_comp[[i]]<- rbind(ts_comp[[i]],t(c(rep(NA,13),log10(REEF_3403_green[[m]]$SiteAbund[1:26]))))
  rownames(ts_comp[[i]])<- c(paste(spp$scientificname,'rvc',sep=":"),paste(spp$scientificname,'REEF',sep=":"))
  colnames(ts_comp[[i]])<- years
  
  Z_1=factor(c('rvc_reef','rvc_reef')) #Common state process for both time-series
  Z_2=factor(c('rvc','reef')) #Different state process for each time-series
  R_matrix[[i]]=matrix(list(0),nrow(ts_comp[[i]]),nrow(ts_comp[[i]])) #create empty Obs. error matrix
  diag(R_matrix[[i]])=c(rep("rvc",times=1),rep("reef",times=1))
  
  #tier 1 - common state process
  tier_1<- list(Z=Z_1,
                Q = 'diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 2 - Diff. U, same process error
  tier_2<- list(Z=Z_2,
                Q ='diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 3
  tier_3<- list(Z=Z_2,
                Q = 'diagonal and unequal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  fit_1<-  MARSS(ts_comp[[i]], model = tier_1,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  fit_2<-  MARSS(ts_comp[[i]], model = tier_2,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  fit_3<-  MARSS(ts_comp[[i]], model = tier_3,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  params.1<- MARSSparamCIs(fit_1)
  params.2<- MARSSparamCIs(fit_2)
  params.3<- MARSSparamCIs(fit_3)
  
  mars_3403[i,1]=spp$scientificname
  mars_3403[i,2]=ifelse(is.null(fit_1$errors)==T,1,0)
  mars_3403[i,3]=ifelse(is.null(fit_2$errors)==T,1,0)  
  mars_3403[i,4]=ifelse(is.null(fit_3$errors)==T,1,0)  
  mars_3403[i,5]=fit_1$AICc
  mars_3403[i,6]=fit_2$AICc  
  mars_3403[i,7]=fit_3$AICc  
  mars_3403[i,8]=params.1$parMean[5]
  mars_3403[i,9]=params.2$parMean[5]
  mars_3403[i,10]=params.3$parMean[5]
  mars_3403[i,11]=params.3$parMean[6]
  mars_3403[i,12]=params.1$parMean[4]
  mars_3403[i,13]=params.2$parMean[3] 
  mars_3403[i,14]=params.2$parMean[4]
  mars_3403[i,15]=params.3$parMean[3]
  mars_3403[i,16]=params.3$parMean[4]
  mars_3403[i,17]=params.1$parMean[2] 
  mars_3403[i,18]=params.1$parMean[3]
  mars_3403[i,19]=params.2$parMean[1] 
  mars_3403[i,20]=params.2$parMean[2]
  mars_3403[i,21]=params.3$parMean[1] 
  mars_3403[i,22]=params.3$parMean[2]
  mars_3403[i,23]=spp$commonname
  mars_3403[i,24]=spp$group_assignment
  mars_3403[i,25]=exp(mean(na.omit(ts_comp[[i]][1,])))
  mars_3403[i,26]=exp(mean(na.omit(ts_comp[[i]][2,])))
  mars_3403[i,27]=sd(exp(na.omit(ts_comp[[i]][1,])))
  mars_3403[i,28]=sd(exp(na.omit(ts_comp[[i]][2,])))
  mars_3403[i,29]=mars_3403[i,5]-min(mars_3403[i,5:7])
  mars_3403[i,30]=mars_3403[i,6]-min(mars_3403[i,5:7])
  mars_3403[i,31]=mars_3403[i,7]-min(mars_3403[i,5:7])
  if(mars_3403[i,29]==0){mars_3403[i,32]=1}
  if(mars_3403[i,30]==0){mars_3403[i,32]=2}
  if(mars_3403[i,31]==0){mars_3403[i,32]=3}
  mars_3403[i,33]=length(na.omit(ts_comp[[i]][1,]))
  mars_3403[i,34]=length(na.omit(ts_comp[[i]][2,]))
  
  
  TS.plot.pdf(x=ts_comp[[i]],m=mars_3403$mod[i],sp=mars_3403$cName[i],GZ='Key Largo')
  dev.off()
  print(i)
}

setwd("C:/Users/14388/Desktop/Scripps - Project 1 RVC and REEF/RVC")

write.csv(mars_3403,'MARSS_model_comparisons_3403_1980_time_matched_expert_spur_groove.csv')

mars_3403<- mars_3403[complete.cases(mars_3403),]

mars_3403$AP<- NA


#####Antipolar trends####
mars_3403_sg_p<- read.csv('./MARSS Output comparisons/MARSS_model_comparisons_3403_1993_expert_sg_isol_cont_siteweight.csv')

mars_3403_sg_p<- mars_3403_sg_p[complete.cases(mars_3403_sg_p),]
for(i in 1:nrow(mars_3403_sg_p)){
  mars_3403_sg_p$AIC.W1[i]<- exp(-0.5*mars_3403_sg_p$dAIC1[i])/sum(c(exp(-0.5*mars_3403_sg_p$dAIC1[i]),exp(-0.5*mars_3403_sg_p$dAIC2[i]),exp(-0.5*mars_3403_sg_p$dAIC3[i])))
  mars_3403_sg_p$AIC.W2[i]<- exp(-0.5*mars_3403_sg_p$dAIC2[i])/sum(c(exp(-0.5*mars_3403_sg_p$dAIC1[i]),exp(-0.5*mars_3403_sg_p$dAIC2[i]),exp(-0.5*mars_3403_sg_p$dAIC3[i])))
  mars_3403_sg_p$AIC.W3[i]<- exp(-0.5*mars_3403_sg_p$dAIC3[i])/sum(c(exp(-0.5*mars_3403_sg_p$dAIC1[i]),exp(-0.5*mars_3403_sg_p$dAIC2[i]),exp(-0.5*mars_3403_sg_p$dAIC3[i])))
  
  if(mars_3403_sg_p$mod[i]==1){mars_3403_sg_p$AP[i]=0}else{
    
    if(abs(mars_3403_sg_p$U2.rvc[i])>=0.005&abs(mars_3403_sg_p$U2.reef[i])>=0.005){
      if(c(sign(mars_3403_sg_p$U2.rvc[i])+sign(mars_3403_sg_p$U2.reef[i]))==0){
        mars_3403_sg_p$AP[i]=1
      }
    }else{mars_3403_sg_p$AP[i]=0}
  }
  }

mars_3403_sg_p$AP
plot(AP~mAbund.rvc,data=mars_3403_sg_p)
plot(AP~mAbund.reef,data=mars_3403_sg_p)
plot(AP~sdAbund.rvc,data=mars_3403_sg_p)
plot(AP~sdAbund.reef,data=mars_3403_sg_p)

mars_3403_sg<- read.csv('./MARSS Output comparisons/MARSS_model_comparisons_3403_1993_expert_spurgroove_siteweight.csv')
mars_3403_sg<- mars_3403_sg[complete.cases(mars_3403_sg),]
for(i in 1:nrow(mars_3403_sg)){
  mars_3403_sg$AIC.W1[i]<- exp(-0.5*mars_3403_sg$dAIC1[i])/sum(c(exp(-0.5*mars_3403_sg$dAIC1[i]),exp(-0.5*mars_3403_sg$dAIC2[i]),exp(-0.5*mars_3403_sg$dAIC3[i])))
  mars_3403_sg$AIC.W2[i]<- exp(-0.5*mars_3403_sg$dAIC2[i])/sum(c(exp(-0.5*mars_3403_sg$dAIC1[i]),exp(-0.5*mars_3403_sg$dAIC2[i]),exp(-0.5*mars_3403_sg$dAIC3[i])))
  mars_3403_sg$AIC.W3[i]<- exp(-0.5*mars_3403_sg$dAIC3[i])/sum(c(exp(-0.5*mars_3403_sg$dAIC1[i]),exp(-0.5*mars_3403_sg$dAIC2[i]),exp(-0.5*mars_3403_sg$dAIC3[i])))
  
  if(mars_3403_sg$mod[i]==1){mars_3403_sg$AP[i]=0}else{
    
    if(abs(mars_3403_sg$U2.rvc[i])>=0.005&abs(mars_3403_sg$U2.reef[i])>=0.005){
      if(c(sign(mars_3403_sg$U2.rvc[i])+sign(mars_3403_sg$U2.reef[i]))==0){
        mars_3403_sg$AP[i]=1
      }
    }else{mars_3403_sg$AP[i]=0}
  }
}

mars_3403_sg$AP
plot(AP~mAbund.rvc,data=mars_3403_sg)
plot(AP~mAbund.reef,data=mars_3403_sg)
plot(AP~sdAbund.rvc,data=mars_3403_sg)
plot(AP~sdAbund.reef,data=mars_3403_sg)


#Differences
m<- match(mars_3403_sg$SP,mars_3403_sg_p$SP)

summary(as.factor(paste(mars_3403_sg$mod,mars_3403_sg_p$mod[m])))

###Ordinal effort filtering
library(ordinal)

test<- clmm2(as.factor(abundance)~1,data=barra,threshold='flexible')


####MARSS batches - 1993 to 2018, expert data####
####3403
REEF.3403.comb.green<- do.call(rbind, lapply(REEF_3403_green, data.frame, stringsAsFactors=FALSE))
REEF.green.3403.sp<- unique(REEF.3403.comb.green$sciName)

#REEF.3403.exp.comb.green<- do.call(rbind, lapply(REEF_3403_exp_green, data.frame, stringsAsFactors=FALSE))
#REEF.green.exp.3403.sp<- unique(REEF.3403.exp.comb.green$sciName)

rvc.3403.comb<- do.call(rbind, lapply(rvc_3403_green, data.frame, stringsAsFactors=FALSE))
rvc.green.3403.sp<- unique(rvc.3403.comb$Species)

ts_comp<- list()
years<- seq(1993,2018,by=1)
n<- NA
R_matrix<- list()





#Region 3404
REEF.3404.comb.green<- do.call(rbind, lapply(REEF_3404_green, data.frame, stringsAsFactors=FALSE))
REEF.green.3404.sp<- unique(REEF.3404.comb.green$comName)

REEF.3404.exp.comb.green<- do.call(rbind, lapply(REEF_3404_exp_green, data.frame, stringsAsFactors=FALSE))
REEF.green.exp.3404.sp<- unique(REEF.3404.exp.comb.green$sciName)

rvc.3404.comb<- do.call(rbind, lapply(rvc_3404_green, data.frame, stringsAsFactors=FALSE))
rvc.green.3404.sp<- unique(rvc.3404.comb$Species)

ts_comp<- list()
years<- seq(1993,2018,by=1)
n<- NA
R_matrix<- list()

setwd("C:/Users/14388/Desktop/Scripps - Project 1 RVC and REEF/RVC/timeseries/1993 expert timeseries - 3404")
mars_3404<- data.frame(SP=NA,conv.1=NA,conv.2=NA,conv.3=NA,tier.1.AICc=NA,tier2.AICc=NA,tier3.AICc=NA,Q1=NA,Q2=NA,Q3.rvc=NA,Q3.reef=NA,U1=NA,U2.rvc=NA,U2.reef=NA,U3.rvc=NA,U3.reef=NA,R1.rvc=NA,R1.reef=NA,R2.rvc=NA,R2.reef=NA,R3.rvc=NA,R3.reef=NA,cName=NA,obs.grp=NA,mAbund.rvc=NA,mAbund.reef=NA,sdAbund.rvc=NA,sdAbund.reef=NA,dAIC1=NA,dAIC2=NA,dAIC3=NA,mod=NA,rvc.years=NA,reef.years=NA)
for(i in 1:length(rvc.green.3404.sp)){
  spp<- filter(fish_reef,rvc_code==rvc.green.3404.sp[i])
  ts_comp[[i]]<- t(log(rvc_3404_green[[i]]$mean_ssu_abundance)) #extract logged time-series
  m<- match(spp$scientificname,REEF.green.exp.3404.sp)
  if(is.na(m)==T){next}
  ts_comp[[i]]<- rbind(ts_comp[[i]],t(REEF_3404_exp_green[[m]]$logAbund[1:26])) #Data from 1995 to 2018
  rownames(ts_comp[[i]])<- c(paste(spp$scientificname,'rvc',sep=":"),paste(spp$scientificname,'REEF',sep=":"))
  colnames(ts_comp[[i]])<- years
  
  Z_1=factor(c('rvc_reef','rvc_reef')) #Common state process for both time-series
  Z_2=factor(c('rvc','reef')) #Different state process for each time-series
  R_matrix[[i]]=matrix(list(0),nrow(ts_comp[[i]]),nrow(ts_comp[[i]])) #create empty Obs. error matrix
  diag(R_matrix[[i]])=c(rep("rvc",times=1),rep("reef",times=1))
  
  #tier 1 - common state process
  tier_1<- list(Z=Z_1,
                Q = 'diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 2 - Diff. U, same process error
  tier_2<- list(Z=Z_2,
                Q ='diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 3
  tier_3<- list(Z=Z_2,
                Q = 'diagonal and unequal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  fit_1<-  MARSS(ts_comp[[i]], model = tier_1,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  fit_2<-  MARSS(ts_comp[[i]], model = tier_2,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  fit_3<-  MARSS(ts_comp[[i]], model = tier_3,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  params.1<- MARSSparamCIs(fit_1)
  params.2<- MARSSparamCIs(fit_2)
  params.3<- MARSSparamCIs(fit_3)
  
  
  mars_3404[i,1]=spp$scientificname
  mars_3404[i,2]=ifelse(is.null(fit_1$errors)==T,1,0)
  mars_3404[i,3]=ifelse(is.null(fit_2$errors)==T,1,0)  
  mars_3404[i,4]=ifelse(is.null(fit_3$errors)==T,1,0)  
  mars_3404[i,5]=fit_1$AICc
  mars_3404[i,6]=fit_2$AICc  
  mars_3404[i,7]=fit_3$AICc  
  mars_3404[i,8]=params.1$parMean[5]
  mars_3404[i,9]=params.2$parMean[5]
  mars_3404[i,10]=params.3$parMean[5]
  mars_3404[i,11]=params.3$parMean[6]
  mars_3404[i,12]=params.1$parMean[4]
  mars_3404[i,13]=params.2$parMean[3] 
  mars_3404[i,14]=params.2$parMean[4]
  mars_3404[i,15]=params.3$parMean[3]
  mars_3404[i,16]=params.3$parMean[4]
  mars_3404[i,17]=params.1$parMean[2] 
  mars_3404[i,18]=params.1$parMean[3]
  mars_3404[i,19]=params.2$parMean[1] 
  mars_3404[i,20]=params.2$parMean[2]
  mars_3404[i,21]=params.3$parMean[1] 
  mars_3404[i,22]=params.3$parMean[2]
  mars_3404[i,23]=spp$commonname
  mars_3404[i,24]=spp$group_assignment
  mars_3404[i,25]=exp(mean(na.omit(ts_comp[[i]][1,])))
  mars_3404[i,26]=exp(mean(na.omit(ts_comp[[i]][2,])))
  mars_3404[i,27]=sd(exp(na.omit(ts_comp[[i]][1,])))
  mars_3404[i,28]=sd(exp(na.omit(ts_comp[[i]][2,])))
  mars_3404[i,29]=mars_3404[i,5]-min(mars_3404[i,5:7])
  mars_3404[i,30]=mars_3404[i,6]-min(mars_3404[i,5:7])
  mars_3404[i,31]=mars_3404[i,7]-min(mars_3404[i,5:7])
  if(mars_3404[i,29]==0){mars_3404[i,32]=1}
  if(mars_3404[i,30]==0){mars_3404[i,32]=2}
  if(mars_3404[i,31]==0){mars_3404[i,32]=3}
  mars_3404[i,33]=length(na.omit(ts_comp[[i]][1,]))
  mars_3404[i,34]=length(na.omit(ts_comp[[i]][2,]))
  
  
  TS.plot.pdf(x=ts_comp[[i]],m=mars_3404$mod[i],sp=mars_3404$cName[i],GZ='Islamorada')
  dev.off()
  
  print(i)
}
mars_3404<- mars_3404[complete.cases(mars_3404$SP),]
setwd("C:/Users/14388/Desktop/Scripps - Project 1 RVC and REEF/RVC")

write.csv(mars_3404,'MARSS_model_comparisons_expert_3404.csv')


#Region 3408
REEF.3408.comb.green<- do.call(rbind, lapply(REEF_3408_green, data.frame, stringsAsFactors=FALSE))
REEF.green.3408.sp<- unique(REEF.3408.comb.green$comName)

REEF.3408.exp.comb.green<- do.call(rbind, lapply(REEF_3408_exp_green, data.frame, stringsAsFactors=FALSE))
REEF.green.exp.3408.sp<- unique(REEF.3408.exp.comb.green$sciName)

rvc.3408.comb<- do.call(rbind, lapply(rvc_3408_green, data.frame, stringsAsFactors=FALSE))
rvc.green.3408.sp<- unique(rvc.3408.comb$Species)

ts_comp<- list()
years<- seq(1993,2018,by=1)
n<- NA
R_matrix<- list()
setwd("C:/Users/14388/Desktop/Scripps - Project 1 RVC and REEF/RVC/timeseries/1993 expert time series - 3408")

mars_3408<- data.frame(SP=NA,conv.1=NA,conv.2=NA,conv.3=NA,tier.1.AICc=NA,tier2.AICc=NA,tier3.AICc=NA,Q1=NA,Q2=NA,Q3.rvc=NA,Q3.reef=NA,U1=NA,U2.rvc=NA,U2.reef=NA,U3.rvc=NA,U3.reef=NA,R1.rvc=NA,R1.reef=NA,R2.rvc=NA,R2.reef=NA,R3.rvc=NA,R3.reef=NA,cName=NA,obs.grp=NA,mAbund.rvc=NA,mAbund.reef=NA,sdAbund.rvc=NA,sdAbund.reef=NA,dAIC1=NA,dAIC2=NA,dAIC3=NA,mod=NA,rvc.years=NA,reef.years=NA)
for(i in 1:length(rvc.green.3408.sp)){
  spp<- filter(fish_reef,rvc_code==rvc.green.3408.sp[i])
  ts_comp[[i]]<- t(log(rvc_3408_green[[i]]$mean_ssu_abundance)) #extract logged time-series
  m<- match(spp$scientificname,REEF.green.exp.3408.sp)
  if(is.na(m)==T){next}
  ts_comp[[i]]<- rbind(ts_comp[[i]],t(REEF_3408_exp_green[[m]]$logAbund[1:26])) #Data from 1995 to 2018
  rownames(ts_comp[[i]])<- c(paste(spp$scientificname,'rvc',sep=":"),paste(spp$scientificname,'REEF',sep=":"))
  colnames(ts_comp[[i]])<- years
  
  Z_1=factor(c('rvc_reef','rvc_reef')) #Common state process for both time-series
  Z_2=factor(c('rvc','reef')) #Different state process for each time-series
  R_matrix[[i]]=matrix(list(0),nrow(ts_comp[[i]]),nrow(ts_comp[[i]])) #create empty Obs. error matrix
  diag(R_matrix[[i]])=c(rep("rvc",times=1),rep("reef",times=1))
  
  #tier 1 - common state process
  tier_1<- list(Z=Z_1,
                Q = 'diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 2 - Diff. U, same process error
  tier_2<- list(Z=Z_2,
                Q ='diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 3
  tier_3<- list(Z=Z_2,
                Q = 'diagonal and unequal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  fit_1<-  MARSS(ts_comp[[i]], model = tier_1,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  fit_2<-  MARSS(ts_comp[[i]], model = tier_2,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  fit_3<-  MARSS(ts_comp[[i]], model = tier_3,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  params.1<- MARSSparamCIs(fit_1)
  params.2<- MARSSparamCIs(fit_2)
  params.3<- MARSSparamCIs(fit_3)
  
  
  mars_3408[i,1]=spp$scientificname
  mars_3408[i,2]=ifelse(is.null(fit_1$errors)==T,1,0)
  mars_3408[i,3]=ifelse(is.null(fit_2$errors)==T,1,0)  
  mars_3408[i,4]=ifelse(is.null(fit_3$errors)==T,1,0)  
  mars_3408[i,5]=fit_1$AICc
  mars_3408[i,6]=fit_2$AICc  
  mars_3408[i,7]=fit_3$AICc  
  mars_3408[i,8]=params.1$parMean[5]
  mars_3408[i,9]=params.2$parMean[5]
  mars_3408[i,10]=params.3$parMean[5]
  mars_3408[i,11]=params.3$parMean[6]
  mars_3408[i,12]=params.1$parMean[4]
  mars_3408[i,13]=params.2$parMean[3] 
  mars_3408[i,14]=params.2$parMean[4]
  mars_3408[i,15]=params.3$parMean[3]
  mars_3408[i,16]=params.3$parMean[4]
  mars_3408[i,17]=params.1$parMean[2] 
  mars_3408[i,18]=params.1$parMean[3]
  mars_3408[i,19]=params.2$parMean[1] 
  mars_3408[i,20]=params.2$parMean[2]
  mars_3408[i,21]=params.3$parMean[1] 
  mars_3408[i,22]=params.3$parMean[2]
  mars_3408[i,23]=spp$commonname
  mars_3408[i,24]=spp$group_assignment
  mars_3408[i,25]=exp(mean(na.omit(ts_comp[[i]][1,])))
  mars_3408[i,26]=exp(mean(na.omit(ts_comp[[i]][2,])))
  mars_3408[i,27]=sd(exp(na.omit(ts_comp[[i]][1,])))
  mars_3408[i,28]=sd(exp(na.omit(ts_comp[[i]][2,])))
  mars_3408[i,29]=mars_3408[i,5]-min(mars_3408[i,5:7])
  mars_3408[i,30]=mars_3408[i,6]-min(mars_3408[i,5:7])
  mars_3408[i,31]=mars_3408[i,7]-min(mars_3408[i,5:7])
  if(mars_3408[i,29]==0){mars_3408[i,32]=1}
  if(mars_3408[i,30]==0){mars_3408[i,32]=2}
  if(mars_3408[i,31]==0){mars_3408[i,32]=3}
  mars_3408[i,33]=length(na.omit(ts_comp[[i]][1,]))
  mars_3408[i,34]=length(na.omit(ts_comp[[i]][2,]))
  
  
  TS.plot.pdf(x=ts_comp[[i]],m=mars_3408$mod[i],sp=mars_3408$cName[i],GZ='Key West')
  dev.off()
  
  print(i)
}

setwd("C:/Users/14388/Desktop/Scripps - Project 1 RVC and REEF/RVC")

write.csv(mars_3408,'MARSS_model_comparisons_expert_3408.csv')


####MARSS batches - 1993 to 2018, expert+nov data####
####3403
REEF.3403.comb.green<- do.call(rbind, lapply(REEF_3403_green, data.frame, stringsAsFactors=FALSE))
REEF.green.3403.sp<- unique(REEF.3403.comb.green$sciName)


rvc.3403.comb<- do.call(rbind, lapply(rvc_3403_green, data.frame, stringsAsFactors=FALSE))
rvc.green.3403.sp<- unique(rvc.3403.comb$Species)

ts_comp<- list()
years<- seq(1993,2018,by=1)
n<- NA
R_matrix<- list()

setwd("C:/Users/14388/Desktop/Scripps - Project 1 RVC and REEF/RVC/timeseries/1993 nov+expert timeseries -3403")
mars_3403<- data.frame(SP=NA,conv.1=NA,conv.2=NA,conv.3=NA,tier.1.AICc=NA,tier2.AICc=NA,tier3.AICc=NA,Q1=NA,Q2=NA,Q3.rvc=NA,Q3.reef=NA,U1=NA,U2.rvc=NA,U2.reef=NA,U3.rvc=NA,U3.reef=NA,R1.rvc=NA,R1.reef=NA,R2.rvc=NA,R2.reef=NA,R3.rvc=NA,R3.reef=NA,cName=NA,obs.grp=NA,mAbund.rvc=NA,mAbund.reef=NA,sdAbund.rvc=NA,sdAbund.reef=NA,dAIC1=NA,dAIC2=NA,dAIC3=NA,mod=NA,years.rvc=NA,years.reef=NA)
for(i in 1:length(rvc.green.3403.sp)){
  spp<- filter(fish_reef,rvc_code==rvc.green.3403.sp[i])
  ts_comp[[i]]<- t(log(rvc_3403_green[[i]]$mean_ssu_abundance)) #extract logged time-series
  m<- match(spp$scientificname,REEF.green.3403.sp)
  if(is.na(m)==T){next}
  ts_comp[[i]]<- rbind(ts_comp[[i]],t(REEF_3403_green[[m]]$logAbund[1:26]))
  rownames(ts_comp[[i]])<- c(paste(spp$scientificname,'rvc',sep=":"),paste(spp$scientificname,'REEF',sep=":"))
  colnames(ts_comp[[i]])<- years
  
  Z_1=factor(c('rvc_reef','rvc_reef')) #Common state process for both time-series
  Z_2=factor(c('rvc','reef')) #Different state process for each time-series
  R_matrix[[i]]=matrix(list(0),nrow(ts_comp[[i]]),nrow(ts_comp[[i]])) #create empty Obs. error matrix
  diag(R_matrix[[i]])=c(rep("rvc",times=1),rep("reef",times=1))
  
  #tier 1 - common state process
  tier_1<- list(Z=Z_1,
                Q = 'diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 2 - Diff. U, same process error
  tier_2<- list(Z=Z_2,
                Q ='diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 3
  tier_3<- list(Z=Z_2,
                Q = 'diagonal and unequal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  fit_1<-  MARSS(ts_comp[[i]], model = tier_1,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  fit_2<-  MARSS(ts_comp[[i]], model = tier_2,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  fit_3<-  MARSS(ts_comp[[i]], model = tier_3,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  params.1<- MARSSparamCIs(fit_1)
  params.2<- MARSSparamCIs(fit_2)
  params.3<- MARSSparamCIs(fit_3)
  
  mars_3403[i,1]=spp$scientificname
  mars_3403[i,2]=ifelse(is.null(fit_1$errors)==T,1,0)
  mars_3403[i,3]=ifelse(is.null(fit_2$errors)==T,1,0)  
  mars_3403[i,4]=ifelse(is.null(fit_3$errors)==T,1,0)  
  mars_3403[i,5]=fit_1$AICc
  mars_3403[i,6]=fit_2$AICc  
  mars_3403[i,7]=fit_3$AICc  
  mars_3403[i,8]=params.1$parMean[5]
  mars_3403[i,9]=params.2$parMean[5]
  mars_3403[i,10]=params.3$parMean[5]
  mars_3403[i,11]=params.3$parMean[6]
  mars_3403[i,12]=params.1$parMean[4]
  mars_3403[i,13]=params.2$parMean[3] 
  mars_3403[i,14]=params.2$parMean[4]
  mars_3403[i,15]=params.3$parMean[3]
  mars_3403[i,16]=params.3$parMean[4]
  mars_3403[i,17]=params.1$parMean[2] 
  mars_3403[i,18]=params.1$parMean[3]
  mars_3403[i,19]=params.2$parMean[1] 
  mars_3403[i,20]=params.2$parMean[2]
  mars_3403[i,21]=params.3$parMean[1] 
  mars_3403[i,22]=params.3$parMean[2]
  mars_3403[i,23]=spp$commonname
  mars_3403[i,24]=spp$group_assignment
  mars_3403[i,25]=exp(mean(na.omit(ts_comp[[i]][1,])))
  mars_3403[i,26]=exp(mean(na.omit(ts_comp[[i]][2,])))
  mars_3403[i,27]=sd(exp(na.omit(ts_comp[[i]][1,])))
  mars_3403[i,28]=sd(exp(na.omit(ts_comp[[i]][2,])))
  mars_3403[i,29]=mars_3403[i,5]-min(mars_3403[i,5:7])
  mars_3403[i,30]=mars_3403[i,6]-min(mars_3403[i,5:7])
  mars_3403[i,31]=mars_3403[i,7]-min(mars_3403[i,5:7])
  if(mars_3403[i,29]==0){mars_3403[i,32]=1}
  if(mars_3403[i,30]==0){mars_3403[i,32]=2}
  if(mars_3403[i,31]==0){mars_3403[i,32]=3}
  mars_3403[i,33]=length(na.omit(ts_comp[[i]][1,]))
  mars_3403[i,34]=length(na.omit(ts_comp[[i]][2,]))
  
#  TS.plot.pdf(x=ts_comp[[i]],m=mars_3403$mod[i],sp=mars_3403$cName[i],GZ='Key Largo')
#  dev.off()
  print(i)
}

setwd("C:/Users/14388/Desktop/Scripps - Project 1 RVC and REEF/RVC")

write.csv(mars_3403,'MARSS_model_comparisons_3403_1993_nov+expert_sg.csv')

####no trend####
mars_3403_nt<- data.frame(SP=NA,conv.1=NA,conv.2=NA,tier.1.AICc=NA,tier2.AICc=NA,Q1=NA,Q2.rvc=NA,Q2.reef=NA,R1.rvc=NA,R1.reef=NA,R2.rvc=NA,R2.reef=NA,cName=NA,obs.grp=NA,mAbund.rvc=NA,mAbund.reef=NA,sdAbund.rvc=NA,sdAbund.reef=NA,dAIC1=NA,dAIC2=NA,mod=NA,years.rvc=NA,years.reef=NA)
for(i in 1:length(rvc.green.3403.sp)){
  spp<- filter(fish_reef,rvc_code==rvc.green.3403.sp[i])
  ts_comp[[i]]<- t(log(rvc_3403_green[[i]]$mean_ssu_abundance)) #extract logged time-series
  m<- match(spp$scientificname,REEF.green.3403.sp)
  if(is.na(m)==T){next}
  ts_comp[[i]]<- rbind(ts_comp[[i]],t(REEF_3403_green[[m]]$logAbund[1:26]))
  rownames(ts_comp[[i]])<- c(paste(spp$scientificname,'rvc',sep=":"),paste(spp$scientificname,'REEF',sep=":"))
  colnames(ts_comp[[i]])<- years
  
  Z_1=factor(c('rvc_reef','rvc_reef')) #Common state process for both time-series
  Z_2=factor(c('rvc','reef')) #Different state process for each time-series
  R_matrix[[i]]=matrix(list(0),nrow(ts_comp[[i]]),nrow(ts_comp[[i]])) #create empty Obs. error matrix
  diag(R_matrix[[i]])=c(rep("rvc",times=1),rep("reef",times=1))
  
  #tier 1 - common state process
  tier_1<- list(Z=Z_1,
                Q = 'diagonal and equal',
                U = 'zero',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 2 - Diff. U, same process error
  tier_2<- list(Z=Z_2,
                Q ='diagonal and unequal',
                R = R_matrix[[i]],
                A = 'scaling',
                U = 'zero',
                x0 = 'unequal',
                tinitx=0)

  
  fit_1<-  MARSS(ts_comp[[i]], model = tier_1,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  fit_2<-  MARSS(ts_comp[[i]], model = tier_2,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  params.1<- MARSSparamCIs(fit_1)
  params.2<- MARSSparamCIs(fit_2)

  mars_3403_nt[i,1]=spp$scientificname
  mars_3403_nt[i,2]=ifelse(is.null(fit_1$errors)==T,1,0)
  mars_3403_nt[i,3]=ifelse(is.null(fit_2$errors)==T,1,0)  
  mars_3403_nt[i,4]=fit_1$AICc
  mars_3403_nt[i,5]=fit_2$AICc  
  mars_3403_nt[i,6]=params.1$parMean[4]
  mars_3403_nt[i,7]=params.2$parMean[3]
  mars_3403_nt[i,8]=params.2$parMean[4]
  mars_3403_nt[i,9]=params.1$parMean[2]
  mars_3403_nt[i,10]=params.1$parMean[3]
  mars_3403_nt[i,11]=params.2$parMean[1]
  mars_3403_nt[i,12]=params.2$parMean[2]
  mars_3403_nt[i,13]=spp$commonname
  mars_3403_nt[i,14]=spp$group_assignment
  mars_3403_nt[i,15]=exp(mean(na.omit(ts_comp[[i]][1,])))
  mars_3403_nt[i,16]=exp(mean(na.omit(ts_comp[[i]][2,])))
  mars_3403_nt[i,17]=sd(exp(na.omit(ts_comp[[i]][1,])))
  mars_3403_nt[i,18]=sd(exp(na.omit(ts_comp[[i]][2,])))
  mars_3403_nt[i,19]=mars_3403_nt[i,4]-min(mars_3403_nt[i,4:5])
  mars_3403_nt[i,20]=mars_3403_nt[i,5]-min(mars_3403_nt[i,4:5])
  if(mars_3403_nt[i,19]==0){mars_3403_nt[i,21]=1}
  if(mars_3403_nt[i,20]==0){mars_3403_nt[i,21]=2}
  mars_3403_nt[i,22]=length(na.omit(ts_comp[[i]][1,]))
  mars_3403_nt[i,23]=length(na.omit(ts_comp[[i]][2,]))
  
  #  TS.plot.pdf(x=ts_comp[[i]],m=mars_3403$mod[i],sp=mars_3403$cName[i],GZ='Key Largo')
  #  dev.off()
  print(i)
}

write.csv(mars_3403,'MARSS_model_comparisons_3403_1993_notrend.csv')



#Region 3404
REEF.3404.comb.green<- do.call(rbind, lapply(REEF_3404_green, data.frame, stringsAsFactors=FALSE))
REEF.green.3404.sp<- unique(REEF.3404.comb.green$sciName)

rvc.3404.comb<- do.call(rbind, lapply(rvc_3404_green, data.frame, stringsAsFactors=FALSE))
rvc.green.3404.sp<- unique(rvc.3404.comb$Species)

ts_comp<- list()
years<- seq(1993,2018,by=1)
n<- NA
R_matrix<- list()

setwd("C:/Users/14388/Desktop/Scripps - Project 1 RVC and REEF/RVC/timeseries/1993 nov+expert timeseries - 3404")
mars_3404<- data.frame(SP=NA,conv.1=NA,conv.2=NA,conv.3=NA,tier.1.AICc=NA,tier2.AICc=NA,tier3.AICc=NA,Q1=NA,Q2=NA,Q3.rvc=NA,Q3.reef=NA,U1=NA,U2.rvc=NA,U2.reef=NA,U3.rvc=NA,U3.reef=NA,R1.rvc=NA,R1.reef=NA,R2.rvc=NA,R2.reef=NA,R3.rvc=NA,R3.reef=NA,cName=NA,obs.grp=NA,mAbund.rvc=NA,mAbund.reef=NA,sdAbund.rvc=NA,sdAbund.reef=NA,dAIC1=NA,dAIC2=NA,dAIC3=NA,mod=NA,rvc.years=NA,reef.years=NA)
for(i in 1:length(rvc.green.3404.sp)){
  spp<- filter(fish_reef,rvc_code==rvc.green.3404.sp[i])
  ts_comp[[i]]<- t(log(rvc_3404_green[[i]]$mean_ssu_abundance)) #extract logged time-series
  m<- match(spp$scientificname,REEF.green.3404.sp)
  if(is.na(m)==T){next}
  ts_comp[[i]]<- rbind(ts_comp[[i]],t(REEF_3404_green[[m]]$logAbund[1:26])) #Data from 1995 to 2018
  rownames(ts_comp[[i]])<- c(paste(spp$scientificname,'rvc',sep=":"),paste(spp$scientificname,'REEF',sep=":"))
  colnames(ts_comp[[i]])<- years
  
  Z_1=factor(c('rvc_reef','rvc_reef')) #Common state process for both time-series
  Z_2=factor(c('rvc','reef')) #Different state process for each time-series
  R_matrix[[i]]=matrix(list(0),nrow(ts_comp[[i]]),nrow(ts_comp[[i]])) #create empty Obs. error matrix
  diag(R_matrix[[i]])=c(rep("rvc",times=1),rep("reef",times=1))
  
  #tier 1 - common state process
  tier_1<- list(Z=Z_1,
                Q = 'diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 2 - Diff. U, same process error
  tier_2<- list(Z=Z_2,
                Q ='diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 3
  tier_3<- list(Z=Z_2,
                Q = 'diagonal and unequal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  fit_1<-  MARSS(ts_comp[[i]], model = tier_1,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  fit_2<-  MARSS(ts_comp[[i]], model = tier_2,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  fit_3<-  MARSS(ts_comp[[i]], model = tier_3,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  params.1<- MARSSparamCIs(fit_1)
  params.2<- MARSSparamCIs(fit_2)
  params.3<- MARSSparamCIs(fit_3)
  
  
  mars_3404[i,1]=spp$scientificname
  mars_3404[i,2]=ifelse(is.null(fit_1$errors)==T,1,0)
  mars_3404[i,3]=ifelse(is.null(fit_2$errors)==T,1,0)  
  mars_3404[i,4]=ifelse(is.null(fit_3$errors)==T,1,0)  
  mars_3404[i,5]=fit_1$AICc
  mars_3404[i,6]=fit_2$AICc  
  mars_3404[i,7]=fit_3$AICc  
  mars_3404[i,8]=params.1$parMean[5]
  mars_3404[i,9]=params.2$parMean[5]
  mars_3404[i,10]=params.3$parMean[5]
  mars_3404[i,11]=params.3$parMean[6]
  mars_3404[i,12]=params.1$parMean[4]
  mars_3404[i,13]=params.2$parMean[3] 
  mars_3404[i,14]=params.2$parMean[4]
  mars_3404[i,15]=params.3$parMean[3]
  mars_3404[i,16]=params.3$parMean[4]
  mars_3404[i,17]=params.1$parMean[2] 
  mars_3404[i,18]=params.1$parMean[3]
  mars_3404[i,19]=params.2$parMean[1] 
  mars_3404[i,20]=params.2$parMean[2]
  mars_3404[i,21]=params.3$parMean[1] 
  mars_3404[i,22]=params.3$parMean[2]
  mars_3404[i,23]=spp$commonname
  mars_3404[i,24]=spp$group_assignment
  mars_3404[i,25]=exp(mean(na.omit(ts_comp[[i]][1,])))
  mars_3404[i,26]=exp(mean(na.omit(ts_comp[[i]][2,])))
  mars_3404[i,27]=sd(exp(na.omit(ts_comp[[i]][1,])))
  mars_3404[i,28]=sd(exp(na.omit(ts_comp[[i]][2,])))
  mars_3404[i,29]=mars_3404[i,5]-min(mars_3404[i,5:7])
  mars_3404[i,30]=mars_3404[i,6]-min(mars_3404[i,5:7])
  mars_3404[i,31]=mars_3404[i,7]-min(mars_3404[i,5:7])
  if(mars_3404[i,29]==0){mars_3404[i,32]=1}
  if(mars_3404[i,30]==0){mars_3404[i,32]=2}
  if(mars_3404[i,31]==0){mars_3404[i,32]=3}
  mars_3404[i,33]=length(na.omit(ts_comp[[i]][1,]))
  mars_3404[i,34]=length(na.omit(ts_comp[[i]][2,]))
  
  
  TS.plot.pdf(x=ts_comp[[i]],m=mars_3404$mod[i],sp=mars_3404$cName[i],GZ='Islamorada')
  dev.off()
  
  print(i)
}
mars_3404<- mars_3404[complete.cases(mars_3404$SP),]
setwd("C:/Users/14388/Desktop/Scripps - Project 1 RVC and REEF/RVC")

write.csv(mars_3404,'MARSS_model_comparisons_nov_expert_3404.csv')


#Region 3408
REEF.3408.comb.green<- do.call(rbind, lapply(REEF_3408_green, data.frame, stringsAsFactors=FALSE))
REEF.green.3408.sp<- unique(REEF.3408.comb.green$comName)

rvc.3408.comb<- do.call(rbind, lapply(rvc_3408_green, data.frame, stringsAsFactors=FALSE))
rvc.green.3408.sp<- unique(rvc.3408.comb$Species)

ts_comp<- list()
years<- seq(1993,2018,by=1)
n<- NA
R_matrix<- list()
setwd("C:/Users/14388/Desktop/Scripps - Project 1 RVC and REEF/RVC/timeseries/1993 expert time series - 3408")

mars_3408<- data.frame(SP=NA,conv.1=NA,conv.2=NA,conv.3=NA,tier.1.AICc=NA,tier2.AICc=NA,tier3.AICc=NA,Q1=NA,Q2=NA,Q3.rvc=NA,Q3.reef=NA,U1=NA,U2.rvc=NA,U2.reef=NA,U3.rvc=NA,U3.reef=NA,R1.rvc=NA,R1.reef=NA,R2.rvc=NA,R2.reef=NA,R3.rvc=NA,R3.reef=NA,cName=NA,obs.grp=NA,mAbund.rvc=NA,mAbund.reef=NA,sdAbund.rvc=NA,sdAbund.reef=NA,dAIC1=NA,dAIC2=NA,dAIC3=NA,mod=NA,rvc.years=NA,reef.years=NA)
for(i in 1:length(rvc.green.3408.sp)){
  spp<- filter(fish_reef,rvc_code==rvc.green.3408.sp[i])
  ts_comp[[i]]<- t(log(rvc_3408_green[[i]]$mean_ssu_abundance)) #extract logged time-series
  m<- match(spp$scientificname,REEF.green.3408.sp)
  if(is.na(m)==T){next}
  ts_comp[[i]]<- rbind(ts_comp[[i]],t(REEF_3408_green[[m]]$logAbund[1:26])) #Data from 1995 to 2018
  rownames(ts_comp[[i]])<- c(paste(spp$scientificname,'rvc',sep=":"),paste(spp$scientificname,'REEF',sep=":"))
  colnames(ts_comp[[i]])<- years
  
  Z_1=factor(c('rvc_reef','rvc_reef')) #Common state process for both time-series
  Z_2=factor(c('rvc','reef')) #Different state process for each time-series
  R_matrix[[i]]=matrix(list(0),nrow(ts_comp[[i]]),nrow(ts_comp[[i]])) #create empty Obs. error matrix
  diag(R_matrix[[i]])=c(rep("rvc",times=1),rep("reef",times=1))
  
  #tier 1 - common state process
  tier_1<- list(Z=Z_1,
                Q = 'diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 2 - Diff. U, same process error
  tier_2<- list(Z=Z_2,
                Q ='diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 3
  tier_3<- list(Z=Z_2,
                Q = 'diagonal and unequal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  fit_1<-  MARSS(ts_comp[[i]], model = tier_1,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  fit_2<-  MARSS(ts_comp[[i]], model = tier_2,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  fit_3<-  MARSS(ts_comp[[i]], model = tier_3,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  params.1<- MARSSparamCIs(fit_1)
  params.2<- MARSSparamCIs(fit_2)
  params.3<- MARSSparamCIs(fit_3)
  
  
  mars_3408[i,1]=spp$scientificname
  mars_3408[i,2]=ifelse(is.null(fit_1$errors)==T,1,0)
  mars_3408[i,3]=ifelse(is.null(fit_2$errors)==T,1,0)  
  mars_3408[i,4]=ifelse(is.null(fit_3$errors)==T,1,0)  
  mars_3408[i,5]=fit_1$AICc
  mars_3408[i,6]=fit_2$AICc  
  mars_3408[i,7]=fit_3$AICc  
  mars_3408[i,8]=params.1$parMean[5]
  mars_3408[i,9]=params.2$parMean[5]
  mars_3408[i,10]=params.3$parMean[5]
  mars_3408[i,11]=params.3$parMean[6]
  mars_3408[i,12]=params.1$parMean[4]
  mars_3408[i,13]=params.2$parMean[3] 
  mars_3408[i,14]=params.2$parMean[4]
  mars_3408[i,15]=params.3$parMean[3]
  mars_3408[i,16]=params.3$parMean[4]
  mars_3408[i,17]=params.1$parMean[2] 
  mars_3408[i,18]=params.1$parMean[3]
  mars_3408[i,19]=params.2$parMean[1] 
  mars_3408[i,20]=params.2$parMean[2]
  mars_3408[i,21]=params.3$parMean[1] 
  mars_3408[i,22]=params.3$parMean[2]
  mars_3408[i,23]=spp$commonname
  mars_3408[i,24]=spp$group_assignment
  mars_3408[i,25]=exp(mean(na.omit(ts_comp[[i]][1,])))
  mars_3408[i,26]=exp(mean(na.omit(ts_comp[[i]][2,])))
  mars_3408[i,27]=sd(exp(na.omit(ts_comp[[i]][1,])))
  mars_3408[i,28]=sd(exp(na.omit(ts_comp[[i]][2,])))
  mars_3408[i,29]=mars_3408[i,5]-min(mars_3408[i,5:7])
  mars_3408[i,30]=mars_3408[i,6]-min(mars_3408[i,5:7])
  mars_3408[i,31]=mars_3408[i,7]-min(mars_3408[i,5:7])
  if(mars_3408[i,29]==0){mars_3408[i,32]=1}
  if(mars_3408[i,30]==0){mars_3408[i,32]=2}
  if(mars_3408[i,31]==0){mars_3408[i,32]=3}
  mars_3408[i,33]=length(na.omit(ts_comp[[i]][1,]))
  mars_3408[i,34]=length(na.omit(ts_comp[[i]][2,]))
  
  
  TS.plot.pdf(x=ts_comp[[i]],m=mars_3408$mod[i],sp=mars_3408$cName[i],GZ='Key West')
  dev.off()
  
  print(i)
}

setwd("C:/Users/14388/Desktop/Scripps - Project 1 RVC and REEF/RVC")

write.csv(mars_3408,'MARSS_model_comparisons_nov_expert_3408.csv')


#### Full data vs. Spur & Groove ####
mars_3403_sg<- read.csv('MARSS_model_comparisons_3403_1993_nov+expert_sg.csv')
mars_3403<-  read.csv('MARSS_model_comparisons_3403_1993_nov+expert.csv')

m<- match(mars_3403_sg$SP,mars_3403$SP)

compare<- data.frame(mars_3403_sg$cName,mars_3403_sg$mod,mars_3403$mod[m])






####Sighting frequency time-series####

#3403
REEF.3403.comb.green<- do.call(rbind, lapply(REEF_3403_green, data.frame, stringsAsFactors=FALSE))
REEF.green.3403.sp<- unique(REEF.3403.comb.green$sciName)

rvc.3403.comb<- do.call(rbind, lapply(rvc_3403_green, data.frame, stringsAsFactors=FALSE))
rvc.green.3403.sp<- unique(rvc.3403.comb$Species)

ts_comp<- list()
years<- seq(1993,2018,by=1)
n<- NA
R_matrix<- list()


setwd("C:/Users/14388/Desktop/Scripps - Project 1 RVC and REEF/RVC/timeseries/1993 occ timeseries - 3403")
mars_3403<- data.frame(SP=NA,conv.1=NA,conv.2=NA,conv.3=NA,tier.1.AICc=NA,tier2.AICc=NA,tier3.AICc=NA,Q1=NA,Q2=NA,Q3.rvc=NA,Q3.reef=NA,U1=NA,U2.rvc=NA,U2.reef=NA,U3.rvc=NA,U3.reef=NA,R1.rvc=NA,R1.reef=NA,R2.rvc=NA,R2.reef=NA,R3.rvc=NA,R3.reef=NA,cName=NA,obs.grp=NA,mAbund.rvc=NA,mAbund.reef=NA,sdAbund.rvc=NA,sdAbund.reef=NA,dAIC1=NA,dAIC2=NA,dAIC3=NA,mod=NA,years.rvc=NA,years.reef=NA)
for(i in 1:length(rvc.green.3403.sp)){
  spp<- filter(fish_reef,rvc_code==rvc.green.3403.sp[i])
  ts_comp[[i]]<- t(log(c(rvc_3403_green[[i]]$ssu.occ/rvc_3403_green[[i]]$nm))) #extract logged time-series
  m<- match(spp$scientificname,REEF.green.3403.sp)
  if(is.na(m)==T){next}
  ts_comp[[i]]<- rbind(ts_comp[[i]],t(log(c(REEF_3403_green[[m]]$n.occ[1:26]/REEF_3403_green[[m]]$tot.surveys[1:26]))))
  rownames(ts_comp[[i]])<- c(paste(spp$scientificname,'rvc',sep=":"),paste(spp$scientificname,'REEF',sep=":"))
  colnames(ts_comp[[i]])<- years
  
  Z_1=factor(c('rvc_reef','rvc_reef')) #Common state process for both time-series
  Z_2=factor(c('rvc','reef')) #Different state process for each time-series
  R_matrix[[i]]=matrix(list(0),nrow(ts_comp[[i]]),nrow(ts_comp[[i]])) #create empty Obs. error matrix
  diag(R_matrix[[i]])=c(rep("rvc",times=1),rep("reef",times=1))
  
  #tier 1 - common state process
  tier_1<- list(Z=Z_1,
                Q = 'diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 2 - Diff. U, same process error
  tier_2<- list(Z=Z_2,
                Q ='diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 3
  tier_3<- list(Z=Z_2,
                Q = 'diagonal and unequal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  fit_1<-  MARSS(ts_comp[[i]], model = tier_1,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  fit_2<-  MARSS(ts_comp[[i]], model = tier_2,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  fit_3<-  MARSS(ts_comp[[i]], model = tier_3,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  params.1<- MARSSparamCIs(fit_1)
  params.2<- MARSSparamCIs(fit_2)
  params.3<- MARSSparamCIs(fit_3)
  
  mars_3403[i,1]=spp$scientificname
  mars_3403[i,2]=ifelse(is.null(fit_1$errors)==T,1,0)
  mars_3403[i,3]=ifelse(is.null(fit_2$errors)==T,1,0)  
  mars_3403[i,4]=ifelse(is.null(fit_3$errors)==T,1,0)  
  mars_3403[i,5]=fit_1$AICc
  mars_3403[i,6]=fit_2$AICc  
  mars_3403[i,7]=fit_3$AICc  
  mars_3403[i,8]=params.1$parMean[5]
  mars_3403[i,9]=params.2$parMean[5]
  mars_3403[i,10]=params.3$parMean[5]
  mars_3403[i,11]=params.3$parMean[6]
  mars_3403[i,12]=params.1$parMean[4]
  mars_3403[i,13]=params.2$parMean[3] 
  mars_3403[i,14]=params.2$parMean[4]
  mars_3403[i,15]=params.3$parMean[3]
  mars_3403[i,16]=params.3$parMean[4]
  mars_3403[i,17]=params.1$parMean[2] 
  mars_3403[i,18]=params.1$parMean[3]
  mars_3403[i,19]=params.2$parMean[1] 
  mars_3403[i,20]=params.2$parMean[2]
  mars_3403[i,21]=params.3$parMean[1] 
  mars_3403[i,22]=params.3$parMean[2]
  mars_3403[i,23]=spp$commonname
  mars_3403[i,24]=spp$group_assignment
  mars_3403[i,25]=exp(mean(na.omit(ts_comp[[i]][1,])))
  mars_3403[i,26]=exp(mean(na.omit(ts_comp[[i]][2,])))
  mars_3403[i,27]=sd(exp(na.omit(ts_comp[[i]][1,])))
  mars_3403[i,28]=sd(exp(na.omit(ts_comp[[i]][2,])))
  mars_3403[i,29]=mars_3403[i,5]-min(mars_3403[i,5:7])
  mars_3403[i,30]=mars_3403[i,6]-min(mars_3403[i,5:7])
  mars_3403[i,31]=mars_3403[i,7]-min(mars_3403[i,5:7])
  if(mars_3403[i,29]==0){mars_3403[i,32]=1}
  if(mars_3403[i,30]==0){mars_3403[i,32]=2}
  if(mars_3403[i,31]==0){mars_3403[i,32]=3}
  mars_3403[i,33]=length(na.omit(ts_comp[[i]][1,]))
  mars_3403[i,34]=length(na.omit(ts_comp[[i]][2,]))
  
  
  TS.plot.pdf(x=ts_comp[[i]],m=mars_3403$mod[i],sp=mars_3403$cName[i],GZ='Key Largo')
  dev.off()
  print(i)
}

setwd("C:/Users/14388/Desktop/Scripps - Project 1 RVC and REEF/RVC")

write.csv(mars_3403,'MARSS_model_comparisons_3403_1993_sighting_occurrences.csv')








#### Analysis of RVC data to 1980 ####
fk_79_98<- read.csv("Greenberg_FL_Pre1999.txt") #File courtesy of J. Blondeau (NOAA)
fk_79_98$len<- as.numeric(ifelse(fk_79_98$len==-9,0,fk_79_98$len)) #Replace -9 with 0

t<- fk_79_98 %>% group_by(station_nr,psu,YEAR,spcode) %>%  #Make old data synonymous with 99 - 2018 data
  summarise(NUM=n(),sum=sum(len),id=ID[1]) %>% arrange(psu,YEAR,station_nr,spcode) #Calculate number of sightings (each row is a sighting)
t$NUM<- ifelse(t$sum==0,0,t$NUM)

#Put back NUM into original data-frame, now each row represents a SSU survey
m<- match(fk_79_98$ID,t$id)
fk_79_98$NUM<- t$NUM[m]
fk_79_98<- fk_79_98[complete.cases(fk_79_98$NUM),]

##3.c Pre-1999 RVC geography ###
fk_79_98$site_ll<- paste(fk_79_98$lat_deg,fk_79_98$lon_deg,sep="_")
fk_79_98_sites<- distinct(fk_79_98,site_ll,.keep_all = T)
set4<- SpatialPoints(cbind(fk_79_98_sites$lat_deg,fk_79_98_sites$lon_deg))
matched_reef_site3<- apply(gDistance(set4, set2, byid=TRUE), 2, which.min)
fk_79_98_sites$REEF_site <- reef_geog$geog[matched_reef_site3]
fk_79_98_sites$REEF_lat<- reef_geog$lat.dd[matched_reef_site3]
fk_79_98_sites$REEF_lon<- reef_geog$lon.dd[matched_reef_site3]
fk_79_98_sites$region.id<- reef_geog$region.id[matched_reef_site3]

fk_79_98[,22:25]<- fk_79_98_sites[match(fk_79_98$site_ll,fk_79_98_sites$site_ll),22:25]

#Synonymize column names
colnames(fk_79_98)<- c(colnames(fk_99_18)[1],colnames(fk_99_18)[3:5],colnames(fk_99_18)[2],colnames(fk_99_18)[6],colnames(fk_99_18)[7:18],colnames(fk_99_18)[20],colnames(fk_99_18)[19],'site_ll',colnames(fk_99_18)[24:27]) 
#Join together so RVC data covers from 1993 to 2018
fk_79_18<- full_join(fk_99_18,fk_79_98)

#Trim down RVC dataset to these species
fk_79_18_trim<- subset(fk_79_18,SPECIES_CD %in% fish_reef$rvc_code)

#Separate out RVC data into the REEF sub-regions
rvc_3403<- subset(fk_79_18_trim,region.id=='3403')
#rvc_3404<- subset(fk_79_18_trim,region.id=='3404')
#rvc_3408<- subset(fk_79_18_trim,region.id=='3408')

#
rvc_trends_3403<- rvc_batch(rvc_3403)
#rvc_trends_3404<- rvc_batch(rvc_3404)
#rvc_trends_3408<- rvc_batch(rvc_3408)

#filter to species seen every year
rvc_3403_green<-  list.filter(rvc_trends_3403,length(na.omit(mean_ssu_abundance))>=24) #71 species
#rvc_3404_green<-  list.filter(rvc_trends_3404,length(na.omit(mean_ssu_abundance))>=28) #0 species
#rvc_3408_green<-  list.filter(rvc_trends_3408,length(na.omit(mean_ssu_abundance))>=30) #0 species


#Region 3403 - expert data
REEF.3403.exp.comb.green<- do.call(rbind, lapply(REEF_3403_exp_green, data.frame, stringsAsFactors=FALSE))
REEF.green.exp.3403.sp<- unique(REEF.3403.exp.comb.green$sciName)

rvc.3403.comb<- do.call(rbind, lapply(rvc_3403_green, data.frame, stringsAsFactors=FALSE))
rvc.green.3403.sp<- unique(rvc.3403.comb$Species)

ts_comp<- list()
years<- seq(1980,2018,by=1)
n<- NA
R_matrix<- list()

setwd("C:/Users/14388/Desktop/Scripps - Project 1 RVC and REEF/RVC/timeseries/1980 timeseries - expert")
mars_3403_1980<- data.frame(SP=NA,conv.1=NA,conv.2=NA,conv.3=NA,tier.1.AICc=NA,tier2.AICc=NA,tier3.AICc=NA,Q1=NA,Q2=NA,Q3.rvc=NA,Q3.reef=NA,U1=NA,U2.rvc=NA,U2.reef=NA,U3.rvc=NA,U3.reef=NA,R1.rvc=NA,R1.reef=NA,R2.rvc=NA,R2.reef=NA,R3.rvc=NA,R3.reef=NA,cName=NA,obs.grp=NA,mAbund.rvc=NA,mAbund.reef=NA,sdAbund.rvc=NA,sdAbund.reef=NA,dAIC1=NA,dAIC2=NA,dAIC3=NA,mod=NA,rvc.years=NA,reef.years=NA)
for(i in 1:length(rvc.green.3403.sp)){
  spp<- filter(fish_reef,rvc_code==rvc.green.3403.sp[i])
  ts_comp[[i]]<- t(log(rvc_3403_green[[i]]$mean_ssu_abundance[2:40])) #extract logged time-series
  m<- match(spp$scientificname,REEF.green.exp.3403.sp)
  if(is.na(m)==T){next}
  ts_comp[[i]]<- rbind(ts_comp[[i]],t(c(rep(NA,13),REEF_3403_exp_green[[m]]$logAbund[1:26])))
  rownames(ts_comp[[i]])<- c(paste(spp$scientificname,'rvc',sep=":"),paste(spp$scientificname,'REEF',sep=":"))
  colnames(ts_comp[[i]])<- years
  
  Z_1=factor(c('rvc_reef','rvc_reef')) #Common state process for both time-series
  Z_2=factor(c('rvc','reef')) #Different state process for each time-series
  R_matrix[[i]]=matrix(list(0),nrow(ts_comp[[i]]),nrow(ts_comp[[i]])) #create empty Obs. error matrix
  diag(R_matrix[[i]])=c(rep("rvc",times=1),rep("reef",times=1))
  
  #tier 1 - common state process
  tier_1<- list(Z=Z_1,
                Q = 'diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 2 - Diff. U, same process error
  tier_2<- list(Z=Z_2,
                Q ='diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 3
  tier_3<- list(Z=Z_2,
                Q = 'diagonal and unequal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
    
    fit_1<-  MARSS(ts_comp[[i]], model = tier_1,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem') 
    fit_2<-  MARSS(ts_comp[[i]], model = tier_2,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem')
    fit_3<-  MARSS(ts_comp[[i]], model = tier_3,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem')
    params.1<- MARSSparamCIs(fit_1)
    params.2<- MARSSparamCIs(fit_2)
    params.3<- MARSSparamCIs(fit_3)
    
    mars_3403_1980[i,1]=spp$scientificname
    mars_3403_1980[i,2]=ifelse(is.null(fit_1$errors)==T,1,0)
    mars_3403_1980[i,3]=ifelse(is.null(fit_2$errors)==T,1,0)  
    mars_3403_1980[i,4]=ifelse(is.null(fit_3$errors)==T,1,0)  
    mars_3403_1980[i,5]=fit_1$AICc
    mars_3403_1980[i,6]=fit_2$AICc  
    mars_3403_1980[i,7]=fit_3$AICc  
    mars_3403_1980[i,8]=params.1$parMean[5]
    mars_3403_1980[i,9]=params.2$parMean[5]
    mars_3403_1980[i,10]=params.3$parMean[5]
    mars_3403_1980[i,11]=params.3$parMean[6]
    mars_3403_1980[i,12]=params.1$parMean[4]
    mars_3403_1980[i,13]=params.2$parMean[3] 
    mars_3403_1980[i,14]=params.2$parMean[4]
    mars_3403_1980[i,15]=params.3$parMean[3]
    mars_3403_1980[i,16]=params.3$parMean[4]
    mars_3403_1980[i,17]=params.1$parMean[2] 
    mars_3403_1980[i,18]=params.1$parMean[3]
    mars_3403_1980[i,19]=params.2$parMean[1] 
    mars_3403_1980[i,20]=params.2$parMean[2]
    mars_3403_1980[i,21]=params.3$parMean[1] 
    mars_3403_1980[i,22]=params.3$parMean[2]
    mars_3403_1980[i,23]=spp$commonname
    mars_3403_1980[i,24]=spp$group_assignment
    mars_3403_1980[i,25]=exp(mean(na.omit(ts_comp[[i]][1,])))
    mars_3403_1980[i,26]=exp(mean(na.omit(ts_comp[[i]][2,])))
    mars_3403_1980[i,27]=sd(exp(na.omit(ts_comp[[i]][1,])))
    mars_3403_1980[i,28]=sd(exp(na.omit(ts_comp[[i]][2,])))
    mars_3403_1980[i,29]=mars_3403_1980[i,5]-min(mars_3403_1980[i,5:7])
    mars_3403_1980[i,30]=mars_3403_1980[i,6]-min(mars_3403_1980[i,5:7])
    mars_3403_1980[i,31]=mars_3403_1980[i,7]-min(mars_3403_1980[i,5:7])
    if(mars_3403_1980[i,29]==0){mars_3403_1980[i,32]=1}
    if(mars_3403_1980[i,30]==0){mars_3403_1980[i,32]=2}
    if(mars_3403_1980[i,31]==0){mars_3403_1980[i,32]=3}
    mars_3403_1980[i,33]=length(na.omit(ts_comp[[i]][1,]))
    mars_3403_1980[i,34]=length(na.omit(ts_comp[[i]][2,]))
    
  TS.plot.pdf(x=ts_comp[[i]],m=mars_3403_1980$mod[i],sp=mars_3403_1980$cName[i],GZ='Key Largo')
    dev.off()
    print(i)
  }


write.csv(mars_3403_1980,'MARSS_model_comparisons_3403_1980_expert_data.csv')

####Site-weighted abundances to 1980####
setwd("C:/Users/14388/Desktop/Scripps - Project 1 RVC and REEF/RVC/timeseries/Site level")

mars_3403_1980_site<- data.frame(SP=NA,conv.1=NA,conv.2=NA,conv.3=NA,tier.1.AICc=NA,tier2.AICc=NA,tier3.AICc=NA,Q1=NA,Q2=NA,Q3.rvc=NA,Q3.reef=NA,U1=NA,U2.rvc=NA,U2.reef=NA,U3.rvc=NA,U3.reef=NA,R1.rvc=NA,R1.reef=NA,R2.rvc=NA,R2.reef=NA,R3.rvc=NA,R3.reef=NA,cName=NA,obs.grp=NA,mAbund.rvc=NA,mAbund.reef=NA,sdAbund.rvc=NA,sdAbund.reef=NA,dAIC1=NA,dAIC2=NA,dAIC3=NA,mod=NA)
for(i in 1:length(rvc.green.3403.sp)){
  spp<- fish_reef[match(rvc.green.3403.sp[i],fish_reef$rvc_code),]
  ts_comp[[i]]<- t(log(rvc_3403_green[[i]]$mean_ssu_abundance[2:40])) #extract logged time-series
  m<- match(spp$scientificname,REEF.green.3403.sp)
  if(is.na(m)==T){next}
  ts_comp[[i]]<- rbind(ts_comp[[i]],t(c(rep(NA,13),REEF_3403_green[[m]]$logSiteAbund[1:26])))
  rownames(ts_comp[[i]])<- c(paste(spp$commonname,'rvc',sep=":"),paste(spp$commonname,'REEF',sep=":"))
  colnames(ts_comp[[i]])<- years
  
  Z_1=factor(c('rvc_reef','rvc_reef')) #Common state process for both time-series
  Z_2=factor(c('rvc','reef')) #Different state process for each time-series
  R_matrix[[i]]=matrix(list(0),nrow(ts_comp[[i]]),nrow(ts_comp[[i]])) #create empty Obs. error matrix
  diag(R_matrix[[i]])=c(rep("rvc",times=1),rep("reef",times=1))
  
  #tier 1 - common state process
  tier_1<- list(Z=Z_1,
                Q = 'diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 2 - Diff. U, same process error
  tier_2<- list(Z=Z_2,
                Q ='diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 3
  tier_3<- list(Z=Z_2,
                Q = 'diagonal and unequal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  fit_1<-  MARSS(ts_comp[[i]], model = tier_1,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem') 
  fit_2<-  MARSS(ts_comp[[i]], model = tier_2,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem')
  fit_3<-  MARSS(ts_comp[[i]], model = tier_3,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem')
  params.1<- MARSSparamCIs(fit_1)
  params.2<- MARSSparamCIs(fit_2)
  params.3<- MARSSparamCIs(fit_3)
  
  mars_3403_1980_site[i,1]=spp$scientificname
  mars_3403_1980_site[i,2]=ifelse(is.null(fit_1$errors)==T,1,0)
  mars_3403_1980_site[i,3]=ifelse(is.null(fit_2$errors)==T,1,0)  
  mars_3403_1980_site[i,4]=ifelse(is.null(fit_3$errors)==T,1,0)  
  mars_3403_1980_site[i,5]=fit_1$AICc
  mars_3403_1980_site[i,6]=fit_2$AICc  
  mars_3403_1980_site[i,7]=fit_3$AICc  
  mars_3403_1980_site[i,8]=params.1$parMean[5]
  mars_3403_1980_site[i,9]=params.2$parMean[5]
  mars_3403_1980_site[i,10]=params.3$parMean[5]
  mars_3403_1980_site[i,11]=params.3$parMean[6]
  mars_3403_1980_site[i,12]=params.1$parMean[4]
  mars_3403_1980_site[i,13]=params.2$parMean[3] 
  mars_3403_1980_site[i,14]=params.2$parMean[4]
  mars_3403_1980_site[i,15]=params.3$parMean[3]
  mars_3403_1980_site[i,16]=params.3$parMean[4]
  mars_3403_1980_site[i,17]=params.1$parMean[2] 
  mars_3403_1980_site[i,18]=params.1$parMean[3]
  mars_3403_1980_site[i,19]=params.2$parMean[1] 
  mars_3403_1980_site[i,20]=params.2$parMean[2]
  mars_3403_1980_site[i,21]=params.3$parMean[1] 
  mars_3403_1980_site[i,22]=params.3$parMean[2]
  mars_3403_1980_site[i,23]=spp$commonname
  mars_3403_1980_site[i,24]=spp$group_assignment
  mars_3403_1980_site[i,25]=exp(mean(na.omit(ts_comp[[i]][1,])))
  mars_3403_1980_site[i,26]=exp(mean(na.omit(ts_comp[[i]][2,])))
  mars_3403_1980_site[i,27]=sd(exp(na.omit(ts_comp[[i]][1,])))
  mars_3403_1980_site[i,28]=sd(exp(na.omit(ts_comp[[i]][2,])))
  mars_3403_1980_site[i,29]=mars_3403_1980_site[i,5]-min(mars_3403_1980_site[i,5:7])
  mars_3403_1980_site[i,30]=mars_3403_1980_site[i,6]-min(mars_3403_1980_site[i,5:7])
  mars_3403_1980_site[i,31]=mars_3403_1980_site[i,7]-min(mars_3403_1980_site[i,5:7])
  if(mars_3403_1980_site[i,29]==0){mars_3403_1980_site[i,32]=1}
  if(mars_3403_1980_site[i,30]==0){mars_3403_1980_site[i,32]=2}
  if(mars_3403_1980_site[i,31]==0){mars_3403_1980_site[i,32]=3}
  
  TS.plot.pdf(x=ts_comp[[i]],m=mars_3403_1980_site$mod[i],sp=mars_3403_1980_site$cName[i],GZ='Key Largo')
  print(i)
}


write.csv(mars_3403_1980_site,'MARSS_model_comparisons_3403_1980_site_weighted.csv')



####Model switches when truncated at 1993 vs 1980
mars_3403_80<- read.csv("MARSS_model_comparisons_3403_1980.csv")
mars_3403_93<- read.csv("MARSS_model_comparisons_3403_1993.csv")
m<- match(mars_3403_80$SP,mars_3403_93$SP)
mars_3403_93_comp<- full_join(mars_3403_80,mars_3403_93,by='SP')
mars_3403_93_comp$mod_match<- paste(mars_3403_93_comp$mod.x,mars_3403_93_comp$mod.y)



#####Simulate time-series####
fk_79_18_dams<- subset(fk_79_18,SPECIES_CD %in% c('SPA AURO','SPH BARR','MIC CHRY'))
dams_3403<- subset(fk_79_18_dams,region.id=='3403')
dams_trends_3403<- rvc_batch(dams_3403)

sp_m<- fish_reef[match(c('SPA AURO','SPH BARR','MIC CHRY'),fish_reef$rvc_code),]
m<- match(sp_m$scientificname,REEF.green.3403.sp)

ts_comp_1<- t(log(dams_trends_3403[[1]]$mean_ssu_abundance[2:40])) #extract logged time-series
ts_comp_1<- rbind(ts_comp_1,t(c(rep(NA,13),REEF_3403_green[[m[1]]]$logSiteAbund[1:26])))
rownames(ts_comp_1)<- c(paste(sp_m$commonname[1],'rvc',sep=":"),paste(sp_m$commonname[1],'REEF',sep=":"))
colnames(ts_comp_1)<- years

Z_2=factor(c('rvc','reef')) #Different state process for each time-series
R_matrix[[i]]=matrix(list(0),nrow(ts_comp_1[[i]]),nrow(ts_comp_1[[i]])) #create empty Obs. error matrix
diag(R_matrix[[i]])=c(rep("rvc",times=1),rep("reef",times=1))

#tier 3
tier_3<- list(Z=Z_2,
              Q = 'diagonal and unequal',
              R = R_matrix[[i]],
              A = 'scaling',
              x0 = 'unequal',
              tinitx=0)

fit_ts<-  MARSS(ts_comp_1, model = tier_3,control=list(maxit=50000,minit=500,conv.test.slope.tol=0.05),method='kem')

params.3<- MARSSparamCIs(fit_ts)
#R
sim_bd<- MARSSsimulate(fit_ts,tSteps=dim(ts_comp_1)[2],nsim=50,miss.loc =ts_comp_1)

sim_frame<- data.frame(sim=seq(1:50),AIC1=NA,AIC2=NA,AIC=NA,mod=NA,R1=NA,R2=NA,U1=NA,U2=NA,Q1=NA,Q2=NA,AIC1.93=NA,AIC2.93=NA,AIC3.93=NA,mod.93=NA)
for(i in 1:50){
  Z_1=factor(c('rvc_reef','rvc_reef')) #Common state process for both time-series
  Z_2=factor(c('rvc','reef')) #Different state process for each time-series
  R_matrix[[i]]=matrix(list(0),nrow(ts_comp[[i]]),nrow(ts_comp[[i]])) #create empty Obs. error matrix
  diag(R_matrix[[i]])=c(rep("rvc",times=1),rep("reef",times=1))
  
  
  #tier 1 - common state process
  tier_1<- list(Z=Z_1,
                Q = 'diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 2 - Diff. U, same process error
  tier_2<- list(Z=Z_2,
                Q ='diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 3
  tier_3<- list(Z=Z_2,
                Q = 'diagonal and unequal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  fit_1<-  MARSS(sim_bd$sim.data[,,i], model = tier_1,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem') 
  fit_2<-  MARSS(sim_bd$sim.data[,,i], model = tier_2,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem')
  fit_3<-  MARSS(sim_bd$sim.data[,,i], model = tier_3,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem')
  params.3<- MARSSparamCIs(fit_3)
  
  #Only 1993_2018
  ts_1993<- sim_bd$sim.data[,,i][1:2,14:39]
  
  fit_1_93<-  MARSS(ts_1993, model = tier_1,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem') 
  fit_2_93<-  MARSS(ts_1993, model = tier_2,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem')
  fit_3_93<-  MARSS(ts_1993, model = tier_3,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem')
  params.3<- MARSSparamCIs(fit_3)
  
  
  sim_frame[i,2]<- fit_1$AICc
  sim_frame[i,3]<- fit_2$AICc
  sim_frame[i,4]<- fit_3$AICc
  sim_frame[i,5]=which.min(sim_frame[i,2:4])
  sim_frame[i,6:11]=params.3$coef[1:6]
  sim_frame[i,12]<- fit_1_93$AICc
  sim_frame[i,13]<- fit_2_93$AICc
  sim_frame[i,14]<- fit_3_93$AICc
  sim_frame[i,15]=which.min(sim_frame[i,12:14])
}
write.csv(sim_frame,'50_simulations_Redband_Parrotfish.csv')

#Barracuda
ts_comp_2<- t(log(dams_trends_3403[[2]]$mean_ssu_abundance[2:40])) #extract logged time-series
ts_comp_2<- rbind(ts_comp_2,t(c(rep(NA,13),REEF_3403_green[[m[2]]]$logSiteAbund[1:26])))
rownames(ts_comp_2)<- c(paste(sp_m$commonname[2],'rvc',sep=":"),paste(sp_m$commonname[2],'REEF',sep=":"))
colnames(ts_comp_2)<- years

Z_2=factor(c('rvc','reef')) #Different state process for each time-series
R_matrix[[i]]=matrix(list(0),nrow(ts_comp[[i]]),nrow(ts_comp[[i]])) #create empty Obs. error matrix
diag(R_matrix[[i]])=c(rep("rvc",times=1),rep("reef",times=1))

#tier 3
tier_3<- list(Z=Z_2,
              Q = 'diagonal and unequal',
              R = R_matrix[[i]],
              A = 'scaling',
              x0 = 'unequal',
              tinitx=0)

fit_ts<-  MARSS(ts_comp_2, model = tier_3,control=list(maxit=50000,minit=500,conv.test.slope.tol=0.05),method='kem')

params.3<- MARSSparamCIs(fit_ts)
#R
sim_bd<- MARSSsimulate(fit_ts,tSteps=dim(ts_comp_2)[2],nsim=50,miss.loc =ts_comp_2)

sim_frame<- data.frame(sim=seq(1:50),AIC1=NA,AIC2=NA,AIC=NA,mod=NA,R1=NA,R2=NA,U1=NA,U2=NA,Q1=NA,Q2=NA,AIC1.93=NA,AIC2.93=NA,AIC3.93=NA,mod.93=NA)
for(i in 1:50){
  Z_1=factor(c('rvc_reef','rvc_reef')) #Common state process for both time-series
  Z_2=factor(c('rvc','reef')) #Different state process for each time-series
  R_matrix[[i]]=matrix(list(0),nrow(ts_comp[[i]]),nrow(ts_comp[[i]])) #create empty Obs. error matrix
  diag(R_matrix[[i]])=c(rep("rvc",times=1),rep("reef",times=1))
  
  
  #tier 1 - common state process
  tier_1<- list(Z=Z_1,
                Q = 'diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 2 - Diff. U, same process error
  tier_2<- list(Z=Z_2,
                Q ='diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 3
  tier_3<- list(Z=Z_2,
                Q = 'diagonal and unequal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  fit_1<-  MARSS(sim_bd$sim.data[,,i], model = tier_1,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem') 
  fit_2<-  MARSS(sim_bd$sim.data[,,i], model = tier_2,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem')
  fit_3<-  MARSS(sim_bd$sim.data[,,i], model = tier_3,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem')
  params.3<- MARSSparamCIs(fit_3)
  
  #Only 1993_2018
  ts_1993<- sim_bd$sim.data[,,i][1:2,14:39]
  
  fit_1_93<-  MARSS(ts_1993, model = tier_1,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem') 
  fit_2_93<-  MARSS(ts_1993, model = tier_2,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem')
  fit_3_93<-  MARSS(ts_1993, model = tier_3,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem')
  params.3<- MARSSparamCIs(fit_3)
  
  
  sim_frame[i,2]<- fit_1$AICc
  sim_frame[i,3]<- fit_2$AICc
  sim_frame[i,4]<- fit_3$AICc
  sim_frame[i,5]=which.min(sim_frame[i,2:4])
  sim_frame[i,6:11]=params.3$coef[1:6]
  sim_frame[i,12]<- fit_1_93$AICc
  sim_frame[i,13]<- fit_2_93$AICc
  sim_frame[i,14]<- fit_3_93$AICc
  sim_frame[i,15]=which.min(sim_frame[i,12:14])
}
write.csv(sim_frame,'50_simulations_Barracuda.csv')

#Yellowtail damsel
ts_comp_3<- t(log(dams_trends_3403[[3]]$mean_ssu_abundance[2:40])) #extract logged time-series
ts_comp_3<- rbind(ts_comp_3,t(c(rep(NA,13),REEF_3403_green[[m[3]]]$logSiteAbund[1:26])))
rownames(ts_comp_3)<- c(paste(sp_m$commonname[3],'rvc',sep=":"),paste(sp_m$commonname[3],'REEF',sep=":"))
colnames(ts_comp_3)<- years

Z_2=factor(c('rvc','reef')) #Different state process for each time-series
R_matrix[[i]]=matrix(list(0),nrow(ts_comp[[i]]),nrow(ts_comp[[i]])) #create empty Obs. error matrix
diag(R_matrix[[i]])=c(rep("rvc",times=1),rep("reef",times=1))

#tier 3
tier_3<- list(Z=Z_2,
              Q = 'diagonal and unequal',
              R = R_matrix[[i]],
              A = 'scaling',
              x0 = 'unequal',
              tinitx=0)

fit_ts<-  MARSS(ts_comp_3, model = tier_3,control=list(maxit=50000,minit=500,conv.test.slope.tol=0.05),method='kem')

params.3<- MARSSparamCIs(fit_ts)
#R
sim_bd<- MARSSsimulate(fit_ts,tSteps=dim(ts_comp_2)[2],nsim=50,miss.loc =ts_comp_2)

sim_frame<- data.frame(sim=seq(1:50),AIC1=NA,AIC2=NA,AIC=NA,mod=NA,R1=NA,R2=NA,U1=NA,U2=NA,Q1=NA,Q2=NA,AIC1.93=NA,AIC2.93=NA,AIC3.93=NA,mod.93=NA)
for(i in 1:50){
  Z_1=factor(c('rvc_reef','rvc_reef')) #Common state process for both time-series
  Z_2=factor(c('rvc','reef')) #Different state process for each time-series
  R_matrix[[i]]=matrix(list(0),nrow(ts_comp[[i]]),nrow(ts_comp[[i]])) #create empty Obs. error matrix
  diag(R_matrix[[i]])=c(rep("rvc",times=1),rep("reef",times=1))
  
  
  #tier 1 - common state process
  tier_1<- list(Z=Z_1,
                Q = 'diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 2 - Diff. U, same process error
  tier_2<- list(Z=Z_2,
                Q ='diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 3
  tier_3<- list(Z=Z_2,
                Q = 'diagonal and unequal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  fit_1<-  MARSS(sim_bd$sim.data[,,i], model = tier_1,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem') 
  fit_2<-  MARSS(sim_bd$sim.data[,,i], model = tier_2,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem')
  fit_3<-  MARSS(sim_bd$sim.data[,,i], model = tier_3,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem')
  params.3<- MARSSparamCIs(fit_3)
  
  #Only 1993_2018
  ts_1993<- sim_bd$sim.data[,,i][1:2,14:39]
  
  fit_1_93<-  MARSS(ts_1993, model = tier_1,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem') 
  fit_2_93<-  MARSS(ts_1993, model = tier_2,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem')
  fit_3_93<-  MARSS(ts_1993, model = tier_3,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem')
  params.3<- MARSSparamCIs(fit_3)
  
  
  sim_frame[i,2]<- fit_1$AICc
  sim_frame[i,3]<- fit_2$AICc
  sim_frame[i,4]<- fit_3$AICc
  sim_frame[i,5]=which.min(sim_frame[i,2:4])
  sim_frame[i,6:11]=params.3$coef[1:6]
  sim_frame[i,12]<- fit_1_93$AICc
  sim_frame[i,13]<- fit_2_93$AICc
  sim_frame[i,14]<- fit_3_93$AICc
  sim_frame[i,15]=which.min(sim_frame[i,12:14])
}
write.csv(sim_frame,'50_simulations_Yellowtail_Damsel.csv')


#Mod 1 example
fk_79_18_chroms<- subset(fk_79_18,SPECIES_CD %in% c('CHR MULT'))
chroms_3403<- subset(fk_79_18_chroms,region.id=='3403')
chroms_trends_3403<- rvc_batch(chroms_3403)

sp_m<- fish_reef[match(c('CHR MULT'),fish_reef$rvc_code),]
m<- match(sp_m$scientificname,REEF.green.3403.sp)

ts_comp_1<- t(log(chroms_trends_3403[[1]]$mean_ssu_abundance[2:40])) #extract logged time-series
ts_comp_1<- rbind(ts_comp_1,t(c(rep(NA,13),REEF_3403_green[[m[1]]]$logSiteAbund[1:26])))
rownames(ts_comp_1)<- c(paste(sp_m$commonname[1],'rvc',sep=":"),paste(sp_m$commonname[1],'REEF',sep=":"))
colnames(ts_comp_1)<- years

Z_2=factor(c('rvc','reef')) #Different state process for each time-series
R_matrix[[i]]=matrix(list(0),nrow(ts_comp[[i]]),nrow(ts_comp[[i]])) #create empty Obs. error matrix
diag(R_matrix[[i]])=c(rep("rvc",times=1),rep("reef",times=1))

#tier 1
tier_1<- list(Z=Z_1,
              Q = 'diagonal and equal',
              R = R_matrix[[i]],
              A = 'scaling',
              x0 = 'unequal',
              tinitx=0)

fit_ts<-  MARSS(ts_comp_1, model = tier_1,control=list(maxit=50000,minit=500,conv.test.slope.tol=0.05),method='kem')

params.3<- MARSSparamCIs(fit_ts)

#
sim_bd<- MARSSsimulate(fit_ts,tSteps=dim(ts_comp_1)[2],nsim=50,miss.loc=ts_comp_1)

sim_frame_c<- data.frame(sim=seq(1:50),AIC1=NA,AIC2=NA,AIC=NA,mod=NA,R1=NA,R2=NA,U1=NA,U2=NA,Q1=NA,Q2=NA,AIC1.93=NA,AIC2.93=NA,AIC3.93=NA,mod.93=NA)
for(i in 1:50){
  Z_1=factor(c('rvc_reef','rvc_reef')) #Common state process for both time-series
  Z_2=factor(c('rvc','reef')) #Different state process for each time-series
  R_matrix[[i]]=matrix(list(0),nrow(ts_comp[[i]]),nrow(ts_comp[[i]])) #create empty Obs. error matrix
  diag(R_matrix[[i]])=c(rep("rvc",times=1),rep("reef",times=1))
  
  
  #tier 1 - common state process
  tier_1<- list(Z=Z_1,
                Q = 'diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 2 - Diff. U, same process error
  tier_2<- list(Z=Z_2,
                Q ='diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 3
  tier_3<- list(Z=Z_2,
                Q = 'diagonal and unequal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  fit_1<-  MARSS(sim_bd$sim.data[,,i], model = tier_1,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem') 
  fit_2<-  MARSS(sim_bd$sim.data[,,i], model = tier_2,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem')
  fit_3<-  MARSS(sim_bd$sim.data[,,i], model = tier_3,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem')
  params.3<- MARSSparamCIs(fit_3)
  
  #Only 1993_2018
  ts_1993<- sim_bd$sim.data[,,i][1:2,14:39]
  
  fit_1_93<-  MARSS(ts_1993, model = tier_1,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem') 
  fit_2_93<-  MARSS(ts_1993, model = tier_2,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem')
  fit_3_93<-  MARSS(ts_1993, model = tier_3,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem')
  params.1<- MARSSparamCIs(fit_1)
  
  
  sim_frame_c[i,2]<- fit_1$AICc
  sim_frame_c[i,3]<- fit_2$AICc
  sim_frame_c[i,4]<- fit_3$AICc
  sim_frame_c[i,5]=which.min(sim_frame_c[i,2:4])
  sim_frame_c[i,6:11]=params.1$coef[1:6]
  sim_frame_c[i,12]<- fit_1_93$AICc
  sim_frame_c[i,13]<- fit_2_93$AICc
  sim_frame_c[i,14]<- fit_3_93$AICc
  sim_frame_c[i,15]=which.min(sim_frame_c[i,12:14])
}
write.csv(sim_frame_c,'50_simulations_Brown_Chromis.csv')

sim_frame_bc<- read.csv('50_simulations_Brown_Chromis.csv')


####Large simulation of 71 series####
rvc_3403<- subset(fk_79_18_trim,region.id=='3403')
rvc_trends_3403<- rvc_batch(rvc_3403)
rvc_3403_green<-list.filter(rvc_trends_3403,length(na.omit(mean_ssu_abundance))>=28) #71 species


REEF.3403.comb.green<- do.call(rbind, lapply(REEF_3403_green, data.frame, stringsAsFactors=FALSE))
REEF.green.3403.sp<- unique(REEF.3403.comb.green$sciName)

rvc.3403.comb<- do.call(rbind, lapply(rvc_3403_green, data.frame, stringsAsFactors=FALSE))
rvc.green.3403.sp<- unique(rvc.3403.comb$Species)

ts_comp<- list()
years<- seq(1980,2018,by=1)
n<- NA
R_matrix<- list()

Params_list<- list()
mod_df<- data.frame(count_1=NA,count_2=NA,count_3=NA)
for(z in 1:500){
  sim_params<- data.frame(SP1=NA,SP2=NA,conv.1=NA,conv.2=NA,conv.3=NA,tier.1.AICc=NA,tier2.AICc=NA,tier3.AICc=NA,dAIC1=NA,dAIC2=NA,dAIC3=NA,mod=NA,wAIC1=NA,wAIC2=NA,wAIC3=NA,sim.run=NA)
  
  for(i in 1:length(rvc.green.3403.sp)){
    ts_comp<- t(log(rvc_3403_green[[i]]$mean_ssu_abundance[2:40])) #extract logged time-series
    sci_n<- fish_reef$scientificname[match(unique(rvc_3403_green[[i]]$Species),fish_reef$rvc_code)]
    spp<- subset(REEF.green.3403.sp,REEF.green.3403.sp != sci_n) #Remove target species from species list
    spp_r<- sample(spp,1) #randomly sample a new species
    m<- match(spp_r,REEF.green.3403.sp) 
    ts_comp<- rbind(ts_comp,t(c(rep(NA,13),REEF_3403_green[[m]]$logAbund[1:26])))
    rownames(ts_comp)<- c(paste(fish_reef$scientificname[match(unique(rvc_3403_green[[i]]$Species),fish_reef$rvc_code)],'rvc',sep=":"),paste(spp_r,'REEF',sep=":"))
    colnames(ts_comp)<- years
    
    Z_1=factor(c('rvc_reef','rvc_reef')) #Common state process for both time-series
    Z_2=factor(c('rvc','reef')) #Different state process for each time-series
    R_matrix=matrix(list(0),nrow(ts_comp),nrow(ts_comp)) #create empty Obs. error matrix
    diag(R_matrix)=c(rep("rvc",times=1),rep("reef",times=1))
    
    #tier 1 - common state process
    tier_1<- list(Z=Z_1,
                  Q = 'diagonal and equal',
                  R = R_matrix,
                  A = 'scaling',
                  x0 = 'unequal',
                  tinitx=0)
    
    #tier 2 - Diff. U, same process error
    tier_2<- list(Z=Z_2,
                  Q ='diagonal and equal',
                  R = R_matrix,
                  A = 'scaling',
                  x0 = 'unequal',
                  tinitx=0)
    
    #tier 3
    tier_3<- list(Z=Z_2,
                  Q = 'diagonal and unequal',
                  R = R_matrix,
                  A = 'scaling',
                  x0 = 'unequal',
                  tinitx=0)
    
    fit_1<-  MARSS(ts_comp, model = tier_1,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem') 
    fit_2<-  MARSS(ts_comp, model = tier_2,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem')
    fit_3<-  MARSS(ts_comp, model = tier_3,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem')
    params.1<- MARSSparamCIs(fit_1)
    params.2<- MARSSparamCIs(fit_2)
    params.3<- MARSSparamCIs(fit_3)
    
    sim_params[i,1]=fish_reef$commonname[match(unique(rvc_3403_green[[i]]$Species),fish_reef$rvc_code)]
    sim_params[i,2]=fish_reef$commonname[match(spp_r,fish_reef$scientificname)]
    sim_params[i,3]=ifelse(is.null(fit_1$errors)==T,1,0)
    sim_params[i,4]=ifelse(is.null(fit_2$errors)==T,1,0)  
    sim_params[i,5]=ifelse(is.null(fit_3$errors)==T,1,0)  
    sim_params[i,6]=fit_1$AICc
    sim_params[i,7]=fit_2$AICc  
    sim_params[i,8]=fit_3$AICc
    sim_params[i,9]=sim_params[i,6]-min(sim_params[i,6:8])
    sim_params[i,10]=sim_params[i,7]-min(sim_params[i,6:8])
    sim_params[i,11]=sim_params[i,8]-min(sim_params[i,6:8])
    if(sim_params[i,9]==0){sim_params[i,12]=1}
    if(sim_params[i,10]==0){sim_params[i,12]=2}
    if(sim_params[i,11]==0){sim_params[i,12]=3}
    sim_params[i,13]=exp(-0.5*sim_params[i,9])/sum(c(exp(-0.5*sim_params[i,9]),exp(-0.5*sim_params[i,10]),exp(-0.5*sim_params[i,11])))
    sim_params[i,14]=exp(-0.5*sim_params[i,10])/sum(c(exp(-0.5*sim_params[i,9]),exp(-0.5*sim_params[i,10]),exp(-0.5*sim_params[i,11])))
    sim_params[i,15]=exp(-0.5*sim_params[i,11])/sum(c(exp(-0.5*sim_params[i,9]),exp(-0.5*sim_params[i,10]),exp(-0.5*sim_params[i,11])))
    sim_params[i,16]=z
    print(i)
  }
  Params_list[[z]]<- list(sim_params)
  mod_df[z,1]=summary(as.factor(sim_params$mod))[1]
  mod_df[z,2]=summary(as.factor(sim_params$mod))[2]
  mod_df[z,3]=summary(as.factor(sim_params$mod))[3]
}

sims_dataframe<- do.call(rbind, lapply(Params_list, data.frame, stringsAsFactors=FALSE))
write.csv(sims_dataframe,'Timeseries_randomized_simulations_Nov.csv')

mod1<- subset(sims_dataframe,mod==1)

####Simulations 1993 data####
#No trend
Params_list<- list()
mod_df<- data.frame(count_1=NA,count_2=NA,count_3=NA)
for(z in 1:500){
  sim_params<- data.frame(SP1=NA,SP2=NA,conv.1=NA,conv.2=NA,tier.1.AICc=NA,tier2.AICc=NA,dAIC1=NA,dAIC2=NA,mod=NA,wAIC1=NA,wAIC2=NA,sim.run=NA)
  
mars_3403_nt<- data.frame(SP=NA,conv.1=NA,conv.2=NA,tier.1.AICc=NA,tier2.AICc=NA,Q1=NA,Q2.rvc=NA,Q2.reef=NA,R1.rvc=NA,R1.reef=NA,R2.rvc=NA,R2.reef=NA,cName=NA,obs.grp=NA,mAbund.rvc=NA,mAbund.reef=NA,sdAbund.rvc=NA,sdAbund.reef=NA,dAIC1=NA,dAIC2=NA,mod=NA,years.rvc=NA,years.reef=NA)
for(i in 1:length(rvc.green.3403.sp)){
  ts_comp<- t(log(rvc_3403_green[[i]]$mean_ssu_abundance[1:26])) #extract logged time-series
  sci_n<- fish_reef$scientificname[match(unique(rvc_3403_green[[i]]$Species),fish_reef$rvc_code)]
  spp<- subset(REEF.green.3403.sp,REEF.green.3403.sp != sci_n) #Remove target species from species list
  spp_r<- sample(spp,1) #randomly sample a new species
  m<- match(spp_r,REEF.green.3403.sp) 
  ts_comp<- rbind(ts_comp,t(REEF_3403_green[[m]]$logAbund[1:26]))
  rownames(ts_comp)<- c(paste(fish_reef$scientificname[match(unique(rvc_3403_green[[i]]$Species),fish_reef$rvc_code)],'rvc',sep=":"),paste(spp_r,'REEF',sep=":"))
  colnames(ts_comp)<- years
  
  Z_1=factor(c('rvc_reef','rvc_reef')) #Common state process for both time-series
  Z_2=factor(c('rvc','reef')) #Different state process for each time-series
  R_matrix=matrix(list(0),nrow(ts_comp),nrow(ts_comp)) #create empty Obs. error matrix
  diag(R_matrix)=c(rep("rvc",times=1),rep("reef",times=1))
  
  #tier 1 - common state process
  tier_1<- list(Z=Z_1,
                Q = 'diagonal and equal',
                U = 'zero',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 2 - Diff. U, same process error
  tier_2<- list(Z=Z_2,
                Q ='diagonal and unequal',
                R = R_matrix[[i]],
                A = 'scaling',
                U = 'zero',
                x0 = 'unequal',
                tinitx=0)
  
  
  fit_1<-  MARSS(ts_comp[[i]], model = tier_1,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  fit_2<-  MARSS(ts_comp[[i]], model = tier_2,control=list(maxit=30000,minit=500,conv.test.slope.tol = 0.1),method='kem') 
  params.1<- MARSSparamCIs(fit_1)
  params.2<- MARSSparamCIs(fit_2)
  
  mars_3403_nt[i,1]=spp$scientificname
  mars_3403_nt[i,2]=ifelse(is.null(fit_1$errors)==T,1,0)
  mars_3403_nt[i,3]=ifelse(is.null(fit_2$errors)==T,1,0)  
  mars_3403_nt[i,4]=fit_1$AICc
  mars_3403_nt[i,5]=fit_2$AICc  
  mars_3403_nt[i,6]=params.1$parMean[4]
  mars_3403_nt[i,7]=params.2$parMean[3]
  mars_3403_nt[i,8]=params.2$parMean[4]
  mars_3403_nt[i,9]=params.1$parMean[2]
  mars_3403_nt[i,10]=params.1$parMean[3]
  mars_3403_nt[i,11]=params.2$parMean[1]
  mars_3403_nt[i,12]=params.2$parMean[2]
  mars_3403_nt[i,13]=spp$commonname
  mars_3403_nt[i,14]=spp$group_assignment
  mars_3403_nt[i,15]=exp(mean(na.omit(ts_comp[[i]][1,])))
  mars_3403_nt[i,16]=exp(mean(na.omit(ts_comp[[i]][2,])))
  mars_3403_nt[i,17]=sd(exp(na.omit(ts_comp[[i]][1,])))
  mars_3403_nt[i,18]=sd(exp(na.omit(ts_comp[[i]][2,])))
  mars_3403_nt[i,19]=mars_3403_nt[i,4]-min(mars_3403_nt[i,4:5])
  mars_3403_nt[i,20]=mars_3403_nt[i,5]-min(mars_3403_nt[i,4:5])
  if(mars_3403_nt[i,19]==0){mars_3403_nt[i,21]=1}
  if(mars_3403_nt[i,20]==0){mars_3403_nt[i,21]=2}
  mars_3403_nt[i,22]=length(na.omit(ts_comp[[i]][1,]))
  mars_3403_nt[i,23]=length(na.omit(ts_comp[[i]][2,]))
  
  #  TS.plot.pdf(x=ts_comp[[i]],m=mars_3403$mod[i],sp=mars_3403$cName[i],GZ='Key Largo')
  #  dev.off()
  print(i)
}

Params_list[[z]]<- list(sim_params)
mod_df[z,1]=summary(as.factor(sim_params$mod))[1]
mod_df[z,2]=summary(as.factor(sim_params$mod))[2]
mod_df[z,3]=summary(as.factor(sim_params$mod))[3]
}

sims_dataframe<- do.call(rbind, lapply(Params_list, data.frame, stringsAsFactors=FALSE))
write.csv(sims_dataframe,'Timeseries_randomized_simulations_1993_Nov.csv')



#Trim down RVC dataset to these species
fk_93_18_trim<- subset(fk_93_18,SPECIES_CD %in% fish_reef$rvc_code)

#Separate out RVC data into the REEF sub-regions
rvc_3403<- subset(fk_93_18_trim,region.id=='3403')

rvc_trends_3403<- rvc_batch(rvc_3403)
rvc_3403_green<-list.filter(rvc_trends_3403,length(na.omit(mean_ssu_abundance))>=18) #101 species


REEF.3403.comb.green<- do.call(rbind, lapply(REEF_3403_green, data.frame, stringsAsFactors=FALSE))
REEF.green.3403.sp<- unique(REEF.3403.comb.green$sciName)

rvc.3403.comb<- do.call(rbind, lapply(rvc_3403_green, data.frame, stringsAsFactors=FALSE))
rvc.green.3403.sp<- unique(rvc.3403.comb$Species)

ts_comp<- list()
years<- seq(1993,2018,by=1)
n<- NA
R_matrix<- list()

Params_list<- list()
mod_df<- data.frame(count_1=NA,count_2=NA)
for(z in 1:500){
  sim_params<- data.frame(SP1=NA,SP2=NA,conv.1=NA,conv.2=NA,tier.1.AICc=NA,tier2.AICc=NA,dAIC1=NA,dAIC2=NA,mod=NA,wAIC1=NA,wAIC2=NA,sim.run=NA)
  
  for(i in 1:length(rvc.green.3403.sp)){
    ts_comp<- t(log(rvc_3403_green[[i]]$mean_ssu_abundance[1:26])) #extract logged time-series
    sci_n<- fish_reef$scientificname[match(unique(rvc_3403_green[[i]]$Species),fish_reef$rvc_code)]
    spp<- subset(REEF.green.3403.sp,REEF.green.3403.sp != sci_n) #Remove target species from species list
    spp_r<- sample(spp,1) #randomly sample a new species
    m<- match(spp_r,REEF.green.3403.sp) 
    ts_comp<- rbind(ts_comp,t(REEF_3403_green[[m]]$logAbund[1:26]))
    rownames(ts_comp)<- c(paste(fish_reef$scientificname[match(unique(rvc_3403_green[[i]]$Species),fish_reef$rvc_code)],'rvc',sep=":"),paste(spp_r,'REEF',sep=":"))
    colnames(ts_comp)<- years
    
    Z_1=factor(c('rvc_reef','rvc_reef')) #Common state process for both time-series
    Z_2=factor(c('rvc','reef')) #Different state process for each time-series
    R_matrix=matrix(list(0),nrow(ts_comp),nrow(ts_comp)) #create empty Obs. error matrix
    diag(R_matrix)=c(rep("rvc",times=1),rep("reef",times=1))
    
    #tier 1 - common state process
    tier_1<- list(Z=Z_1,
                  Q = 'diagonal and equal',
                  U = 'zero',
                  R = R_matrix,
                  A = 'scaling',
                  x0 = 'unequal',
                  tinitx=0)
    
    #tier 2 - Diff. U, same process error
    tier_2<- list(Z=Z_2,
                  Q ='diagonal and unequal',
                  R = R_matrix,
                  A = 'scaling',
                  U = 'zero',
                  x0 = 'unequal',
                  tinitx=0)
    
    fit_1<-  MARSS(ts_comp, model = tier_1,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem') 
    fit_2<-  MARSS(ts_comp, model = tier_2,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem')
    params.1<- MARSSparamCIs(fit_1)
    params.2<- MARSSparamCIs(fit_2)
    
    sim_params[i,1]=fish_reef$commonname[match(unique(rvc_3403_green[[i]]$Species),fish_reef$rvc_code)]
    sim_params[i,2]=fish_reef$commonname[match(spp_r,fish_reef$scientificname)]
    sim_params[i,3]=ifelse(is.null(fit_1$errors)==T,1,0)
    sim_params[i,4]=ifelse(is.null(fit_2$errors)==T,1,0)  
    sim_params[i,5]=fit_1$AICc
    sim_params[i,6]=fit_2$AICc  
    sim_params[i,7]=sim_params[i,5]-min(sim_params[i,5:6])
    sim_params[i,8]=sim_params[i,6]-min(sim_params[i,5:6])
    if(sim_params[i,8]==0){sim_params[i,9]=1}
    if(sim_params[i,8]==0){sim_params[i,9]=2}
    sim_params[i,10]=exp(-0.5*sim_params[i,7])/sum(c(exp(-0.5*sim_params[i,7]),exp(-0.5*sim_params[i,8])))
    sim_params[i,11]=exp(-0.5*sim_params[i,8])/sum(c(exp(-0.5*sim_params[i,7]),exp(-0.5*sim_params[i,8])))
    sim_params[i,12]=z
    print(i)
  }
  
  Params_list[[z]]<- list(sim_params)
  mod_df[z,1]=summary(as.factor(sim_params$mod))[1]
  mod_df[z,2]=summary(as.factor(sim_params$mod))[2]
}

sims_dataframe<- do.call(rbind, lapply(Params_list, data.frame, stringsAsFactors=FALSE))
write.csv(sims_dataframe,'Timeseries_randomized_simulations_1993_Nov.csv')

mod1<- subset(sims_dataframe,mod==1)



####Different abundance metric comparisons####
#Region 3403
REEF.3403.comb.green<- do.call(rbind, lapply(REEF_3403_green, data.frame, stringsAsFactors=FALSE))
REEF.green.3403.sp<- unique(REEF.3403.comb.green$sciName)

REEF.3403.exp.comb.green<- do.call(rbind, lapply(REEF_3403_exp_green, data.frame, stringsAsFactors=FALSE))
REEF.green.exp.3403.sp<- unique(REEF.3403.exp.comb.green$sciName)


rvc.3403.comb<- do.call(rbind, lapply(rvc_3403_green, data.frame, stringsAsFactors=FALSE))
rvc.green.3403.sp<- unique(rvc.3403.comb$Species)

setwd("C:/Users/14388/Desktop/Scripps - Project 1 RVC and REEF/RVC/timeseries_comp")
load("C:/Users/14388/Desktop/Scripps - Project 1 RVC and REEF/RVC/timeseries_comp/ts_plots.RData")
for(i in 1:50){
  par(xpd=T)
  pdf(file=paste(paste(fish_reef$commonname[match(unique(REEF_3403_green[[i]]$sciName),fish_reef$scientificname)],'Key Largo',sep="-"),'.pdf',sep=''),width=8,height=6)
  plot(c(scale(logAbund))~year,data=REEF_3403_green[[i]],ylab='Scaled Abundance (std devs.)',xlab='Year',bty='l',type='n',main=paste(fish_reef$commonname[match(unique(REEF_3403_green[[i]]$sciName),fish_reef$scientificname)],'Key Largo',sep="-"),ylim=c(-3.5,3.5))
  lines(c(scale(logAbund))~year,data=REEF_3403_green[[i]],lwd=2,col='navy')
  points(c(scale(logAbund))~year,data=REEF_3403_green[[i]],pch=21,col='white',bg='navy',cex=1.5)
  lines(c(scale(logSiteAbund))~year,data=REEF_3403_green[[i]],lwd=2,col='darkcyan')
  points(c(scale(logSiteAbund))~year,data=REEF_3403_green[[i]],pch=21,col='white',bg='darkcyan',cex=1.5)
  
  m<- match(unique(REEF_3403_green[[i]]$sciName),REEF.green.exp.3403.sp)
  if(is.na(m)==F){
    lines(c(scale(logAbund))~year,data=REEF_3403_exp_green[[m]],lwd=2,col='coral4')
    points(c(scale(logAbund))~year,data=REEF_3403_exp_green[[m]],pch=21,col='white',bg='coral4',cex=1.5)
    lines(c(scale(logSiteAbund))~year,data=REEF_3403_exp_green[[m]],lwd=2,col='darkgoldenrod4')
    points(c(scale(logSiteAbund))~year,data=REEF_3403_exp_green[[m]],pch=21,col='white',bg='darkgoldenrod4',cex=1.5)
    
  }
  
  m2<- match(unique(REEF_3403_green[[i]]$sciName),fish_reef$scientificname)
  rvc_code<- fish_reef$rvc_code[m]
  m3<- match(rvc_code,rvc.green.3403.sp)
  if(is.na(m3)==F){
    lines(scale(log(mean_ssu_abundance))~YEAR,data=rvc_3403_green[[m3]],lwd=2,col='darkgray')
    points(scale(log(mean_ssu_abundance))~YEAR,data=rvc_3403_green[[m3]],pch=21,col='white',bg='darkgray',cex=1.5)
  }
  legend(x=2012,y=4,c('Nov+Exp survey abund.','Nov+Exp site abund.','Exp abund.','Exp site abund.','RVC'),text.col=c('navy','darkcyan','coral4','darkgoldenrod4','darkgray'),bty='n')
  dev.off(paste(paste(fish_reef$commonname[match(unique(REEF_3403_green[[i]]$sciName),fish_reef$scientificname)],'Key Largo',sep="-"),'.pdf',sep=''))
}

for(i in 51:100){
  par(xpd=T)
  pdf(file=paste(paste(fish_reef$commonname[match(unique(REEF_3403_green[[i]]$sciName),fish_reef$scientificname)],'Key Largo',sep="-"),'.pdf',sep=''),width=8,height=6)
  plot(c(scale(logAbund))~year,data=REEF_3403_green[[i]],ylab='Scaled Abundance (std devs.)',xlab='Year',bty='l',type='n',main=paste(fish_reef$commonname[match(unique(REEF_3403_green[[i]]$sciName),fish_reef$scientificname)],'Key Largo',sep="-"),ylim=c(-3.5,3.5))
  lines(c(scale(logAbund))~year,data=REEF_3403_green[[i]],lwd=2,col='navy')
  points(c(scale(logAbund))~year,data=REEF_3403_green[[i]],pch=21,col='white',bg='navy',cex=1.5)
  lines(c(scale(logSiteAbund))~year,data=REEF_3403_green[[i]],lwd=2,col='darkcyan')
  points(c(scale(logSiteAbund))~year,data=REEF_3403_green[[i]],pch=21,col='white',bg='darkcyan',cex=1.5)
  
  m<- match(unique(REEF_3403_green[[i]]$sciName),REEF.green.exp.3403.sp)
  if(is.na(m)==F){
    lines(c(scale(logAbund))~year,data=REEF_3403_exp_green[[m]],lwd=2,col='coral4')
    points(c(scale(logAbund))~year,data=REEF_3403_exp_green[[m]],pch=21,col='white',bg='coral4',cex=1.5)
    lines(c(scale(logSiteAbund))~year,data=REEF_3403_exp_green[[m]],lwd=2,col='darkgoldenrod4')
    points(c(scale(logSiteAbund))~year,data=REEF_3403_exp_green[[m]],pch=21,col='white',bg='darkgoldenrod4',cex=1.5)
    
  }
  
  m2<- match(unique(REEF_3403_green[[i]]$sciName),fish_reef$scientificname)
  rvc_code<- fish_reef$rvc_code[m]
  m3<- match(rvc_code,rvc.green.3403.sp)
  if(is.na(m3)==F){
    lines(scale(log(mean_ssu_abundance))~YEAR,data=rvc_3403_green[[m3]],lwd=2,col='darkgray')
    points(scale(log(mean_ssu_abundance))~YEAR,data=rvc_3403_green[[m3]],pch=21,col='white',bg='darkgray',cex=1.5)
  }
  legend(x=2012,y=4,c('Nov+Exp survey abund.','Nov+Exp site abund.','Exp abund.','Exp site abund.','RVC'),text.col=c('navy','darkcyan','coral4','darkgoldenrod4','darkgray'),bty='n')
  dev.off(paste(paste(fish_reef$commonname[match(unique(REEF_3403_green[[i]]$sciName),fish_reef$scientificname)],'Key Largo',sep="-"),'.pdf',sep=''))
}

for(i in 101:150){
  par(xpd=T)
  pdf(file=paste(paste(fish_reef$commonname[match(unique(REEF_3403_green[[i]]$sciName),fish_reef$scientificname)],'Key Largo',sep="-"),'.pdf',sep=''),width=8,height=6)
  plot(c(scale(logAbund))~year,data=REEF_3403_green[[i]],ylab='Scaled Abundance (std devs.)',xlab='Year',bty='l',type='n',main=paste(fish_reef$commonname[match(unique(REEF_3403_green[[i]]$sciName),fish_reef$scientificname)],'Key Largo',sep="-"),ylim=c(-3.5,3.5))
  lines(c(scale(logAbund))~year,data=REEF_3403_green[[i]],lwd=2,col='navy')
  points(c(scale(logAbund))~year,data=REEF_3403_green[[i]],pch=21,col='white',bg='navy',cex=1.5)
  lines(c(scale(logSiteAbund))~year,data=REEF_3403_green[[i]],lwd=2,col='darkcyan')
  points(c(scale(logSiteAbund))~year,data=REEF_3403_green[[i]],pch=21,col='white',bg='darkcyan',cex=1.5)
  
  m<- match(unique(REEF_3403_green[[i]]$sciName),REEF.green.exp.3403.sp)
  if(is.na(m)==F){
    lines(c(scale(logAbund))~year,data=REEF_3403_exp_green[[m]],lwd=2,col='coral4')
    points(c(scale(logAbund))~year,data=REEF_3403_exp_green[[m]],pch=21,col='white',bg='coral4',cex=1.5)
    lines(c(scale(logSiteAbund))~year,data=REEF_3403_exp_green[[m]],lwd=2,col='darkgoldenrod4')
    points(c(scale(logSiteAbund))~year,data=REEF_3403_exp_green[[m]],pch=21,col='white',bg='darkgoldenrod4',cex=1.5)
    
  }
  
  m2<- match(unique(REEF_3403_green[[i]]$sciName),fish_reef$scientificname)
  rvc_code<- fish_reef$rvc_code[m]
  m3<- match(rvc_code,rvc.green.3403.sp)
  if(is.na(m3)==F){
    lines(scale(log(mean_ssu_abundance))~YEAR,data=rvc_3403_green[[m3]],lwd=2,col='darkgray')
    points(scale(log(mean_ssu_abundance))~YEAR,data=rvc_3403_green[[m3]],pch=21,col='white',bg='darkgray',cex=1.5)
  }
  legend(x=2012,y=4,c('Nov+Exp survey abund.','Nov+Exp site abund.','Exp abund.','Exp site abund.','RVC'),text.col=c('navy','darkcyan','coral4','darkgoldenrod4','darkgray'),bty='n')
  dev.off(paste(paste(fish_reef$commonname[match(unique(REEF_3403_green[[i]]$sciName),fish_reef$scientificname)],'Key Largo',sep="-"),'.pdf',sep=''))
}

for(i in 151:length(REEF.green.3403.sp)){
  par(xpd=T)
  pdf(file=paste(paste(fish_reef$commonname[match(unique(REEF_3403_green[[i]]$sciName),fish_reef$scientificname)],'Key Largo',sep="-"),'.pdf',sep=''),width=8,height=6)
  plot(c(scale(logAbund))~year,data=REEF_3403_green[[i]],ylab='Scaled Abundance (std devs.)',xlab='Year',bty='l',type='n',main=paste(fish_reef$commonname[match(unique(REEF_3403_green[[i]]$sciName),fish_reef$scientificname)],'Key Largo',sep="-"),ylim=c(-3.5,3.5))
  lines(c(scale(logAbund))~year,data=REEF_3403_green[[i]],lwd=2,col='navy')
  points(c(scale(logAbund))~year,data=REEF_3403_green[[i]],pch=21,col='white',bg='navy',cex=1.5)
  lines(c(scale(logSiteAbund))~year,data=REEF_3403_green[[i]],lwd=2,col='darkcyan')
  points(c(scale(logSiteAbund))~year,data=REEF_3403_green[[i]],pch=21,col='white',bg='darkcyan',cex=1.5)
  
  m<- match(unique(REEF_3403_green[[i]]$sciName),REEF.green.exp.3403.sp)
  if(is.na(m)==F){
    lines(c(scale(logAbund))~year,data=REEF_3403_exp_green[[m]],lwd=2,col='coral4')
    points(c(scale(logAbund))~year,data=REEF_3403_exp_green[[m]],pch=21,col='white',bg='coral4',cex=1.5)
    lines(c(scale(logSiteAbund))~year,data=REEF_3403_exp_green[[m]],lwd=2,col='darkgoldenrod4')
    points(c(scale(logSiteAbund))~year,data=REEF_3403_exp_green[[m]],pch=21,col='white',bg='darkgoldenrod4',cex=1.5)
    
  }
  
  m2<- match(unique(REEF_3403_green[[i]]$sciName),fish_reef$scientificname)
  rvc_code<- fish_reef$rvc_code[m]
  m3<- match(rvc_code,rvc.green.3403.sp)
  if(is.na(m3)==F){
    lines(scale(log(mean_ssu_abundance))~YEAR,data=rvc_3403_green[[m3]],lwd=2,col='darkgray')
    points(scale(log(mean_ssu_abundance))~YEAR,data=rvc_3403_green[[m3]],pch=21,col='white',bg='darkgray',cex=1.5)
  }
  legend(x=2012,y=4,c('Nov+Exp survey abund.','Nov+Exp site abund.','Exp abund.','Exp site abund.','RVC'),text.col=c('navy','darkcyan','coral4','darkgoldenrod4','darkgray'),bty='n')
  dev.off(paste(paste(fish_reef$commonname[match(unique(REEF_3403_green[[i]]$sciName),fish_reef$scientificname)],'Key Largo',sep="-"),'.pdf',sep=''))
}

####Occurrences####
setwd("C:/Users/14388/Desktop/Scripps - Project 1 RVC and REEF/RVC/occupancy timeseries")
load("C:/Users/14388/Desktop/Scripps - Project 1 RVC and REEF/RVC/timeseries_comp/ts_plots.RData")

for(i in 1:50){
  par(xpd=T)
  pdf(file=paste(paste(fish_reef$commonname[match(unique(REEF_3403_green[[i]]$sciName),fish_reef$scientificname)],'Key Largo',sep="-"),'.pdf',sep=''),width=8,height=6)
  plot(c(n.occ/tot.surveys)~year,data=REEF_3403_green[[i]],ylab='P(Occurring)',xlab='Year',bty='l',type='n',main=paste(fish_reef$commonname[match(unique(REEF_3403_green[[i]]$sciName),fish_reef$scientificname)],'Key Largo',sep="-"),ylim=c(0,1))
  lines(c(n.occ/tot.surveys)~year,data=REEF_3403_green[[i]],lwd=2,col='navy')
  points(c(n.occ/tot.surveys)~year,data=REEF_3403_green[[i]],pch=21,col='white',bg='navy',cex=1.5)
  lines(c(occ.sites/sites.surveyed)~year,data=REEF_3403_green[[i]],lwd=2,col='darkcyan')
  points(c(occ.sites/sites.surveyed)~year,data=REEF_3403_green[[i]],pch=21,col='white',bg='darkcyan',cex=1.5)
  
  m<- match(unique(REEF_3403_green[[i]]$sciName),REEF.green.exp.3403.sp)
  if(is.na(m)==F){
    lines(c(n.occ/tot.surveys)~year,data=REEF_3403_exp_green[[m]],lwd=2,col='coral4')
    points(c(n.occ/tot.surveys)~year,data=REEF_3403_exp_green[[m]],pch=21,col='white',bg='coral4',cex=1.5)
    lines(c(occ.sites/sites.surveyed)~year,data=REEF_3403_exp_green[[m]],lwd=2,col='darkgoldenrod4')
    points(c(occ.sites/sites.surveyed)~year,data=REEF_3403_exp_green[[m]],pch=21,col='white',bg='darkgoldenrod4',cex=1.5)
    
  }
  
  m2<- match(unique(REEF_3403_green[[i]]$sciName),fish_reef$scientificname)
  rvc_code<- fish_reef$rvc_code[m]
  m3<- match(rvc_code,rvc.green.3403.sp)
  if(is.na(m3)==F){
    lines(c(ssu.occ/nm)~YEAR,data=rvc_3403_green[[m3]],lwd=2,col='darkgray')
    points(c(ssu.occ/nm)~YEAR,data=rvc_3403_green[[m3]],pch=21,col='white',bg='darkgray',cex=1.5)
    lines(c(psu.occ/n)~YEAR,data=rvc_3403_green[[m3]],lwd=2,col='antiquewhite3')
    points(c(psu.occ/n)~YEAR,data=rvc_3403_green[[m3]],pch=21,col='white',bg='antiquewhite3',cex=1.5)
    
  }
  legend(x=2012,y=1.05,c('Nov+Exp surv freq.','Nov+Exp site freq.','Exp surv freq.','Exp site freq.','RVC surv freq.','RVC site freq.'),text.col=c('navy','darkcyan','coral4','darkgoldenrod4','darkgray','antiquewhite3'),bty='n',y.intersp = 0.7)
  dev.off(paste(paste(fish_reef$commonname[match(unique(REEF_3403_green[[i]]$sciName),fish_reef$scientificname)],'Key Largo',sep="-"),'.pdf',sep=''))
}

for(i in 51:100){
  par(xpd=T)
  pdf(file=paste(paste(fish_reef$commonname[match(unique(REEF_3403_green[[i]]$sciName),fish_reef$scientificname)],'Key Largo',sep="-"),'.pdf',sep=''),width=8,height=6)
  plot(c(n.occ/tot.surveys)~year,data=REEF_3403_green[[i]],ylab='P(Occurring)',xlab='Year',bty='l',type='n',main=paste(fish_reef$commonname[match(unique(REEF_3403_green[[i]]$sciName),fish_reef$scientificname)],'Key Largo',sep="-"),ylim=c(0,1))
  lines(c(n.occ/tot.surveys)~year,data=REEF_3403_green[[i]],lwd=2,col='navy')
  points(c(n.occ/tot.surveys)~year,data=REEF_3403_green[[i]],pch=21,col='white',bg='navy',cex=1.5)
  lines(c(occ.sites/sites.surveyed)~year,data=REEF_3403_green[[i]],lwd=2,col='darkcyan')
  points(c(occ.sites/sites.surveyed)~year,data=REEF_3403_green[[i]],pch=21,col='white',bg='darkcyan',cex=1.5)
  
  m<- match(unique(REEF_3403_green[[i]]$sciName),REEF.green.exp.3403.sp)
  if(is.na(m)==F){
    lines(c(n.occ/tot.surveys)~year,data=REEF_3403_exp_green[[m]],lwd=2,col='coral4')
    points(c(n.occ/tot.surveys)~year,data=REEF_3403_exp_green[[m]],pch=21,col='white',bg='coral4',cex=1.5)
    lines(c(occ.sites/sites.surveyed)~year,data=REEF_3403_exp_green[[m]],lwd=2,col='darkgoldenrod4')
    points(c(occ.sites/sites.surveyed)~year,data=REEF_3403_exp_green[[m]],pch=21,col='white',bg='darkgoldenrod4',cex=1.5)
    
  }
  
  m2<- match(unique(REEF_3403_green[[i]]$sciName),fish_reef$scientificname)
  rvc_code<- fish_reef$rvc_code[m]
  m3<- match(rvc_code,rvc.green.3403.sp)
  if(is.na(m3)==F){
    lines(c(ssu.occ/nm)~YEAR,data=rvc_3403_green[[m3]],lwd=2,col='darkgray')
    points(c(ssu.occ/nm)~YEAR,data=rvc_3403_green[[m3]],pch=21,col='white',bg='darkgray',cex=1.5)
    lines(c(psu.occ/n)~YEAR,data=rvc_3403_green[[m3]],lwd=2,col='antiquewhite3')
    points(c(psu.occ/n)~YEAR,data=rvc_3403_green[[m3]],pch=21,col='white',bg='antiquewhite3',cex=1.5)
    
  }
  legend(x=2012,y=1.05,c('Nov+Exp surv freq.','Nov+Exp site freq.','Exp surv freq.','Exp site freq.','RVC surv freq.','RVC site freq.'),text.col=c('navy','darkcyan','coral4','darkgoldenrod4','darkgray','antiquewhite3'),bty='n',y.intersp = 0.7)
  dev.off(paste(paste(fish_reef$commonname[match(unique(REEF_3403_green[[i]]$sciName),fish_reef$scientificname)],'Key Largo',sep="-"),'.pdf',sep=''))
}
for(i in 101:150){
  par(xpd=T)
  pdf(file=paste(paste(fish_reef$commonname[match(unique(REEF_3403_green[[i]]$sciName),fish_reef$scientificname)],'Key Largo',sep="-"),'.pdf',sep=''),width=8,height=6)
  plot(c(n.occ/tot.surveys)~year,data=REEF_3403_green[[i]],ylab='P(Occurring)',xlab='Year',bty='l',type='n',main=paste(fish_reef$commonname[match(unique(REEF_3403_green[[i]]$sciName),fish_reef$scientificname)],'Key Largo',sep="-"),ylim=c(0,1))
  lines(c(n.occ/tot.surveys)~year,data=REEF_3403_green[[i]],lwd=2,col='navy')
  points(c(n.occ/tot.surveys)~year,data=REEF_3403_green[[i]],pch=21,col='white',bg='navy',cex=1.5)
  lines(c(occ.sites/sites.surveyed)~year,data=REEF_3403_green[[i]],lwd=2,col='darkcyan')
  points(c(occ.sites/sites.surveyed)~year,data=REEF_3403_green[[i]],pch=21,col='white',bg='darkcyan',cex=1.5)
  
  m<- match(unique(REEF_3403_green[[i]]$sciName),REEF.green.exp.3403.sp)
  if(is.na(m)==F){
    lines(c(n.occ/tot.surveys)~year,data=REEF_3403_exp_green[[m]],lwd=2,col='coral4')
    points(c(n.occ/tot.surveys)~year,data=REEF_3403_exp_green[[m]],pch=21,col='white',bg='coral4',cex=1.5)
    lines(c(occ.sites/sites.surveyed)~year,data=REEF_3403_exp_green[[m]],lwd=2,col='darkgoldenrod4')
    points(c(occ.sites/sites.surveyed)~year,data=REEF_3403_exp_green[[m]],pch=21,col='white',bg='darkgoldenrod4',cex=1.5)
    
  }
  
  m2<- match(unique(REEF_3403_green[[i]]$sciName),fish_reef$scientificname)
  rvc_code<- fish_reef$rvc_code[m]
  m3<- match(rvc_code,rvc.green.3403.sp)
  if(is.na(m3)==F){
    lines(c(ssu.occ/nm)~YEAR,data=rvc_3403_green[[m3]],lwd=2,col='darkgray')
    points(c(ssu.occ/nm)~YEAR,data=rvc_3403_green[[m3]],pch=21,col='white',bg='darkgray',cex=1.5)
    lines(c(psu.occ/n)~YEAR,data=rvc_3403_green[[m3]],lwd=2,col='antiquewhite3')
    points(c(psu.occ/n)~YEAR,data=rvc_3403_green[[m3]],pch=21,col='white',bg='antiquewhite3',cex=1.5)
    
  }
  legend(x=2012,y=1.05,c('Nov+Exp surv freq.','Nov+Exp site freq.','Exp surv freq.','Exp site freq.','RVC surv freq.','RVC site freq.'),text.col=c('navy','darkcyan','coral4','darkgoldenrod4','darkgray','antiquewhite3'),bty='n',y.intersp = 0.7)
  dev.off(paste(paste(fish_reef$commonname[match(unique(REEF_3403_green[[i]]$sciName),fish_reef$scientificname)],'Key Largo',sep="-"),'.pdf',sep=''))
}

for(i in 151:length(REEF.green.3403.sp)){
  par(xpd=T)
  pdf(file=paste(paste(fish_reef$commonname[match(unique(REEF_3403_green[[i]]$sciName),fish_reef$scientificname)],'Key Largo',sep="-"),'.pdf',sep=''),width=8,height=6)
  plot(c(n.occ/tot.surveys)~year,data=REEF_3403_green[[i]],ylab='P(Occurring)',xlab='Year',bty='l',type='n',main=paste(fish_reef$commonname[match(unique(REEF_3403_green[[i]]$sciName),fish_reef$scientificname)],'Key Largo',sep="-"),ylim=c(0,1))
  lines(c(n.occ/tot.surveys)~year,data=REEF_3403_green[[i]],lwd=2,col='navy')
  points(c(n.occ/tot.surveys)~year,data=REEF_3403_green[[i]],pch=21,col='white',bg='navy',cex=1.5)
  lines(c(occ.sites/sites.surveyed)~year,data=REEF_3403_green[[i]],lwd=2,col='darkcyan')
  points(c(occ.sites/sites.surveyed)~year,data=REEF_3403_green[[i]],pch=21,col='white',bg='darkcyan',cex=1.5)
  
  m<- match(unique(REEF_3403_green[[i]]$sciName),REEF.green.exp.3403.sp)
  if(is.na(m)==F){
    lines(c(n.occ/tot.surveys)~year,data=REEF_3403_exp_green[[m]],lwd=2,col='coral4')
    points(c(n.occ/tot.surveys)~year,data=REEF_3403_exp_green[[m]],pch=21,col='white',bg='coral4',cex=1.5)
    lines(c(occ.sites/sites.surveyed)~year,data=REEF_3403_exp_green[[m]],lwd=2,col='darkgoldenrod4')
    points(c(occ.sites/sites.surveyed)~year,data=REEF_3403_exp_green[[m]],pch=21,col='white',bg='darkgoldenrod4',cex=1.5)
    
  }
  
  m2<- match(unique(REEF_3403_green[[i]]$sciName),fish_reef$scientificname)
  rvc_code<- fish_reef$rvc_code[m]
  m3<- match(rvc_code,rvc.green.3403.sp)
  if(is.na(m3)==F){
    lines(c(ssu.occ/nm)~YEAR,data=rvc_3403_green[[m3]],lwd=2,col='darkgray')
    points(c(ssu.occ/nm)~YEAR,data=rvc_3403_green[[m3]],pch=21,col='white',bg='darkgray',cex=1.5)
    lines(c(psu.occ/n)~YEAR,data=rvc_3403_green[[m3]],lwd=2,col='antiquewhite3')
    points(c(psu.occ/n)~YEAR,data=rvc_3403_green[[m3]],pch=21,col='white',bg='antiquewhite3',cex=1.5)
    
  }
  legend(x=2012,y=1.05,c('Nov+Exp surv freq.','Nov+Exp site freq.','Exp surv freq.','Exp site freq.','RVC surv freq.','RVC site freq.'),text.col=c('navy','darkcyan','coral4','darkgoldenrod4','darkgray','antiquewhite3'),bty='n',y.intersp = 0.7)
  dev.off(paste(paste(fish_reef$commonname[match(unique(REEF_3403_green[[i]]$sciName),fish_reef$scientificname)],'Key Largo',sep="-"),'.pdf',sep=''))
}


#####Occupancy time-series####




#Rank.plot.pdf<- function(x,m){
#  table<- t(as.matrix(D[,10:14]))
#  colnames(table)<- D$Year
#  pdf(file=paste(paste(unique(D$comName),loc_full$geog[match(GZ,loc_full$geogid)],sep='-'),'.pdf',sep=''),width=8,height=6)
#  par(xpd=T,mar=c(4,4,3.5,3))
# barplot(table,ylab='Proportion of sightings in each category',xlab='Year',col=C,main=paste(unique(D$comName),loc_full$geog[match(GZ,loc_full$geogid)],sep=' - '))
# text(y=rep(1.02,nrow(D)),x=seq(0.7,31.8,length.out = 27),D$N.Surveys,cex=0.7)
 # text(y=1.02,x=34.5,'No. surveys',cex=0.8)
#text(y=0.91,x=34.5,'Rank\n Abundance',cex=0.8)
 # legend(x=32.5,y=0.9,rev(c('0','1','2','3','4')),cex=1.1,bty='n',text.col=rev(C))
  #dev.off(paste(paste(unique(D$comName),loc_full$geog[match(GZ,loc_full$geogid)],sep='-'),'.pdf',sep=''))
#}


for(i in 1:length(rvc.green.3403.sp)){
  spp<- filter(fish_reef,rvc_code==rvc.green.3403.sp[i])
  ts_comp[[i]]<- t(log(rvc_3403_green[[i]]$mean_ssu_abundance)) #extract logged time-series
  m<- match(spp$scientificname,REEF.green.3403.sp)
  if(is.na(m)==T){next}
  ts_comp[[i]]<- rbind(ts_comp[[i]],t(c(rep(NA,13),REEF_3403_green[[m]]$logAbund[1:26])))
  rownames(ts_comp[[i]])<- c(paste(spp$scientificname,'rvc',sep=":"),paste(spp$scientificname,'REEF',sep=":"))
  colnames(ts_comp[[i]])<- years
  
  Z_1=factor(c('rvc_reef','rvc_reef')) #Common state process for both time-series
  Z_2=factor(c('rvc','reef')) #Different state process for each time-series
  R_matrix[[i]]=matrix(list(0),nrow(ts_comp[[i]]),nrow(ts_comp[[i]])) #create empty Obs. error matrix
  diag(R_matrix[[i]])=c(rep("rvc",times=1),rep("reef",times=1))
  
  #tier 1 - common state process
  tier_1<- list(Z=Z_1,
                Q = 'diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 2 - Diff. U, same process error
  tier_2<- list(Z=Z_2,
                Q ='diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 3
  tier_3<- list(Z=Z_2,
                Q = 'diagonal and unequal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  fit_1<-  MARSS(ts_comp[[i]], model = tier_1,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem') 
  fit_2<-  MARSS(ts_comp[[i]], model = tier_2,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem')
  fit_3<-  MARSS(ts_comp[[i]], model = tier_3,control=list(maxit=30000,minit=500,conv.test.slope.tol=0.05),method='kem')
  
  
#pdf(file=paste('./timeseries/',fish_reef$commonname[match(rvc.green.3403.sp[i],fish_reef$rvc_code)],'.pdf',sep=''),width=8,height=6)
  plot(ts_comp[[i]][1,]~colnames(ts_comp[[i]]),type='b',col='navy',pch=21,bg='navy',bty='l',ylab='log(abundance)',xlab='Year',main=fish_reef$commonname[match(rvc.green.3403.sp[i],fish_reef$rvc_code)],ylim=c(min(na.omit(ts_comp[[i]][1,])),max(na.omit(ts_comp[[i]][2,]))))
  lines(c(ts_comp[[i]][2,])~c(colnames(ts_comp[[i]])),col='darkred')
  points(ts_comp[[i]][2,]~colnames(ts_comp[[i]]),pch=21,bg='darkred')
  if(which.min(c(fit_1$AICc,fit_2$AICc,fit_3$AICc))==1){
    lines(c(fit_1$states)~colnames(ts_comp[[i]]),col='darkgray',lty=5,lwd=3)
  }
  if(which.min(c(fit_1$AICc,fit_2$AICc,fit_3$AICc))==2){
    lines(c(fit_2$states[1,])~colnames(ts_comp[[i]]),col='navy',lty=5,lwd=3)
    lines(c(fit_2$states[2,])~colnames(ts_comp[[i]]),col='darkred',lty=5,lwd=3)
  }
  if(which.min(c(fit_1$AICc,fit_2$AICc,fit_3$AICc))==3){
    lines(c(fit_3$states[1,])~colnames(ts_comp[[i]]),col='navy',lty=5,lwd=3)
    lines(c(fit_3$states[2,])~colnames(ts_comp[[i]]),col='darkred',lty=5,lwd=3)
  }
#dev.off(paste('./timeseries/',fish_reef$commonname[match(rvc.green.3403.sp[i],fish_reef$rvc_code)],'.pdf',sep=''))
  print(i)
}

#Region 3408
REEF.3408.comb.green<- do.call(rbind, lapply(REEF_3408_green, data.frame, stringsAsFactors=FALSE))
REEF.green.3408.sp<- unique(REEF.3408.comb.green$comName)

rvc.3408.comb<- do.call(rbind, lapply(rvc_3408_green, data.frame, stringsAsFactors=FALSE))
rvc.green.3408.sp<- unique(rvc.3408.comb$Species)

ts_comp<- list()
years<- seq(1980,2018,by=1)
n<- NA
R_matrix<- list()

mars_3408<- data.frame(SP=NA,conv.1=NA,conv.2=NA,conv.3=NA,tier.1.AICc=NA,tier2.AICc=NA,tier3.AICc=NA,Q1=NA,Q2=NA,Q3.rvc=NA,Q3.reef=NA,U1=NA,U2.rvc=NA,U2.reef=NA,U3.rvc=NA,U3.reef=NA,R1.rvc=NA,R1.reef=NA,R2.rvc=NA,R2.reef=NA,R3.rvc=NA,R3.reef=NA,cName=NA,obs.grp=NA,mAbund.rvc=NA,mAbund.reef=NA,sdAbund.rvc=NA,sdAbund.reef=NA,dAIC1=NA,dAIC2=NA,dAIC3=NA)
for(i in 1:length(rvc.green.3408.sp)){
  spp<- filter(fish_reef,rvc_code==rvc.green.3408.sp[i])
  ts_comp[[i]]<- t(log(rvc_3408_green[[i]]$mean_ssu_abundance)) #extract logged time-series
  m<- match(spp$scientificname,REEF.green.3408.sp)
  if(is.na(m)==T){next}
  ts_comp[[i]]<- rbind(ts_comp[[i]],t(c(rep(NA,13),REEF_3403_green[[m]]$logAbund[1:26])))
  rownames(ts_comp[[i]])<- c(paste(spp$scientificname,'rvc',sep=":"),paste(spp$scientificname,'REEF',sep=":"))
  colnames(ts_comp[[i]])<- years
  
  Z_1=factor(c('rvc_reef','rvc_reef')) #Common state process for both time-series
  Z_2=factor(c('rvc','reef')) #Different state process for each time-series
  R_matrix[[i]]=matrix(list(0),nrow(ts_comp[[i]]),nrow(ts_comp[[i]])) #create empty Obs. error matrix
  diag(R_matrix[[i]])=c(rep("rvc",times=1),rep("reef",times=1))
  
  #tier 1 - common state process
  tier_1<- list(Z=Z_1,
                Q = 'diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 2 - Diff. U, same process error
  tier_2<- list(Z=Z_2,
                Q ='diagonal and equal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  #tier 3
  tier_3<- list(Z=Z_2,
                Q = 'diagonal and unequal',
                R = R_matrix[[i]],
                A = 'scaling',
                x0 = 'unequal',
                tinitx=0)
  
  
  fit_1<-  MARSS(ts_comp[[i]], model = tier_1,control=list(maxit=10000,minit=500,conv.test.slope.tol=0.05),method='kem') 
  fit_2<-  MARSS(ts_comp[[i]], model = tier_2,control=list(maxit=10000,minit=500,conv.test.slope.tol=0.05),method='kem') 
  fit_3<-  MARSS(ts_comp[[i]], model = tier_3,control=list(maxit=10000,minit=500,conv.test.slope.tol=0.05),method='kem') 
  params.1<- MARSSparamCIs(fit_1)
  params.2<- MARSSparamCIs(fit_2)
  params.3<- MARSSparamCIs(fit_3)
  
  
  mars_3408[i,1]=spp$scientificname
  mars_3408[i,2]=ifelse(is.null(fit_1$errors)==T,1,0)
  mars_3408[i,3]=ifelse(is.null(fit_2$errors)==T,1,0)  
  mars_3408[i,4]=ifelse(is.null(fit_3$errors)==T,1,0)  
  mars_3408[i,5]=fit_1$AICc
  mars_3408[i,6]=fit_2$AICc  
  mars_3408[i,7]=fit_3$AICc  
  mars_3408[i,8]=params.1$parMean[5]
  mars_3408[i,9]=params.2$parMean[5]
  mars_3408[i,10]=params.3$parMean[5]
  mars_3408[i,11]=params.3$parMean[6]
  mars_3408[i,12]=params.1$parMean[4]
  mars_3408[i,13]=params.2$parMean[3] 
  mars_3408[i,14]=params.2$parMean[4]
  mars_3408[i,15]=params.3$parMean[3]
  mars_3408[i,16]=params.3$parMean[4]
  mars_3408[i,17]=params.1$parMean[2] 
  mars_3408[i,18]=params.1$parMean[3]
  mars_3408[i,19]=params.2$parMean[1] 
  mars_3408[i,20]=params.2$parMean[2]
  mars_3408[i,21]=params.3$parMean[1] 
  mars_3408[i,22]=params.3$parMean[2]
  mars_3408[i,23]=spp$commonname
  mars_3408[i,24]=spp$group_assignment
  mars_3408[i,25]=exp(mean(na.omit(ts_comp[[i]][1,])))
  mars_3408[i,26]=exp(mean(na.omit(ts_comp[[i]][2,])))
  mars_3408[i,27]=sd(exp(na.omit(ts_comp[[i]][1,])))
  mars_3408[i,28]=sd(exp(na.omit(ts_comp[[i]][2,])))
  mars_3408[i,29]=mars_3408[i,5]-min(mars_3408[i,5:7])
  mars_3408[i,30]=mars_3408[i,6]-min(mars_3408[i,5:7])
  mars_3408[i,31]=mars_3408[i,7]-min(mars_3408[i,5:7])
  print(i)
}
mars_3408<- mars_3408[complete.cases(mars_3408),]
mars_3408$mod<- NA
for(i in 1:nrow(mars_3408)){
  if(mars_3408[i,29]==0){mars_3408$mod[i]=1}
  if(mars_3408[i,30]==0){mars_3408$mod[i]=2}
  if(mars_3408[i,31]==0){mars_3408$mod[i]=3}
}

write.csv(mars_3408,'MARSS_model_comparisons_3408_1980.csv')


####Comparison plots###
mars_3403<- read.csv('MARSS_model_comparisons_3403_1980.csv')

mars_3403_m1<- subset(mars_3403,mod==1) #19
mars_3403_m2<- subset(mars_3403,mod==2) #27
mars_3403_m3<- subset(mars_3403,mod==3) #11

#U - 1980_2018 Key Largo
plot(seq(-0.11,0.11,by=0.01)~c(seq(-0.11,0.11,by=0.01)),type='n',ylab='U (Reef)',xlab='U (RVC)',bty='l',main='Population trends - Key Largo (1980-2018)',ylim=c(-0.1,0.06))
abline(h=0,lty=5)
abline(v=0,lty=5)
points(mars_3403_m1$U1~c(mars_3403_m1$U1),pch=21,bg='navy',col='white',cex=1.5)
points(mars_3403_m2$U2.reef~mars_3403_m2$U2.rvc,pch=21,bg='goldenrod',col='white',cex=1.5)
points(mars_3403_m3$U3.reef~mars_3403_m3$U3.rvc,pch=21,bg='darkred',col='white',cex=1.5)
text(y=jitter(mars_3403_m1$U1,150),x=c(mars_3403_m1$U1),mars_3403_m1$cName,cex=0.8)
text(y=jitter(mars_3403_m2$U2.reef,150),x=c(mars_3403_m2$U2.rvc),mars_3403_m2$cName,cex=0.8)
text(y=jitter(mars_3403_m3$U3.reef,150),x=c(mars_3403_m3$U3.rvc),mars_3403_m3$cName,cex=0.8)

#correlation
U.rvc<- c(mars_3403_m1$U1,mars_3403_m2$U2.rvc,mars_3403_m3$U3.rvc)
U.reef<- c(mars_3403_m1$U1,mars_3403_m2$U2.reef,mars_3403_m3$U3.reef)
cor.test(U.rvc,U.reef)

#Q - 1980_2018 Key Largo
plot(seq(0,0.2,by=0.01)~c(seq(0,0.2,by=0.01)),type='n',ylab='Q (Reef)',xlab='Q (RVC)',bty='l',main='Process Error - Key Largo (1980-2018)',xlim=c(-0.01,0.12))
points(mars_3403_m1$Q1~c(mars_3403_m1$Q1),pch=21,bg='navy',col='white',cex=1.5)
points(mars_3403_m2$Q2~c(mars_3403_m2$Q2),pch=21,bg='goldenrod',col='white',cex=1.5)
points(mars_3403_m3$Q3.reef~mars_3403_m3$Q3.rvc,pch=21,bg='darkred',col='white',cex=1.5)
text(y=jitter(mars_3403_m1$Q1,80),x=c(mars_3403_m1$Q1),mars_3403_m1$cName,cex=0.8)
text(y=jitter(mars_3403_m2$Q2,80),x=c(mars_3403_m2$Q2),mars_3403_m2$cName,cex=0.8)
text(y=jitter(mars_3403_m3$Q3.reef,80),x=c(mars_3403_m3$Q3.rvc),mars_3403_m3$cName,cex=0.8)

Q.rvc<- c(mars_3403_m1$Q,mars_3403_m2$Q.rvc,mars_3403_m3$Q3.rvc)
Q.reef<- c(mars_3403_m1$Q,mars_3403_m2$Q,mars_3403_m3$Q3.reef)
cor.test(Q.rvc,Q.reef)

##1993 to 2018 data
mars_3403_93<- read.csv('MARSS_model_comparisons_3403_1993.csv')

mars_3403_93_m1<- subset(mars_3403_93,mod==1) #19
mars_3403_93_m2<- subset(mars_3403_93,mod==2) #27
mars_3403_93_m3<- subset(mars_3403_93,mod==3) #11

#U - 1980_2018 Key Largo
plot(seq(-0.18,0.15,by=0.01)~c(seq(-0.18,0.15,by=0.01)),type='n',ylab='U (Reef)',xlab='U (RVC)',bty='l',main='Population trends - Key Largo (1993-2018)')
abline(h=0,lty=5)
abline(v=0,lty=5)
points(mars_3403_93_m1$U1~c(mars_3403_93_m1$U1),pch=21,bg='navy',col='white',cex=1.5)
points(mars_3403_93_m2$U2.reef~mars_3403_93_m2$U2.rvc,pch=21,bg='goldenrod',col='white',cex=1.5)
points(mars_3403_93_m3$U3.reef~mars_3403_93_m3$U3.rvc,pch=21,bg='darkred',col='white',cex=1.5)
text(y=jitter(mars_3403_93_m1$U1,50),x=c(mars_3403_93_m1$U1),mars_3403_m1$cName,cex=0.8)
text(y=jitter(mars_3403_93_m2$U2.reef,50),x=c(mars_3403_93_m2$U2.rvc),mars_3403_93_m2$cName,cex=0.8)
text(y=jitter(mars_3403_93_m3$U3.reef,50),x=c(mars_3403_93_m3$U3.rvc),mars_3403_93_m3$cName,cex=0.8)

#correlation
U.rvc<- c(mars_3403_93_m1$U1,mars_3403_93_m2$U2.rvc,mars_3403_93_m3$U3.rvc)
U.reef<- c(mars_3403_93_m1$U1,mars_3403_93_m2$U2.reef,mars_3403_93_m3$U3.reef)
cor.test(U.rvc,U.reef)

#Q - 1980_2018 Key Largo
plot(seq(0,0.3,by=0.01)~c(seq(0,0.3,by=0.01)),type='n',ylab='Q (Reef)',xlab='Q (RVC)',bty='l',main='Process Error - Key Largo (1980-2018)',xlim=c(-0.01,0.3))
points(mars_3403_93_m1$Q1~c(mars_3403_93_m1$Q1),pch=21,bg='navy',col='white',cex=1.5)
points(mars_3403_93_m2$Q2~c(mars_3403_93_m2$Q2),pch=21,bg='goldenrod',col='white',cex=1.5)
points(mars_3403_93_m3$Q3.reef~mars_3403_93_m3$Q3.rvc,pch=21,bg='darkred',col='white',cex=1.5)
text(y=jitter(mars_3403_93_m1$Q1,80),x=c(mars_3403_93_m1$Q1),mars_3403_93_m1$cName,cex=0.8)
text(y=jitter(mars_3403_93_m2$Q2,80),x=c(mars_3403_93_m2$Q2),mars_3403_93_m2$cName,cex=0.8)
text(y=jitter(mars_3403_93_m3$Q3.reef,80),x=c(mars_3403_93_m3$Q3.rvc),mars_3403_93_m3$cName,cex=0.8)

Q.rvc<- c(mars_3403_93_m1$Q,mars_3403_93_m2$Q.rvc,mars_3403_93_m3$Q3.rvc)
Q.reef<- c(mars_3403_93_m1$Q,mars_3403_93_m2$Q,mars_3403_93_m3$Q3.reef)
cor.test(Q.rvc,Q.reef)

##Islamorada 1993 to 2018 data
mars_3404_93<- read.csv('MARSS_model_comparisons_3404.csv')

mars_3404_93_m1<- subset(mars_3404_93,mod==1) #19
mars_3404_93_m2<- subset(mars_3404_93,mod==2) #27
mars_3404_93_m3<- subset(mars_3404_93,mod==3) #11

#U - 1980_2018 Islamorada
plot(seq(-0.12,0.18,by=0.01)~c(seq(-0.12,0.18,by=0.01)),type='n',ylab='U (Reef)',xlab='U (RVC)',bty='l',main='Population trends - Islamorada (1993-2018)')
abline(h=0,lty=5)
abline(v=0,lty=5)
points(mars_3404_93_m1$U1~c(mars_3404_93_m1$U1),pch=21,bg='navy',col='white',cex=1.5)
points(mars_3404_93_m2$U2.reef~mars_3404_93_m2$U2.rvc,pch=21,bg='goldenrod',col='white',cex=1.5)
points(mars_3404_93_m3$U3.reef~mars_3404_93_m3$U3.rvc,pch=21,bg='darkred',col='white',cex=1.5)
text(y=jitter(mars_3404_93_m1$U1,50),x=c(mars_3404_93_m1$U1),mars_3404_93_m1$cName,cex=0.8)
text(y=jitter(mars_3404_93_m2$U2.reef,50),x=c(mars_3404_93_m2$U2.rvc),mars_3404_93_m2$cName,cex=0.8)
text(y=jitter(mars_3404_93_m3$U3.reef,50),x=c(mars_3404_93_m3$U3.rvc),mars_3404_93_m3$cName,cex=0.8)

#correlation
U.rvc<- c(mars_3404_93_m1$U1,mars_3404_93_m2$U2.rvc,mars_3404_93_m3$U3.rvc)
U.reef<- c(mars_3404_93_m1$U1,mars_3404_93_m2$U2.reef,mars_3404_93_m3$U3.reef)
cor.test(U.rvc,U.reef)

#Q - 1980_2018 Key Largo
plot(seq(0,0.25,by=0.01)~c(seq(0,0.25,by=0.01)),type='n',ylab='Q (Reef)',xlab='Q (RVC)',bty='l',main='Process Error - Islamorada (1993-2018)',xlim=c(-0.01,0.25))
points(mars_3404_93_m1$Q1~c(mars_3404_93_m1$Q1),pch=21,bg='navy',col='white',cex=1.5)
points(mars_3404_93_m2$Q2~c(mars_3404_93_m2$Q2),pch=21,bg='goldenrod',col='white',cex=1.5)
points(mars_3404_93_m3$Q3.reef~mars_3404_93_m3$Q3.rvc,pch=21,bg='darkred',col='white',cex=1.5)
text(y=jitter(mars_3404_93_m1$Q1,1),x=c(mars_3404_93_m1$Q1),mars_3404_93_m1$cName,cex=0.8)
text(y=jitter(mars_3404_93_m2$Q2,1),x=c(mars_3404_93_m2$Q2),mars_3404_93_m2$cName,cex=0.8)
text(y=jitter(mars_3404_93_m3$Q3.reef,1),x=c(mars_3404_93_m3$Q3.rvc),mars_3404_93_m3$cName,cex=0.8)

Q.rvc<- c(mars_3404_93_m1$Q1,mars_3404_93_m2$Q2,mars_3404_93_m3$Q3.rvc)
Q.reef<- c(mars_3404_93_m1$Q,mars_3404_93_m2$Q2,mars_3404_93_m3$Q3.reef)
cor.test(Q.rvc,Q.reef)

##Key West 1993 to 2018 data
mars_3408_93<- read.csv('MARSS_model_comparisons_3408.csv')

mars_3408_93_m1<- subset(mars_3408_93,mod==1) #19
mars_3408_93_m2<- subset(mars_3408_93,mod==2) #27
mars_3408_93_m3<- subset(mars_3408_93,mod==3) #11

#U - 1980_2018 Islamorada
plot(seq(-0.13,0.13,by=0.01)~c(seq(-0.13,0.13,by=0.01)),type='n',ylab='U (Reef)',xlab='U (RVC)',bty='l',main='Population trends - Key West (1993-2018)')
abline(h=0,lty=5)
abline(v=0,lty=5)
points(mars_3408_93_m1$U1~c(mars_3408_93_m1$U1),pch=21,bg='navy',col='white',cex=1.5)
points(mars_3408_93_m2$U2.reef~mars_3408_93_m2$U2.rvc,pch=21,bg='goldenrod',col='white',cex=1.5)
points(mars_3408_93_m3$U3.reef~mars_3408_93_m3$U3.rvc,pch=21,bg='darkred',col='white',cex=1.5)
text(y=jitter(mars_3408_93_m1$U1,50),x=c(mars_3408_93_m1$U1),mars_3408_93_m1$cName,cex=0.8)
text(y=jitter(mars_3408_93_m2$U2.reef,50),x=c(mars_3408_93_m2$U2.rvc),mars_3408_93_m2$cName,cex=0.8)
text(y=jitter(mars_3408_93_m3$U3.reef,50),x=c(mars_3408_93_m3$U3.rvc),mars_3408_93_m3$cName,cex=0.8)

#correlation
U.rvc<- c(mars_3408_93_m1$U1,mars_3408_93_m2$U2.rvc,mars_3408_93_m3$U3.rvc)
U.reef<- c(mars_3408_93_m1$U1,mars_3408_93_m2$U2.reef,mars_3408_93_m3$U3.reef)
cor.test(U.rvc,U.reef)

#Q - 1980_2018 Key Largo
plot(seq(0,0.5,by=0.01)~c(seq(0,0.5,by=0.01)),type='n',ylab='Q (Reef)',xlab='Q (RVC)',bty='l',main='Process Error - Key West (1993-2018)',xlim=c(-0.03,0.5))
points(mars_3408_93_m1$Q1~c(mars_3408_93_m1$Q1),pch=21,bg='navy',col='white',cex=1.5)
points(mars_3408_93_m2$Q2~c(mars_3408_93_m2$Q2),pch=21,bg='goldenrod',col='white',cex=1.5)
points(mars_3408_93_m3$Q3.reef~mars_3408_93_m3$Q3.rvc,pch=21,bg='darkred',col='white',cex=1.5)
text(y=jitter(mars_3408_93_m1$Q1,1),x=c(mars_3408_93_m1$Q1),mars_3408_93_m1$cName,cex=0.8)
text(y=jitter(mars_3408_93_m2$Q2,1),x=c(mars_3408_93_m2$Q2),mars_3408_93_m2$cName,cex=0.8)
text(y=jitter(mars_3408_93_m3$Q3.reef,1),x=c(mars_3408_93_m3$Q3.rvc),mars_3408_93_m3$cName,cex=0.8)

Q.rvc<- c(mars_3408_93_m1$Q1,mars_3408_93_m2$Q2,mars_3408_93_m3$Q3.rvc)
Q.reef<- c(mars_3408_93_m1$Q1,mars_3408_93_m2$Q2,mars_3408_93_m3$Q3.reef)
cor.test(Q.rvc,Q.reef)



##

plot(seq(-0.15,0.1,by=0.01)~c(seq(-0.15,0.1,by=0.01)),type='n',ylab='U (Reef)',xlab='U (RVC)',bty='l',main='Population trends - Islamorada',ylim=c(-0.05,0.08),xlim=c(-0.07,0.1))
points(mars_3404_m1$U1~c(mars_3404_m1$U1),pch=21,bg='navy',col='white',cex=1.5)
points(mars_3404_m2$U2.reef~mars_3404_m2$U2.rvc,pch=21,bg='goldenrod',col='white',cex=1.5)
points(mars_3404_m3$U3.reef~mars_3404_m3$U3.rvc,pch=21,bg='darkred',col='white',cex=1.5)
text(y=jitter(mars_3404_m1$U1,20),x=c(mars_3404_m1$U1),mars_3404_m1$cName,cex=0.8)
text(y=jitter(mars_3404_m2$U2.reef,20),x=c(mars_3404_m2$U2.rvc),mars_3404_m2$cName,cex=0.8)
text(y=jitter(mars_3404_m3$U3.reef,20),x=c(mars_3404_m3$U3.rvc),mars_3404_m3$cName,cex=0.8)

plot(seq(0,0.3,by=0.01)~c(seq(0,0.3,by=0.01)),type='n',ylab='Q (Reef)',xlab='Q (RVC)',bty='l',main='Process Error - Islamorada',xlim=c(-0.01,0.3))
points(mars_3404_m1$Q1~c(mars_3404_m1$Q1),pch=21,bg='navy',col='white',cex=1.5)
points(mars_3404_m2$Q2~c(mars_3404_m2$Q2),pch=21,bg='goldenrod',col='white',cex=1.5)
points(mars_3404_m3$Q3.reef~mars_3404_m3$Q3.rvc,pch=21,bg='darkred',col='white',cex=1.5)
text(y=jitter(mars_3404_m1$Q1,1),x=c(mars_3404_m1$Q1),mars_3404_m1$cName,cex=0.8)
text(y=jitter(mars_3404_m2$Q2,1),x=c(mars_3404_m2$Q2),mars_3404_m2$cName,cex=0.8)


mars_3408<- read.csv('MARSS_model_comparisons_3408_1994.csv')
mars_3408$mod<- NA
for(i in 1:nrow(mars_3408)){
  if(mars_3408[i,30]==0){mars_3408$mod[i]=1}
  if(mars_3408[i,31]==0){mars_3408$mod[i]=2}
  if(mars_3408[i,32]==0){mars_3408$mod[i]=3}
}

mars_3408_m1<- subset(mars_3408,mod==1) #35
mars_3408_m2<- subset(mars_3408,mod==2) #26
mars_3408_m3<- subset(mars_3408,mod==3) #9

plot(seq(-0.15,0.1,by=0.01)~c(seq(-0.15,0.1,by=0.01)),type='n',ylab='U (Reef)',xlab='U (RVC)',bty='l',main='Population trends - Key West',ylim=c(-0.15,0.12),xlim=c(-0.15,0.12))
points(mars_3408_m1$U1~c(mars_3408_m1$U1),pch=21,bg='navy',col='white',cex=1.5)
points(mars_3408_m2$U2.reef~mars_3408_m2$U2.rvc,pch=21,bg='goldenrod',col='white',cex=1.5)
points(mars_3408_m3$U3.reef~mars_3408_m3$U3.rvc,pch=21,bg='darkred',col='white',cex=1.5)
text(y=jitter(mars_3408_m1$U1,20),x=c(mars_3408_m1$U1),mars_3408_m1$cName,cex=0.8)
text(y=jitter(mars_3408_m2$U2.reef,20),x=c(mars_3408_m2$U2.rvc),mars_3408_m2$cName,cex=0.8)
text(y=jitter(mars_3408_m3$U3.reef,20),x=c(mars_3408_m3$U3.rvc),mars_3408_m3$cName,cex=0.8)


plot(seq(0,0.4,by=0.01)~c(seq(0,0.4,by=0.01)),type='n',ylab='Q (Reef)',xlab='Q (RVC)',bty='l',main='Process Error - Key West',xlim=c(-0.01,0.4))
points(mars_3408_m1$Q1~c(mars_3408_m1$Q1),pch=21,bg='navy',col='white',cex=1.5)
points(mars_3408_m2$Q2~c(mars_3408_m2$Q2),pch=21,bg='goldenrod',col='white',cex=1.5)
points(mars_3408_m3$Q3.reef~mars_3408_m3$Q3.rvc,pch=21,bg='darkred',col='white',cex=1.5)
text(y=jitter(mars_3408_m1$Q1,20),x=c(mars_3408_m1$Q1),mars_3408_m1$cName,cex=0.8)
text(y=jitter(mars_3408_m2$Q2,20),x=c(mars_3408_m2$Q2),mars_3408_m2$cName,cex=0.8)
text(y=jitter(mars_3408_m3$Q3.reef,20),x=c(mars_3408_m3$Q3.rvc),mars_3408_m3$cName,cex=0.8)









for(i in 1:length(RVC.green.4303.sp)){
    #Data from RVC time-series 
    m_rvc<- match(rvc.green.4303.sp[i],fish_reef$rvc_code) #Find species
    ts_comp[[i]]<- t(log(rvc_3403_green[[i]]$mean_ssu_abundance)) #extract logged time-series
    #Find observation group for that species to include in MARSS
    obs<- fish_reef$group_assignment[m_rvc]
    fish_obs<- filter(fish_reef, group_assignment %in% obs) #Find all species in observation group
    fish_obs<- filter(fish_obs, rvc_code!=rvc.green.4303.sp[i]) #Only keep non-target species
    m_obs_rvc<- na.omit(match(fish_obs$rvc_code,rvc.green.4303.sp)) #Match other species in green list
    if(length(m_obs_rvc)>5){
      m_obs_rvc<- m_obs_rvc[sample(length(m_obs_rvc),5)]
    }
    for(x in 1:length(m_obs_rvc)){ #add each other time-series from the observation group
      ts_comp[[i]]<- rbind(ts_comp[[i]],t(log(rvc_3403_green[[m_obs_rvc[x]]]$mean_ssu_abundance)))
    }#add rvc time-series for obs error
    
    ##add REEF time-series
    m_reef<- match(fish_reef$scientificname[m_rvc],REEF.green.4303.sp) #Find species in reef data
    ts_comp[[i]]<- rbind(ts_comp[[i]],t(REEF_3403_green[[m_reef]]$logAbund[1:26])) #extract ts from 1993 to 2018
    colnames(ts_comp[[i]])<- years #columns = years
    
    #Add observation groups from REEF
    m_obs_reef<- na.omit(match(fish_obs$scientificname,REEF.green.4303.sp)) #Match other species in REEF green list
    if(length(m_obs_reef)>5){
      m_obs_reef<- m_obs_reef[sample(length(m_obs_reef),5)]
    }
    
    for(z in 1:length(m_obs_reef)){
      ts_comp[[i]]<- rbind(ts_comp[[i]],t(REEF_3403_green[[m_obs_reef[z]]]$logAbund[1:26])) #add rvc time-series for obs error
    }
    rownames(ts_comp[[i]])<- c(paste(fish_reef$scientificname[m_rvc],'rvc',sep=":"),rep('obs_grp:rvc',length(m_obs_rvc)),paste(fish_reef$scientificname[m_rvc],'reef',sep=":"),rep('obs_grp:reef',length(m_obs_reef)))
    n[i] = nrow(ts_comp[[i]]) - 1 # obtain n
    
    R_matrix[[i]]=matrix(list(0),nrow(ts_comp[[i]]),nrow(ts_comp[[i]])) #create empty Obs. error matrix
    diag(R_matrix[[i]])=c(rep("rvc",times=length(m_obs_rvc)+1),rep("reef",times=length(m_obs_reef)+1))
 #   U_matrix_1[[i]]=matrix(c(paste('u',seq(1:nrow(ts_comp[[i]])),sep='')),nrow(ts_comp[[i]]),1) #Trends are different for reef/rvc
    
  #  U_matrix_2[[i]]=matrix(list(0),nrow(ts_comp[[i]]),nrow(ts_comp[[i]])) #create empty Obs. error matrix
  #  diag(U_matrix_2[[i]])<- c(paste('u',seq(1,nrow(ts_comp[[i]]),by=1),sep=''))
  #  U_matrix_2[[i]][length(m_obs_rvc)+2,length(m_obs_rvc)+2]<- 'u1' #Set REEF time-series to have same U as RVC
    
 #   U_matrix_2[[i]]=matrix(c('u1',paste('u',seq(2,c(1+length(m_obs_rvc)),by=1),sep=''),'u1',paste('u',c(c(length(m_obs_rvc)+3):nrow(ts_comp[[i]])),sep='')),nrow(ts_comp[[i]]),1) #trends are the same for REEF/RVC
    Z_1=factor(c('rvc',paste('rvc',seq(1,length(m_obs_rvc),by=1),sep=""),'reef',paste('reef',seq(1,length(m_obs_reef),by=1),sep=""))) #Different process error for each time-series
    Z_2=factor(c('rvc_reef',paste('rvc',seq(1,length(m_obs_rvc),by=1),sep=""),'rvc_reef',paste('reef',seq(1,length(m_obs_reef),by=1),sep=""))) #Common process error for both time-series

    
    
    #tier 1
    mod_2<- list(Z=Z_2,
                 Q = 'diagonal and equal',
                 R = R_matrix[[i]],
                 A = 'scaling',
                 x0 = 'unequal',
                 tinitx=0)
    
    #tier 2
    mod_1a<- list(Z=Z_1,
                  Q ='diagonal and equal',
                  R = R_matrix[[i]],
                  A = 'scaling',
                  x0 = 'unequal',
                  tinitx=0)
    
    #tier 3
    mod_1<- list(Z=Z_1,
                 Q = 'diagonal and unequal',
                 R = R_matrix[[i]],
                 A = 'scaling',
                 x0 = 'unequal',
                 tinitx=0)
    
    
    
    #Model 1: Different U, different Q
    mod_1<- list(Z=Z_1,
                 Q = 'diagonal and unequal',
                 R = R_matrix[[i]],
                 A = 'scaling',
                 x0 = 'unequal',
                 tinitx=0)
    
    #Model 1a: Different U, all same Q
    mod_1a<- list(Z=Z_1,
                 Q ='diagonal and equal',
                 R = R_matrix[[i]],
                 A = 'scaling',
                 x0 = 'unequal',
                 tinitx=0)
    
    #Model 2: Same trends, same process error
    mod_2<- list(Z=Z_2,
                 Q = 'diagonal and unequal',
                 R = R_matrix[[i]],
                 A = 'scaling',
                 x0 = 'unequal',
                 tinitx=0)
    
    #Model 2a: Different U, all same Q
    mod_2a<- list(Z=Z_2,
                  Q = 'diagonal and equal',
                  R = R_matrix[[i]],
                  A = 'scaling',
                  x0 = 'unequal',
                  tinitx=0)
    
 
    fit_1<- MARSS(ts_comp[[i]], model = mod_1,control=list(maxit=50000)) 
    fit_2<- MARSS(ts_comp[[i]], model = mod_2,control=list(maxit=50000)) 
    fit_1a<- MARSS(ts_comp[[i]], model = mod_1a,control=list(maxit=50000)) 
    fit_2a<- MARSS(ts_comp[[i]], model = mod_2a,control=list(maxit=50000)) 
    
    
    ts_comp_trim<- rbind(ts_comp[[i]][1,],ts_comp[[i]][7,])
    Z_1=factor(c('rvc','reef')) #Different process error for each time-series
    Z_2=factor(c('rvc_reef','rvc_reef')) #Common process error for both time-series
    R_matrix[[i]]=matrix(list(0),nrow(ts_comp_trim),nrow(ts_comp_trim)) #create empty Obs. error matrix
    diag(R_matrix[[i]])=c(rep("rvc",times=1),rep("reef",times=1))
    
    fit_1t<- MARSS(ts_comp_trim, model = mod_1,control=list(maxit=50000)) 
    fit_2t<- MARSS(ts_comp_trim, model = mod_2,control=list(maxit=50000)) 
    fit_1a_t<- MARSS(ts_comp_trim, model = mod_1a,control=list(maxit=50000)) 
    fit_2a_t<- MARSS(ts_comp_trim, model = mod_2a,control=list(maxit=50000)) 
    
  
    fit_3<- MARSS(ts_comp[[i]], model = mod_3,control=list(maxit=20000)) 
    fit_4<- MARSS(ts_comp[[i]], model = mod_4,control=list(maxit=20000)) 
    
                         
}

R_matrix=matrix(list(0),10,10)

mod_list[[6]] <- list(Z = factor(c("PrincessParrotfish:Loc1","QueenParrotfish:Loc1",
                                   "RedbandParrotfish:Loc1", "StoplightParrotfish:Loc1",
                                   "StripedParrotfish:Loc1","PrincessParrotfish:Loc2","QueenParrotfish:Loc2",
                                   "RedbandParrotfish:Loc2", "StoplightParrotfish:Loc2",
                                   "StripedParrotfish:Loc2")), #all species at each location come frome different state processes
                      Q = "diagonal and unequal", 
                      R = R_matrix,
                      A = "scaling",
                      U = "zero",
                      x0 = "unequal", 
                      tinitx = 0) 


REEF.4304.comb.green<- do.call(rbind, lapply(REEF_3404_green, data.frame, stringsAsFactors=FALSE))
REEF.trend.spp.4304<- unique(REEF.4304.comb.green$comName)
qp.4304<- filter(REEF.4304.comb.green, comName=='Scarus vetula')
acf(qp.4304$meanAbund,main='Queen Parrotfish - 4304')

REEF.4308.comb.green<- do.call(rbind, lapply(REEF_3408_green, data.frame, stringsAsFactors=FALSE))
REEF.trend.spp.4308<- unique(REEF.4308.comb.green$comName)
qp.4308<- filter(REEF.4308.comb.green, comName=='Scarus vetula')
acf(qp.4308$meanAbund,main='Queen Parrotfish - 4308')


###Comparisons
mars_3403_93_novexp<- read.csv('MARSS_model_comparisons_3403_1993_nov+expert.csv')
mars_3403_93_exp<-read.csv('MARSS_model_comparisons_3403_1993_expert.csv')
m<- match(mars_3403_93_exp$SP,mars_3403_93_novexp$SP)  
comp_93_q<- cbind(mars_3403_93_exp$Q3.reef,mars_3403_93_novexp$Q3.reef[m])
plot(comp_93_q[,1],comp_93_q[,2])
nrow(subset(comp_93_q,comp_93_q[,1]>0))
nrow(subset(comp_93_q,comp_93_q[,2]>0))

mars_3403_80_novexp<- read.csv('MARSS_model_comparisons_3403_1980.csv')
mars_3403_80_exp<-read.csv('MARSS_model_comparisons_3403_1980_expert_data.csv')
names(mars_3403_80_novexp)<- names(mars_3403_80_exp)
m2<- match(mars_3403_80_novexp$SP,mars_3403_80_exp$SP)

comp_80_q<- cbind(mars_3403_80_novexp$Q3.reef,mars_3403_80_exp$Q3.reef[m2])
plot(comp_80_q[,1],comp_80_q[,2])

nrow(subset(comp_80_q,comp_80_q[,1]>0))
nrow(subset(comp_80_q,comp_80_q[,2]>0))

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


###Stan models
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
  real beta;
  //effect of covariates
  
  //deviations from intercept
  real a_hab[N_hab]; //deviation between habitats
  real a_yr[N_yr]; //deviation between years
 
  
  //st dev on the deviations
  real<lower = 0> sigma_hab; //sigma on habitat
  real<lower = 0> sigma_yr; //sigma on year
}

model{
  //priors
  alpha ~ normal(0,10);
  beta ~ normal(0,2);

  //standard deviations
  sigma_hab ~ cauchy(0, 5);
  sigma_yr ~ cauchy(0, 5);
  
  //varying intercepts
  a_hab ~ normal(0, sigma_hab);
  a_yr ~student_t(5, 0, sigma_yr);
  
  for(i in 1:N){
    y[i] ~ bernoulli_logit(alpha + a_yr[year[i]] + a_hab[hab_class[i]]+ X[i,]*beta);
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
   int<lower=0> site_hab[N_site]; //site within habitat
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
  real a_yr[N_yr]; //deviation between years
  real a_dv[N_dv]; //deviation between divers
  real dev_hab[N_hab]; //deviation between habitats
  real dev_site[N_site]; ///deviation between sites
  
  //covariates for effort variables
  vector[K] beta;
  
  //st dev on the deviations
  real<lower = 0> sigma_hab;
  real<lower = 0> sigma_yr;
  real<lower = 0> sigma_site;
  real<lower = 0> sigma_dv;
}

transformed parameters{
  real a_hab[N_hab]; //deviation between habitats
  real a_site[N_site]; //deviation between sites
  
 for(h in 1:N_hab){
    a_hab[h] = alpha + dev_hab[h];
  }
  
  for(s in 1:N_site){
    a_site[s] = a_hab[site_hab[s]] + dev_site[s];
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
  dev_hab ~ normal(0, sigma_hab);
  a_yr ~ student_t(5, 0, sigma_yr);
  dev_site ~ normal(0, sigma_site);
  a_dv ~ normal(0, sigma_dv);
  
  for(i in 1:N){
    y[i] ~ bernoulli_logit(a_site[site[i]] + a_yr[year[i]] +  a_dv[diver[i]] + X[i,]*beta);
  }
}
"


####Batch through species and retain the mean estimate (intercept) from both models
for(i in 1:nrow(fish_reef)){
  spp_rvc<- rvc_occs[[i]]
  spp_reef<- reef_occs[[i]]
  
  X_rvc<- matrix(data=c(scale(as.numeric(spp_rvc$DEPTH))))
  X_reef<- matrix(data=c(scale(as.numeric(spp_reef$btime)),scale(as.numeric(spp_reef$visibility)),scale(as.numeric(spp_reef$current)),scale(as.numeric(spp_reef$averagedepth))),ncol=4,nrow=nrow(spp_reef))
  
  
  rvc_mod<- rstan::stan(model_code = logit_test_rvc, data = list(y = spp_rvc$occ, 
                                                                   N = nrow(spp_rvc),
                                                                   N_hab = length(unique(spp_rvc$HAB_CD2)),
                                                                   hab_class=as.numeric(factor(spp_rvc$HAB_CD2)),
                                                                   N_yr =length(unique(spp_rvc$YEAR)),
                                                                   year=as.numeric(factor(spp_rvc$YEAR)),
                                                                   X=X_rvc,
                                                                   K=ncol(X_rvc)),
                          pars = c("alpha", "a_hab",'sigma_hab','sigma_yr','a_yr','beta'),
                          control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 200, chains = 4, iter = 700, thin = 1)
  
  reef_mod<- rstan::stan(model_code = logit_test_reef, data = list(y = spp_reef$occ, 
                                                                      N = nrow(spp_reef),
                                                                      N_hab = length(unique(spp_reef$hab_class2)),
                                                                      hab_class=as.numeric(factor(spp_reef$hab_class2)),
                                                                      N_yr =length(unique(spp_reef$year)),
                                                                      year=as.numeric(factor(spp_reef$year)),
                                                                      site=as.numeric(factor(spp_reef$geogr)),
                                                                      N_site=length(unique(spp_reef$geogr)),
                                                                      site_hab=as.numeric(factor(distinct(spp_reef,geogr,.keep_all = T)$hab_class2)),
                                                                      diver=as.numeric(factor(spp_reef$fish_memberid)),
                                                                      N_dv=length(unique(spp_reef$fish_memberid)),
                                                                      K=ncol(X_reef),
                                                                      X=X_reef),
                            pars = c("alpha",'a_yr','dev_hab','sigma_hab','sigma_yr','sigma_dv','sigma_site','beta'),
                            control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 200, chains = 4, iter = 700, thin = 1)
  
}

posterior <- as.array(reef_mod)
dimnames(posterior)
color_scheme_set("viridis")
mcmc_areas(
  posterior,
  pars = c("alpha", paste('a_hab[',seq(1:4),']',sep=''),'sigma_hab','sigma_yr','sigma_dv','sigma_site'),
  prob = 0.95, # 80% intervals
  prob_outer = 1, # 99%
  point_est = "median"
)
mcmc_dens_overlay(posterior, c("alpha", paste('a_yr[',seq(1:23),']',sep=''),'sigma_hab','sigma_yr'))
mcmc_pairs(posterior, pars = c("alpha", paste('a_hab[',seq(1:4),']',sep=''),'sigma_hab',paste('a_yr[',seq(1:4),']',sep='')),
           off_diag_args = list(size = 1.5))
mcmc_trace(posterior, pars = c(paste('a_hab[',seq(1:4),']',sep='')))

           
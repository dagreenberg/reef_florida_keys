rm(list=ls())
setwd("C:/Users/14388/Desktop/reef_florida_keys_data")
library(rstan);library(loo);library(bayesplot);library(dplyr);library(tidyverse);library(stringr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

####Functions####
rvc_occurrence = function(x,GZ,sp){
  x$SSU_YEAR<- paste(x$PRIMARY_SAMPLE_UNIT,x$STATION_NR,x$YEAR,sep='_')
  x1= x %>% subset(region.id==GZ)
  x1= complete(x1,SSU_YEAR,nesting(SPECIES_CD),fill=list(NUM=0))
  x2= subset(x1,SPECIES_CD==sp)
  x3=  x2 %>% dplyr::group_by(SSU_YEAR) %>%
    dplyr::summarise(NUM.total=sum(NUM),occ=NA) %>% #Sums up the number of counts
    mutate(occ=ifelse(NUM.total>0,1,0)) %>% arrange(SSU_YEAR) #Also scores presence/absence at the SSU level
  x3[,4:31]<- x2[match(x3$SSU_YEAR,x2$SSU_YEAR),2:29]
  x3<- transform(x3,psu_id=match(LAT_LON,unique(LAT_LON)))
  
}

occ_ts = function(x){
  x %>% group_by(YEAR) %>% summarize(n.occ=sum(occ),n.surv=n(),p.occ=n.occ/n.surv)
}

comp_plot = function(ts,mod_mat,sp,GZ){
  plot(ts$p.occ~ts$YEAR,type='n',xlab='Year',ylab='Probability of occurrence',ylim=c(0,1),bty='l',main=paste(sp,GZ,sep=' - '))
  lines(ts$p.occ~ts$YEAR,lwd=2,col='darkblue')
  points(ts$p.occ~ts$YEAR,pch=21,col='darkblue',bg='white',cex=1.5,lwd=1.5)
  points(mod_mat$median~ts$YEAR,pch=21,col=adjustcolor('dodgerblue3',alpha.f = 0.8),cex=1.2,lwd=1.2)
  lines(mod_mat$median~ts$YEAR,lwd=1.5,col=adjustcolor('dodgerblue3',alpha.f = 0.8))
  x.polygon <- c(ts$YEAR, rev(ts$YEAR)) # Define a polygon x value for adding to a plot
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
#
fk_99_18<- read.csv("./RVC/Florida_Keys_RVC.csv") #RVC survey data from 1999 to 2018
fk_79_98<- read.csv("./RVC/Pre1999_edit.txt") #RVC survey dating from 1979 to 1998

#Synonymize column names
colnames(fk_79_98)<- c(colnames(fk_99_18)[1],colnames(fk_99_18)[3:5],colnames(fk_99_18)[2],colnames(fk_99_18)[6],colnames(fk_99_18)[7:18],colnames(fk_99_18)[19],colnames(fk_99_18)[20]) 
#Join together datasets
fk_79_18<- full_join(fk_99_18,fk_79_98)


library(sp)
reef_geog<- read.csv("./REEF Tropical Western Atlantic/TWAgeog.csv")
#site_surveys<-REEF %>% group_by(geogr) %>% summarize(n=n_distinct(formid))
#m<- match(reef_geog$geogid,site_surveys$geogr)
#reef_geog$no.surveys<- site_surveys$n[m]
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

fk_93_18<- subset(fk_79_18,YEAR>=1993) #Subset for the dataset from 1993 to match the first year of REEF surveys



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
logit_test_yr2<-"data{
  int<lower=1> N;//number of observations (SSU surveys)
  int y[N]; //presence or absence on each survey
  int<lower=0> N_hab; //number of habitat classes
  int<lower=1,upper=N_hab> hab_class[N]; // vector of habitat class identities
  int<lower=0> N_yr; //number of years
  int<lower=1,upper=N_yr> year[N]; // vector of year
}
parameters {
  //global intercept
  real alpha;
  
  //
  
  //deviations from intercept
  real a_hab[N_hab]; //deviation between habitats
  real a_yr[N_yr]; //deviation between months
  
  //st dev on the deviations
  real<lower = 0> sigma_hab;
  real<lower = 0> sigma_yr;
}

transformed parameters{
  vector[N] eta;
  
  for(n in 1:N){
    eta[n] = alpha + a_hab[hab_class[n]] + a_yr[year[n]];
  }
}

model{
  //priors
  alpha ~ normal(0,10);
  
  //standard deviations
  sigma_hab ~ cauchy(0, 5);
  sigma_yr ~ cauchy(0, 5);
  
  //varying intercepts
  a_hab ~ normal(0, sigma_hab);
  a_yr ~ student_t(5, 0, sigma_yr);
 
  for(i in 1:N){
    y[i] ~ bernoulli_logit(eta[i]);
  }
}
"

logit_test_hab<-"data{
  int<lower=1> N;//number of observations (SSU surveys)
  int y[N]; //presence or absence on each survey
  int<lower=0> N_hab; //number of habitat classes
  int<lower=1,upper=N_hab> hab_class[N]; // vector of habitat class identities
}
parameters {
  //global intercept
  real alpha;
  
  //deviations from intercept
  real a_hab[N_hab]; //deviation between habitats

  //st dev on the deviations
  real<lower = 0> sigma_hab;
}

transformed parameters{
  vector[N] eta;
  
  for(n in 1:N){
    eta[n] = alpha + a_hab[hab_class[n]];
  }
}

model{
  //priors
  alpha ~ normal(0,10);
  
  //standard deviations
  sigma_hab ~ cauchy(0, 5);
  
  //varying intercepts
  a_hab ~ normal(0, sigma_hab);
 
  for(i in 1:N){
    y[i] ~ bernoulli_logit(eta[i]);
  }
}
"

logit_test_psu<-"data{
  int<lower=1> N;//number of observations (SSU surveys)
  int y[N]; //presence or absence on each survey
  int<lower=0> N_psu; //number of primary sample units
  int<lower=0> N_hab; //number of habitat classes
  int<lower=1,upper=N_psu> psu_id[N]; // vector of psu identities
  int<lower=1,upper=N_hab> hab_class[N]; // vector of habitat class identities
}
parameters {
  //global intercept
  real alpha;
  
  //
  
  //deviations from intercept
  real a_hab[N_hab]; //deviation between habitats
  real a_psu[N_psu]; //deviation between months
  
  //st dev on the deviations
  real<lower = 0> sigma_hab;
  real<lower = 0> sigma_psu;
}

transformed parameters{
  vector[N] eta;
  
  for(n in 1:N){
    eta[n] = alpha + a_hab[hab_class[n]] + a_psu[psu_id[n]];
  }
}

model{
  //priors
  alpha ~ normal(0,10);
  
  //standard deviations
  sigma_hab ~ cauchy(0, 5);
  sigma_psu ~ cauchy(0, 5);
  
  //varying intercepts
  a_hab ~ normal(0, sigma_hab);
  a_psu ~ normal(0, sigma_psu);
 
  for(i in 1:N){
    y[i] ~ bernoulli_logit(eta[i]);
  }
}
"
psu_hab<- distinct(blue_angel,psu_id,.keep_all=T)
psu_in_hab<- as.numeric(factor(psu_hab$HABITAT_CD))


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

test_1_y<- rstan::stan(model_code = logit_test_yr, data = list(y = blue_angel$occ, 
                                                     N = nrow(blue_angel),
                                                     N_hab = length(unique(blue_angel$HABITAT_CD)),
                                                     hab_class=as.numeric(factor(blue_angel$HABITAT_CD)),
                                                     N_mth = length(unique(blue_angel$MONTH)), 
                                                     month=as.numeric(factor(blue_angel$MONTH)),
                                                     N_yr =length(unique(blue_angel$YEAR)),
                                                     year=as.numeric(factor(blue_angel$YEAR))),
                                                    pars = c("alpha", "a_hab","a_mth",'a_yr','sigma_hab','sigma_mth','sigma_yr'),
                                                  control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 400, chains = 4, iter = 2000, thin = 1)

blue_angel<- rvc_occurrence(fk_93_18,sp='HOL BERM',GZ=3403) #blue angelfish in key largo


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

comp_plot(ts=occ_ts(blue_angel),mod_mat = a_mat,sp='Blue Angelfish',GZ='Key Largo')

gray_angel<- rvc_occurrence(fk_93_18,sp='POM PARU',GZ=3403) #blue angelfish in key largo

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

comp_plot(ts=occ_ts(gray_angel),mod_mat = a_mat,sp='Gray Angelfish',GZ='Key Largo')

gray_angel<- rvc_occurrence(fk_93_18,sp='POM PARU',GZ=3408) #blue angelfish in key west
test_1_y2_g2<- rstan::stan(model_code = logit_test_yr2, data = list(y = gray_angel$occ, 
                                                                   N = nrow(gray_angel),
                                                                   N_hab = length(unique(gray_angel$HABITAT_CD)),
                                                                   hab_class=as.numeric(factor(gray_angel$HABITAT_CD)),
                                                                   N_yr =length(unique(gray_angel$YEAR)),
                                                                   year=as.numeric(factor(gray_angel$YEAR))),
                          pars = c("alpha", "a_hab",'a_yr','sigma_hab','sigma_yr'),
                          control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 400, chains = 4, iter = 2000, thin = 1)

ga_params<- rstan::extract(test_1_y2_g2)
a_mat<- data.frame(median=NA,l.95=NA,u.95=NA)
for(i in 1:23){
  a_mat[i,1]=median(plogis(ga_params$alpha+ga_params$a_yr[,i]))
  a_mat[i,2]=quantile(plogis(ga_params$alpha+ga_params$a_yr[,i]),0.975)
  a_mat[i,3]=quantile(plogis(ga_params$alpha+ga_params$a_yr[,i]),0.025)
}

comp_plot(ts=occ_ts(gray_angel),mod_mat = a_mat,sp='Gray Angelfish',GZ='Key West')

posterior <- as.array(test_1_y2_g)
dimnames(posterior)
color_scheme_set("viridis")
mcmc_areas(
  posterior,
  pars = c("alpha", paste('a_yr[',seq(1:23),']',sep=''),'sigma_hab','sigma_yr'),
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

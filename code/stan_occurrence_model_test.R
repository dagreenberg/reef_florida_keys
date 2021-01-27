logit_test<-"data{
  int<lower=1> N;//number of observations (SSU surveys)
  vector y[N]; //presence or absence on each survey
  int<lower=0> N_psu; //number of primary sample units
  int<lower=0> N_hab; //number of habitat classes
  int<lower=1,upper=N_psu> psu_id[N] // vector of psu identities
  int<lower=1,upper=h>    psu_in_hab[N_psu] // index of psu in each habitat class
  int<lower=1,upper=N_hab> hab_class[N] // vector of habitat class identities
  int<lower=1,upper=N_years> year[N] // vector of year identities

}
parameters {
  //global intercept
  real alpha;
  
  //deviations from intercept
  real dev_psu[N_psu]; //deviation between psus
  real dev_hab[N_hab]; //deviation between psus
  
  //st dev on the deviations
  real<lower = 0> sigma_psu;
  real<lower = 0> sigma_hab;   
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
  //priors
  alpha ~ normal(0,10);
  
  //standard deviations
  sigma_hab ~ cauchy(0, 2.5)
  sigma_psu ~ cauchy(0, 2.5)
  
  //varying intercepts
  dev_hab ~ normal(0, sigma_hab)
  dev_psu ~ normal(0, sigma_psu)
 
  for(i in 1:N){
    y[i] ~ bernoulli_logit(alpha_p[psu_id[i]]);
  }
}
"

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
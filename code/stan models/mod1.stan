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
  int<lower=1> N1;//number of observations (REEF surveys)
  int<lower=1> N2;//number of observations (REEF surveys)
  int y1[N1]; //abundance category for each survey
  int y2[N2]; //abundance category for each survey
  int<lower=0> N_psu; //number of habitat classes
  int<lower=1,upper=N_psu> psu_yr[N1]; // vector of RVC primary sample units
  int<lower=0> N_hab1; //number of habitat classes
  int<lower=1,upper=N_hab1> hab_class1[N1]; // vector of habitat class identities
  int<lower=0> N_strat1; //number of strata in RVC
  int<lower=1,upper=N_strat1> stratum1[N1]; // vector of RVC stratum identities
  int<lower=0> N_site; //number of sites in REEF
  int<lower=1,upper=N_site> site[N2]; // vector of site identities
  int<lower=0> N_hab2; //number of habitat classes - REEF
  int<lower=1,upper=N_hab2> hab_class2[N2]; // vector of habitat class identities
  int<lower=0> N_strat2; //number of strata - REEF
  int<lower=1,upper=N_strat2> stratum2[N2]; // vector of REEF stratum identities
  int<lower=0> N_dv; //number of divers
  int<lower=1,upper=N_dv> diver[N2]; // vector of diver identities
  int<lower=0> N_dmy; //number of site day clusters
  int<lower=1,upper=N_dmy> dmy[N2]; // vector of site day cluster identities
  int<lower=0> N_my; //number of monthly clusters
  int<lower=1,upper=N_my> my[N2]; // vector of monthly survey cluster identities
  int Z1; // columns in the covariate matrix
  int Z2; // columns in the covariate matrix
  matrix[N1,Z1] X1; // design matrix X for RVC
  matrix[N2,Z2] X2; // design matrix X for REEF
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
  real<lower=1e-9> recip_phi; //inverse of overdispersion parameter
  real a; //scalar for time-series 2

  //deviations from intercept
  vector[Z1] beta1; //effort coefficients - RVC
  vector[Z2] beta2; //effort coefficients - REEF
  vector[N_psu] a_psu; //deviation among habitats - REEF
  vector[N_hab1] a_hab1; //deviation between habitats - RVC
  vector[N_hab2] a_hab2; //deviation among habitats - REEF
  vector[N_strat1] a_strat1; //deviations among strata - RVC
  vector[N_strat2] a_strat2; //deviations among strata - REEF
  vector[N_site] a_site; //deviation between sites
  vector[N_dv] a_dv; //deviation between divers
  vector[N_dmy] a_dmy; //deviation between site day clusters
  vector[N_my] a_my; //deviation between site day clusters
 
  //variance on the deviance components
  real<lower = 0> sd_psu;
  real<lower = 0> sd_strat1;
  real<lower = 0> sd_strat2;
  real<lower = 0> sd_hab1;
  real<lower = 0> sd_hab2;
  real<lower = 0> sd_site;
  real<lower = 0> sd_dv;
  real<lower = 0> sd_dmy;
  real<lower = 0> sd_my;
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
  
  x[1] = x0 + pro_dev[1]*sd_q;
   
  for(t in 2:TT){
    x[t] = x[t-1] + pro_dev[t]*sd_q;
  }
   
  for(i in 1:N_yr1){
    a_yr1[i] = x[yr_index1[i]] + obs_dev1[i]*sd_r1; 
  }
    for(i in 1:N_yr2){
      a_yr2[i] = x[yr_index2[i]] + a + obs_dev2[i]*sd_r2; 
}
  
}  

model{
  //priors
  c ~ induced_dirichlet(rep_vector(1, K), 0); //prior on ordered cut-points
  beta1 ~ normal(0,2); //covariates - rvc
  beta2 ~ normal(0,2); //covariates - reef
  recip_phi ~ cauchy(0,5); //reciprocal of overdispersion
  x0 ~ normal(0,5); //initial state
  a ~ normal(0,5); //scalar
 
  //variance terms
  sd_psu ~ gamma(2,3);
  sd_strat1 ~ gamma(2,3);
  sd_hab1 ~ gamma(2,3);
  sd_strat2 ~ gamma(2,3);
  sd_hab2 ~ gamma(2,3);
  sd_q ~ gamma(2,10);
  sd_r1 ~ gamma(2,10);
  sd_r2 ~ gamma(2,10);
  sd_site ~ gamma(2,3);
  sd_dv ~ gamma(2,3);
  sd_dmy ~ gamma(2,3);
  sd_my ~ gamma(2,3);
  
  //varying intercepts
  a_psu ~ std_normal();
  a_strat1 ~ std_normal();
  a_hab1 ~ std_normal();
  a_hab2 ~ std_normal();
  a_strat2 ~ std_normal();
  a_site ~ std_normal();
  a_dv ~ std_normal();
  a_dmy ~ std_normal();
  a_my ~ std_normal();
  
  for(t in 1:N_yr1){
   obs_dev1[t] ~ std_normal();
  }
  
  for(t in 1:TT){
    pro_dev[t] ~ std_normal();
    obs_dev2[t] ~ std_normal();
  }

  y1 ~ neg_binomial_2_log(a_yr1[year_id1] + a_psu[psu_yr]*sd_psu + a_hab1[hab_class1]*sd_hab1 + a_strat1[stratum1]*sd_strat1 + X1*beta1,phi);
  
  y2 ~ ordered_logistic(a_yr2[year_id2]+a_site[site]*sd_site + a_hab2[hab_class2]*sd_hab2 + a_strat2[stratum2]*sd_strat2 +a_dv[diver]*sd_dv+a_dmy[dmy]*sd_dmy+a_my[my]*sd_my+X2*beta2,c);
  
}
  generated quantities{
  vector[N1+N2] log_lik;
  vector[N1] rvc_y_rep;
  vector[N2] reef_y_rep;
  for (i in 1:N1){ log_lik[i] = neg_binomial_2_log_lpmf(y1[i]|a_yr1[year_id1[i]] + a_psu[psu_yr[i]]*sd_psu + a_hab1[hab_class1[i]]*sd_hab1 + a_strat1[stratum1[i]]*sd_strat1 + X1[i,]*beta1,phi);
  rvc_y_rep[i] = neg_binomial_2_log_rng(a_yr1[year_id1[i]] + a_psu[psu_yr[i]]*sd_psu+ a_hab1[hab_class1[i]]*sd_hab1 + a_strat1[stratum1[i]]*sd_strat1 + X1[i,]*beta1,phi);
}
for (z in 1:N2){ log_lik[N1+z] = ordered_logistic_lpmf(y2[z]|a_yr2[year_id2[z]]+a_site[site[z]]*sd_site+a_hab2[hab_class2[z]]*sd_hab2+ a_strat2[stratum2[z]]*sd_strat2 +a_dv[diver[z]]*sd_dv+a_dmy[dmy[z]]*sd_dmy+a_my[my[z]]*sd_my+X2[z,]*beta2, c);
 reef_y_rep[z] = ordered_logistic_rng(a_yr2[year_id2[z]]+a_site[site[z]]*sd_site+a_hab2[hab_class2[z]]*sd_hab2+a_dv[diver[z]]*sd_dv+ a_strat2[stratum2[z]]*sd_strat2+a_dmy[dmy[z]]*sd_dmy+a_my[my[z]]*sd_my+X2[z,]*beta2, c);
}
} 

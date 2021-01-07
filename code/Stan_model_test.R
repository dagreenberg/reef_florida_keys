library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())




mod_1<-"data{
  int<lower=0> TT; //Timespan
  int<lower=0> N; // num of ts; rows of y
  int<lower=0> n_pos; // number of non-NA values
  int<lower=0> col_indx_pos[n_pos]; // col index of non-NA vals
  int<lower=0> row_indx_pos[n_pos]; // row index of non-NA vals
  int<lower=0> ts_id[n_pos]; //id for the scaling parameter (a)
  vector[n_pos] y;
}
parameters{
  real x0;
  real a;
  vector[TT] pro_dev;
  real<lower=0> sd_q;
  real<lower=0> sd_r[N];
}
transformed parameters{
  vector[TT] x;
  x[1] = x0 + pro_dev[1];
  for(t in 2:TT) {
    x[t] = x[t-1] + pro_dev[t];
  }
}
model{
  x0 ~ normal(y[1],10);
  a ~ normal(0,10);
  sd_q ~ cauchy(0,1);
  
  for(i in 1:N){
    sd_r[i] ~ cauchy(0,1);
  }
  
  for(t in 1:TT){
    pro_dev[t] ~ normal(0, sd_q);
  }
  
  
  for(i in 1:n_pos){
    y[i] ~ normal(x[col_indx_pos[i]]+a*ts_id[i], sd_r[row_indx_pos[i]]);
  }
}
generated quantities{
  vector[n_pos] log_lik;
  for (i in 1:n_pos) log_lik[i] = normal_lpdf(y[i]|x[col_indx_pos[i]]+a*ts_id[i], sd_r);
}
"
mod_1.1<-"data{
  int<lower=0> TT; //Timespan
  int<lower=0> N; // num of ts; rows of y
  int<lower=0> n_pos; // number of non-NA values
  int<lower=0> col_indx_pos[n_pos]; // col index of non-NA vals
  int<lower=0> row_indx_pos[n_pos]; // row index of non-NA vals
  int<lower=0> ts_id[n_pos]; //id for the scaling parameter (a)
  vector[n_pos] y;
}
parameters{
  real x0;
  real a;
  real u;
  vector[TT] pro_dev;
  real<lower=0> sd_q;
  real<lower=0> sd_r[N];
}
transformed parameters{
  vector[TT] x;
  x[1] = x0 + u + pro_dev[1];
  for(t in 2:TT) {
    x[t] = x[t-1] + u + pro_dev[t];
  }
}
model{
  x0 ~ normal(y[1],10);
  a ~ normal(0,10);
  sd_q ~ cauchy(0,1);
  u ~ normal(0,2);
  for(i in 1:N){
    sd_r[i] ~ cauchy(0,1);
  }
  
  for(t in 1:TT){
    pro_dev[t] ~ normal(0, sd_q);
  }
  
  
  for(i in 1:n_pos){
    y[i] ~ normal(x[col_indx_pos[i]]+a*ts_id[i], sd_r[row_indx_pos[i]]);
  }
}
generated quantities{
  vector[n_pos] log_lik;
  for (i in 1:n_pos) log_lik[i] = normal_lpdf(y[i]|x[col_indx_pos[i]]+a*ts_id[i], sd_r);
}
"


##
mod<- 
"data {
  int<lower=0> TT; // length of ts
  int<lower=0> N; // num of ts; rows of y
  int<lower=0> n_pos; // number of non-NA values in y
  int<lower=0> col_indx_pos[n_pos]; // col index of non-NA vals
  int<lower=0> row_indx_pos[n_pos]; // row index of non-NA vals
  int<lower=0> ts_pos[n_pos] // binary for rvc and reef timeseries
  vector[n_pos] y;
}
parameters {
  vector x0; // initial state
  real u;
  real a;
  vector pro_dev[TT]; // refed as pro_dev[TT]
  real<lower=0> sd_q;
  real<lower=0> sd_r[N]; // obs variances are different
}
transformed parameters {
  vector[TT] x; // refed as x[TT,N]
    x[1] = x0 + a*ts_pos + pro_dev[1];
    for(t in 2:TT){
      x[t] = x[t-1] + a*ts_pos + pro_dev[t];
    }
  }
}
model {
  sd_q ~ cauchy(0,5);
 )
  for(i in 1:N){
    x0[i] ~ normal(y[i],10); // assume no missing y[1]
    sd_r[i] ~ cauchy(0,5);
    a ~ normal(0,10
    for(t in 1:TT){
    pro_dev[t,i] ~ normal(0, sd_q);
    }
  }
  u ~ normal(0,2);
  for(i in 1:n_pos){
    y[i] ~ normal(x[col_indx_pos[i], row_indx_pos[i]]+a, sd_r[row_indx_pos[i]]);
  }
}
generated quantities {
  vector[n_pos] log_lik;
  for (n in 1:n_pos) log_lik[n] = normal_lpdf(y[n] | x[col_indx_pos[n], row_indx_pos[n]], sd_r[row_indx_pos[n]]);
}
"

mod_2<-"
data {
  int<lower=0> TT; // length of ts
  int<lower=0> N; // num of ts; rows of y
  int<lower=0> n_pos; // number of non-NA values in y
  int<lower=0> col_indx_pos[n_pos]; // col index of non-NA vals
  int<lower=0> row_indx_pos[n_pos]; // row index of non-NA vals
  vector[n_pos] y;
}
parameters {
  vector[N] x0; // initial states
  vector[N] pro_dev[TT]; // refed as pro_dev[TT,N]
  real<lower=0> sd_q;
  real<lower=0> sd_r[N]; // obs variances are different
}
transformed parameters {
  vector[N] x[TT]; // refed as x[TT,N]
  for(i in 1:N){
    x[1,i] = x0[i] + pro_dev[1,i];
    for(t in 2:TT) {
      x[t,i] = x[t-1,i] + pro_dev[t,i];
    }
  }
}
model {
  sd_q ~ cauchy(0,1);
  for(i in 1:N){
    x0[i] ~ normal(y[i],10); // assume no missing y[1]
    sd_r[i] ~ cauchy(0,1);
    for(t in 1:TT){
    pro_dev[t,i] ~ normal(0, sd_q);
    }
  }
 
  for(i in 1:n_pos){
    y[i] ~ normal(x[col_indx_pos[i], row_indx_pos[i]], sd_r[row_indx_pos[i]]);
  }
}
generated quantities {
  vector[n_pos] log_lik;
  for (n in 1:n_pos) log_lik[n] = normal_lpdf(y[n] | x[col_indx_pos[n], row_indx_pos[n]], sd_r[row_indx_pos[n]]);
}
"

mod_2.1<-"
data {
  int<lower=0> TT; // length of ts
  int<lower=0> N; // num of ts to compare; first 2 rows of y
  int<lower=0> J; // num of ts for regularizing; rows of y
  int<lower=0> n_pos; // number of non-NA values in y
  int<lower=0> col_indx_pos[n_pos]; // col index of non-NA vals
  int<lower=0> row_indx_pos[n_pos]; // row index of non-NA vals
  vector[n_pos] y;
}
parameters {
  vector[N] x0; // initial states
  vector[N] pro_dev[TT]; // refed as pro_dev[TT,N]
  real<lower=0> sd_q; //shared proc variance
  real<lower=0> sd_r[N]; // obs variances are different
  vector[N] u; //trend
}
transformed parameters {
  vector[N] x[TT]; // refed as x[TT,N]
  for(i in 1:N){
    x[1,i] = x0[i] + u[i] + pro_dev[1,i];
    for(t in 2:TT) {
      x[t,i] = x[t-1,i] + u[i] + pro_dev[t,i];
    }
  }
}
model {
  sd_q ~ cauchy(0,1);
  for(i in 1:N){
    x0[i] ~ normal(y[i],10); // assume no missing y[1]
    sd_r[i] ~ cauchy(0,5);
    u[i] ~ normal(0,2);
    for(t in 1:TT){
    pro_dev[t,i] ~ normal(0, sd_q);
    }
  }
 
  for(i in 1:n_pos){
    y[i] ~ normal(x[col_indx_pos[i], row_indx_pos[i]], sd_r[row_indx_pos[i]]);
  }
}
generated quantities {
  vector[n_pos] log_lik;
  for (n in 1:n_pos) log_lik[n] = normal_lpdf(y[n] | x[col_indx_pos[n], row_indx_pos[n]], sd_r[row_indx_pos[n]]);
}
"

mod_3<-"
data {
  int<lower=0> TT; // length of ts
  int<lower=0> N; // num of ts; rows of y
  int<lower=0> n_pos; // number of non-NA values in y
  int<lower=0> col_indx_pos[n_pos]; // col index of non-NA vals
  int<lower=0> row_indx_pos[n_pos]; // row index of non-NA vals
  vector[n_pos] y;
}
parameters {
  vector[N] x0; // initial states
  vector[N] pro_dev[TT]; // refed as pro_dev[TT,N]
  vector[N] u;
  real<lower=0> sd_q; // process variances are equal
  real<lower=0> sd_r[N]; // obs variances are different
}
transformed parameters {
  vector[N] x[TT]; // refed as x[TT,N]
  for(i in 1:N){
    x[1,i] = x0[i] + u[i] + pro_dev[1,i];
    for(t in 2:TT) {
      x[t,i] = x[t-1,i] + u[i] + pro_dev[t,i];
    }
  }
}
model {
  sd_q ~ cauchy(0,5);
  for(i in 1:N){
    x0[i] ~ normal(y[i],10); // assume no missing y[1]
    sd_r[i] ~ cauchy(0,5);
    u[i] ~ normal(0,2);
    for(t in 1:TT){
    pro_dev[t,i] ~ normal(0, sd_q);
   
    }
  }
 
  for(i in 1:n_pos){
    y[i] ~ normal(x[col_indx_pos[i], row_indx_pos[i]], sd_r[row_indx_pos[i]]);
  }
}
generated quantities {
  vector[n_pos] log_lik;
  for (n in 1:n_pos) log_lik[n] = normal_lpdf(y[n] | x[col_indx_pos[n], row_indx_pos[n]], sd_r[row_indx_pos[n]]);
}
"

Y<- ts_comp[[1]]

ypos <- Y[!is.na(Y)]
n_pos <- length(ypos)  #number on non-NA ys
indx_pos <- which(!is.na(Y), arr.ind = TRUE)  #index on the non-NAs
col_indx_pos <- as.vector(indx_pos[, "col"])
row_indx_pos <- as.vector(indx_pos[, "row"])
ts_pos<- ifelse(row_indx_pos==1,0,1)

test_1 <- rstan::stan(model_code = mod_1, data = list(y = ypos, 
                                                   TT = ncol(Y), N = nrow(Y), n_pos = n_pos, ts_id=ts_pos, col_indx_pos = col_indx_pos, 
                                                   row_indx_pos = row_indx_pos), pars = c("sd_q", "x", "sd_r",'a',
                                                                                          "x0"),control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 1000, chains = 4, iter = 2000, thin = 1)
test_1.1 <- rstan::stan(model_code = mod_1.1, data = list(y = ypos, 
                                                      TT = ncol(Y), N = nrow(Y), n_pos = n_pos, ts_id=ts_pos, col_indx_pos = col_indx_pos, 
                                                      row_indx_pos = row_indx_pos), pars = c("sd_q", "x", "sd_r",'a',"u",
                                                                                             "x0",'log_lik'),control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 1000, chains = 4, iter = 2000, thin = 1)


test_2 <- rstan::stan(model_code = mod_2, data = list(y = ypos, 
                                                      TT = ncol(Y), N = nrow(Y), n_pos = n_pos, col_indx_pos = col_indx_pos, 
                                                      row_indx_pos = row_indx_pos), pars = c("sd_q", "x", "sd_r", 
                                                                                             "x0"),control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 1000, chains = 4, iter = 2000, thin = 1)
test_2.1 <- rstan::stan(model_code = mod_2.1, data = list(y = ypos, 
                                                      TT = ncol(Y), N = nrow(Y), n_pos = n_pos, col_indx_pos = col_indx_pos, 
                                                      row_indx_pos = row_indx_pos), pars = c("sd_q", "x", "sd_r", "u", 
                                                                                             "x0"),control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 1000, chains = 4, iter = 2000, thin = 1)

test_3<- rstan::stan(model_code = mod_3, data = list(y = ypos, 
                                                      TT = ncol(Y), N = nrow(Y), n_pos = n_pos, col_indx_pos = col_indx_pos, 
                                                      row_indx_pos = row_indx_pos), pars = c("sd_q", "x", "sd_r", 
                                                                                             "u", "x0","log_lik"),control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 1000, chains = 4, iter = 2000, thin = 1)
loo(test_3)


library(shinystan)
launch_shinystan(test_2)

m2<- extract(test_2)
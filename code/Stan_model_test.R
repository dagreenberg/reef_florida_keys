library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

mod_1<-"data{
  int<lower=0> TT; //Timespan
  int<lower=0> N; // num of ts; rows of y
  int<lower=0> n_pos; // number of non-NA values
  int<lower=0> indx_pos[n_pos]; // index of the non-NA value
  int<lower=0> id; //id for the scaling parameter (a)
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
  x[1] = x0 + a*id + pro_dev[1];
  for(i in 2:TT) {
    x[i] = x[i-1] + a*id + pro_dev[i];
  }
}
model{
  x0 ~ normal(y[1],10);
  a ~ normal(0,10);
  sd_q ~ cauchy(0,5);
  sd_r ~ cauchy(0,5);
  pro_dev ~ normal(0, sd_q);
  
  for(i in 1:n_pos){
    y[i] ~ normal(x[indx_pos[i]], sd_r);
  }
}
generated quantities{
  vector[n_pos] log_lik;
  for (i in 1:n_pos) log_lik[i] = normal_lpdf(y[i] | x[indx_pos[i]], sd_r);
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
  vector[n_pos] y;
}
parameters {
  vector x0; // initial state
  real u;
  real a;
  vector[N] pro_dev[TT]; // refed as pro_dev[TT,N]
  real<lower=0> sd_q;
  real<lower=0> sd_r[N]; // obs variances are different
}
transformed parameters {
  vector[TT] x; // refed as x[TT,N]
  for(i in 1:N){
    x[1,i] = x0[i] + u + pro_dev[1,i];
    for(t in 2:TT) {
      x[t,i] = x[t-1,i] + u + pro_dev[t,i];
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
  sd_q ~ cauchy(0,5);
  for(i in 1:N){
    x0[i] ~ normal(y[i],10); // assume no missing y[1]
    sd_r[i] ~ cauchy(0,5);
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
    for(t in 1:TT){
    pro_dev[t,i] ~ normal(0, sd_q);
    u[i] ~ normal(0,2);
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

ypos <- Y[!is.na(Y)]
n_pos <- length(ypos)  #number on non-NA ys
indx_pos <- which(!is.na(Y), arr.ind = TRUE)  #index on the non-NAs
col_indx_pos <- as.vector(indx_pos[, "col"])
row_indx_pos <- as.vector(indx_pos[, "row"])


test_1 <- rstan::stan(model_code = mod_1, data = list(y = ypos, 
                                                   TT = ncol(Y), N = nrow(Y), n_pos = n_pos, col_indx_pos = col_indx_pos, 
                                                   row_indx_pos = row_indx_pos), pars = c("sd_q", "x", "sd_r", 
                                                                                          "u", "x0"),control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 1000, chains = 4, iter = 2000, thin = 1)
test_2 <- rstan::stan(model_code = mod_2, data = list(y = ypos, 
                                                      TT = ncol(Y), N = nrow(Y), n_pos = n_pos, col_indx_pos = col_indx_pos, 
                                                      row_indx_pos = row_indx_pos), pars = c("sd_q", "x", "sd_r", 
                                                                                            "x0"),control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 1000, chains = 4, iter = 2000, thin = 1)

test_3<- rstan::stan(model_code = mod_3, data = list(y = ypos, 
                                                      TT = ncol(Y), N = nrow(Y), n_pos = n_pos, col_indx_pos = col_indx_pos, 
                                                      row_indx_pos = row_indx_pos), pars = c("sd_q", "x", "sd_r", 
                                                                                             "u", "x0"),control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 1000, chains = 4, iter = 2000, thin = 1)



library(shinystan)
launch_shinystan(test_2)

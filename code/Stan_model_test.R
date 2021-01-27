rm(list=ls())

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#####
TS_plot_stan<- function(ts,sp,GZ,mod){
  #pdf(paste(paste(sp,GZ,sep='_'),'.pdf',sep=''),width=8,height=6)
  plot(ts[1,]~c(seq(1993,2018)),type='n',ylim=c(min(na.omit(c(ts))),max(na.omit(c(ts)))),col='darkblue',bty='l',ylab=expression('log'[10]*' (Counts per survey)'),xlab='Year',main=paste(sp,GZ,sep=', '))
  lines(ts[1,]~c(seq(1993,2018)),col='navy',lwd=2)
  points(ts[1,]~c(seq(1993,2018)),col='white',pch=21,bg='navy',cex=1.2)
  lines(ts[2,]~c(seq(1993,2018)),col='darkred',lwd=2)
  points(ts[2,]~c(seq(1993,2018)),col='white',pch=21,bg='darkred',cex=1.2)
  if(mod==1){
    x=NA
    for(z in 1:ncol(params_1$x)){
      x[z]=median(params_1$x[,z])
    }
    lines(x~c(seq(1993,2018)),lty=5,lwd=2,col='slategray')
    lines(x+median(params_1$a)~c(seq(1993,2018)),lty=5,lwd=2,col='slategray')
  }
  if(mod==2){
    x1=NA
    x2=NA
    for(z in 1:ncol(params_2$x[,,1])){
      x1[z]=median(params_2$x[,,1][,z])
      x2[z]=median(params_2$x[,,2][,z])
    }
    
    lines(x1~c(seq(1993,2018)),lty=5,lwd=2,col='steelblue')
    lines(x2~c(seq(1993,2018)),lty=5,lwd=2,col='firebrick')
  }
  if(mod==3){
    x1=NA
    x2=NA
    for(z in 1:ncol(params_3$x[,,1])){
      x1[z]=median(params_3$x[,,1][,z])
      x2[z]=median(params_3$x[,,2][,z])
    }
    
    lines(x1~c(seq(1993,2018)),lty=5,lwd=2,col='steelblue')
    lines(x2~c(seq(1993,2018)),lty=5,lwd=2,col='firebrick')
  }
  
  legend(2016,c(max(na.omit(c(ts)))*1.05),c('RVC','REEF'),text.col=c('navy','darkred'),bty='n')
 # dev.off(paste(paste(sp,GZ,sep='_'),'.pdf',sep=''))
}
####




load("C:/Users/14388/Desktop/reef_florida_keys_data/stan_model_testing.RData")


mod_1<-"data{
  int<lower=0> TT; //Timespan
  int<lower=0> N; // num of ts; rows of y
  int<lower=0> n_pos; // number of non-NA values
  int<lower=0> col_indx_pos[n_pos]; // col index of non-NA vals
  int<lower=0> row_indx_pos[n_pos]; // row index of non-NA vals
  int<lower=0> ts_id[n_pos]; // id for the scaling parameter (a)
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
  x0 ~ normal(0,5);
  a ~ normal(0,5);
  sd_q ~ cauchy(0,1);
  
  for(i in 1:N){
    sd_r[i] ~ cauchy(0,1);
  }
  
  for(t in 1:TT){
    pro_dev[t] ~ normal(0, sd_q);
  }
  
  
  for(i in 1:n_pos){
    y[i] ~ normal(x[col_indx_pos[i]]+(a*ts_id[i]), sd_r[row_indx_pos[i]]);
  }
}
generated quantities{
  vector[n_pos] log_lik;
  for (i in 1:n_pos) log_lik[i] = normal_lpdf(y[i]|x[col_indx_pos[i]]+(a*ts_id[i]), sd_r[row_indx_pos[i]]);
}
"

mod_1.2<-"data{
  int<lower=0> TT; //Timespan
  int<lower=0> N; // num of ts; rows of y
  int<lower=0> n_pos; // number of non-NA values
  int<lower=0> col_indx_pos[n_pos]; // col index of non-NA vals
  int<lower=0> row_indx_pos[n_pos]; // row index of non-NA vals
  int<lower=0> ts_id[n_pos]; // id for the scaling parameter (a)
  vector[n_pos] y;
}
parameters{
  real x0;
  real a;
  vector[TT] x;
  real<lower=0> sd_q;
  real<lower=0> sd_r[N];
}
model{
  x0 ~ normal(y[1],10);
  a ~ normal(0,10);
  sd_q ~ cauchy(0,.1);
  
  for(i in 1:N){
    sd_r[i] ~ cauchy(0,.1);
  }
  
   for(t in 2:TT) {
      x[t] ~ normal(x[t-1],sd_q);
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
mod_3.2<-"
data {
  int<lower=0> TT; // length of ts
  int<lower=0> N; // num of ts; rows of y
  int<lower=0> n_pos; // number of non-NA values in y
  int<lower=0> col_indx_pos[n_pos]; // col index of non-NA vals
  int<lower=0> row_indx_pos[n_pos]; // row index of non-NA vals
  vector[n_pos] y;
}
parameters {
  vector[N] x[TT]; // refed as x[TT,N];
  vector[N] x0; // initial states
  real u[N]; //trend
  real<lower=0> sd_q[N];
  real<lower=0> sd_r[N]; // obs variances are different
}

model{
sd_q ~ cauchy(0,1);
  for(i in 1:N){
    x0[i] ~ normal(y[i],10); // assume no missing y[1]
    sd_r[i] ~ cauchy(0,1);
    
    x[1,i] ~ normal(x0[i],sd_q[i]);
    for(t in 2:TT) {
      x[t,i] ~ normal(x[t-1,i]+u[i],sd_q[i]);
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

mod_2_nm_scalar<-"
data {
  int<lower=0> TT; // length of ts
  int<lower=0> N; // num of ts; rows of y
  int<lower=0> n_pos; // number of non-NA values in y
  int<lower=0> col_indx_pos[n_pos]; // col index of non-NA vals
  int<lower=0> row_indx_pos[n_pos]; // row index of non-NA vals
  vector[n_pos] nm;
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
    y[i] ~ normal(x[col_indx_pos[i], row_indx_pos[i]], sd_r[row_indx_pos[i]]/nm[i]);
  }
}
generated quantities {
  vector[n_pos] log_lik;
  for (n in 1:n_pos) log_lik[n] = normal_lpdf(y[n] | x[col_indx_pos[n], row_indx_pos[n]], sd_r[row_indx_pos[n]]/nm[n]);
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

mod_3.1<-"
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
  real<lower=0> sd_q[N];
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
  for(i in 1:N){
    x0[i] ~ normal(y[i],10); // assume no missing y[1]
    sd_r[i] ~ cauchy(0,1);
    sd_q[i] ~ cauchy(0,1);
    for(t in 1:TT){
    pro_dev[t,i] ~ normal(0, sd_q[i]);
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

mod_3.2<-"
data {
  int<lower=0> TT; // length of ts
  int<lower=0> N; // num of ts; rows of y
  int<lower=0> n_pos; // number of non-NA values in y
  int<lower=0> col_indx_pos[n_pos]; // col index of non-NA vals
  int<lower=0> row_indx_pos[n_pos]; // row index of non-NA vals
  vector[n_pos] y;
}
parameters {
  vector[N] x[TT]; // refed as x[TT,N];
  vector[N] x0; // initial states
  real u[N]; //trend
  real<lower=0> sd_q[N];
  real<lower=0> sd_r[N]; // obs variances are different
}

model{
  for(i in 1:N){
    x0[i] ~ normal(y[i],10); // assume no missing y[1]
    sd_r[i] ~ cauchy(0,1);
    sd_q[i] ~ cauchy(0,1);
    u[i] ~ normal(0,2);
    x[1,i] ~ normal(x0[i],sd_q[i]);
    for(t in 2:TT) {
      x[t,i] ~ normal(x[t-1,i]+u[i],sd_q[i]);
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



s0[1] ~ normal(init_s0,sigma_surv_pro);
for(i in 2:N)
{
  s0[i] ~ normal(s0[i-1],sigma_surv_pro);
}



Y<- ts_comp[[1]]

ypos <- Y[!is.na(Y)]
n_pos <- length(ypos)  #number on non-NA ys
indx_pos <- which(!is.na(Y), arr.ind = TRUE)  #index on the non-NAs
col_indx_pos <- as.vector(indx_pos[, "col"])
row_indx_pos <- as.vector(indx_pos[, "row"])
ts_pos<- ifelse(row_indx_pos==1,0,1)

test_1<- rstan::stan(model_code = mod_1, data = list(y = ypos, 
                                                   TT = ncol(Y), N = nrow(Y), n_pos = n_pos, ts_id=ts_pos, col_indx_pos = col_indx_pos, 
                                                   row_indx_pos = row_indx_pos), pars = c("sd_q", "x", "sd_r",'a',
                                                                                         "x0"),control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 1000, chains = 4, iter = 2000, thin = 1)
test_1.1 <- rstan::stan(model_code = mod_1.1, data = list(y = ypos, 
                                                      TT = ncol(Y), N = nrow(Y), n_pos = n_pos, ts_id=ts_pos, col_indx_pos = col_indx_pos, 
                                                      row_indx_pos = row_indx_pos), pars = c("sd_q", "x", "sd_r",'a',
                                                                                             "x0",'u'),control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 1000, chains = 4, iter = 2000, thin = 1)


test_2 <- rstan::stan(model_code = mod_2, 
                      data = list(y = ypos, TT = ncol(Y), N = nrow(Y), n_pos = n_pos, col_indx_pos = col_indx_pos, 
                      row_indx_pos = row_indx_pos), pars = c("sd_q", "x", "sd_r","x0",'log_lik'),
                      control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 1000, chains = 4, iter = 2000, thin = 1)

test_3.1 <- rstan::stan(model_code = mod_3.1, 
                        data = list(y = ypos, TT = ncol(Y), N = nrow(Y), n_pos = n_pos, col_indx_pos = col_indx_pos, 
                                    row_indx_pos = row_indx_pos), pars = c("sd_q", "x", "sd_r","x0",'u'),
                        control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 1000, chains = 4, iter = 2000, thin = 1)



test_3.2 <- rstan::stan(model_code = mod_3.2, 
                      data = list(y = ypos, TT = ncol(Y), N = nrow(Y), n_pos = n_pos, col_indx_pos = col_indx_pos, 
                                  row_indx_pos = row_indx_pos), pars = c("sd_q", "x", "sd_r","x0",'u'),
                      control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 1000, chains = 4, iter = 2000, thin = 1)

                                                                                             
                                                                                             
                                                                                                                                                                                   
test_2.1 <- rstan::stan(model_code = mod_2.1, data = list(y = ypos, 
                                                      TT = ncol(Y), N = nrow(Y), n_pos = n_pos, col_indx_pos = col_indx_pos, 
                                                      row_indx_pos = row_indx_pos), pars = c("sd_q", "x", "sd_r", "u", 
                                                                                             "x0"),control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 1000, chains = 4, iter = 2000, thin = 1)

test_3<- rstan::stan(model_code = mod_3, data = list(y = ypos, 
                                                      TT = ncol(Y), N = nrow(Y), n_pos = n_pos, col_indx_pos = col_indx_pos, 
                                                      row_indx_pos = row_indx_pos), pars = c("sd_q", "x", "sd_r", 
                                                                                             "x0"),control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 1000, chains = 4, iter = 2000, thin = 1)


#library(shinystan)
shinystan::launch_shinystan(test_1.1)

#m2<- extract(test_2)



mars_3403_stan<- data.frame(SP=NA,loo1=NA,loo1.se=NA,loo2=NA,loo2.se=NA,loo3=NA,loo3.se=NA,
                            Q1=NA,Q2=NA,Q3.rvc=NA,Q3.reef=NA,R1.rvc=NA,R1.reef=NA,R2.rvc=NA,R2.reef=NA,R3.rvc=NA,R3.reef=NA,
                            best.mod=NA,waic1=NA,p_waic1=NA,waic2=NA,p_waic2=NA,waic3=NA,p_waic3=NA)

setwd("C:/Users/14388/Desktop/Scripps - Project 1 RVC and REEF/RVC/timeseries/sg stan test")
for(i in 1:length(rvc.green.3403.sp)){
  spp<- dplyr::filter(fish_reef,commonname==rvc.green.3403.sp[i])
  ts_comp[[i]]<- t(log10(rvc_3403_green[[i]]$mean_ssu_abundance)) #extract logged time-series
  m<- match(spp$commonname,REEF.green.3403.sp)
  if(is.na(m)==T){next}
  ts_comp[[i]]<- rbind(ts_comp[[i]],t(log10(REEF_3403_green[[m]]$mean_site_abund[1:26])))
  rownames(ts_comp[[i]])<- c(paste(spp$commonname,'rvc',sep=":"),paste(spp$commonname,'REEF',sep=":"))
  colnames(ts_comp[[i]])<- years
  
  Y<- ts_comp[[i]]
  
  ypos <- Y[!is.na(Y)]
  n_pos <- length(ypos)  #number on non-NA ys
  indx_pos <- which(!is.na(Y), arr.ind = TRUE)  #index on the non-NAs
  col_indx_pos <- as.vector(indx_pos[, "col"])
  row_indx_pos <- as.vector(indx_pos[, "row"])
  
  ts_pos<- ifelse(row_indx_pos==1,0,1)
  
  test_1<- rstan::stan(model_code = mod_1, data = list(y = ypos, 
                                                       TT = ncol(Y), N = nrow(Y), n_pos = n_pos, ts_id=ts_pos, col_indx_pos = col_indx_pos, 
                                                       row_indx_pos = row_indx_pos), pars = c("sd_q", "x", "sd_r",'a',
                                                                                              "x0",'log_lik'),control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 1000, chains = 4, iter = 2000, thin = 1)
  test_2 <- rstan::stan(model_code = mod_2, 
                        data = list(y = ypos, TT = ncol(Y), N = nrow(Y), n_pos = n_pos, col_indx_pos = col_indx_pos, 
                                    row_indx_pos = row_indx_pos), pars = c("sd_q", "x", "sd_r","x0",'log_lik'),
                        control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 1000, chains = 4, iter = 2000, thin = 1)
  
  test_3 <- rstan::stan(model_code = mod_3.1, 
                        data = list(y = ypos, TT = ncol(Y), N = nrow(Y), n_pos = n_pos, col_indx_pos = col_indx_pos, 
                                    row_indx_pos = row_indx_pos), pars = c("sd_q", "x", "sd_r","x0",'log_lik'),
                        control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 1000, chains = 4, iter = 2000, thin = 1)
  
  params_1<- extract(test_1)
  params_2<- extract(test_2)
  params_3<- extract(test_3)
  
 loo1= loo(test_1)
  loo2= loo(test_2)
  loo3=loo(test_3)
  
  ll1<-  loo::extract_log_lik(test_1)
  ll2<-  loo::extract_log_lik(test_2)
  ll3<-  loo::extract_log_lik(test_3)
  waic1<- loo::waic(ll1)
  waic2<- loo::waic(ll2)
  waic3<- loo::waic(ll3)
  
#  waic_comp<- loo:compare(waic1,waic2,waic3)
  
  loo_comp<- loo::loo_compare(loo1,loo2,loo3)
  
  
  mars_3403_stan[i,1]=spp$commonname
  mars_3403_stan[i,2]=loo_comp["model1",1]
  mars_3403_stan[i,3]=loo_comp["model1",2]  
  mars_3403_stan[i,4]=loo_comp["model2",1]
  mars_3403_stan[i,5]=loo_comp["model2",2]
  mars_3403_stan[i,6]=loo_comp["model3",1]
  mars_3403_stan[i,7]=loo_comp["model3",2]
  mars_3403_stan[i,8]=median(params_1$sd_q) 
#  mars_3403_stan[i,9]=sort(params_1$sd_q)[length(params_1$sd_q)*0.975]
  mars_3403_stan[i,9]=median(params_2$sd_q) 
 # mars_3403_stan[i,11]=sort(params_2$sd_q)[length(params_1$sd_q)*0.975]
  mars_3403_stan[i,10]=median(params_3$sd_q[,1]) 
#  mars_3403_stan[i,13]=sort(params_3$sd_q[,1])[length(params_3$sd_q[,1])*0.975]
  mars_3403_stan[i,11]=median(params_3$sd_q[,2]) 
 # mars_3403_stan[i,15]=sort(params_3$sd_q[,2])[length(params_3$sd_q[,2])*0.975]
  mars_3403_stan[i,12]=median(params_1$sd_r[,1]) 
  mars_3403_stan[i,13]=median(params_1$sd_r[,2]) 
  mars_3403_stan[i,14]=median(params_2$sd_r[,1]) 
  mars_3403_stan[i,15]=median(params_2$sd_r[,2]) 
  mars_3403_stan[i,16]=median(params_3$sd_r[,1]) 
  mars_3403_stan[i,17]=median(params_3$sd_r[,2]) 
  if(mars_3403_stan[i,2]==0){mars_3403_stan[i,18]=1}
  if(mars_3403_stan[i,4]==0){mars_3403_stan[i,18]=2}
  if(mars_3403_stan[i,6]==0){mars_3403_stan[i,18]=3}
  mars_3403_stan[i,19]=waic1$elpd_waic
  mars_3403_stan[i,20]=waic1$p_waic
  mars_3403_stan[i,21]=waic2$elpd_waic
  mars_3403_stan[i,22]=waic2$p_waic
  mars_3403_stan[i,23]=waic3$elpd_waic
  mars_3403_stan[i,24]=waic3$p_waic
  

  TS_plot_stan(ts=ts_comp[[i]],sp=spp$commonname,GZ='Key Largo',mod=3)
  dev.off()
  print(i)
}

TS_plot<- function(ts,sp,GZ){
  plot(ts[1,]~c(seq(1993,2018)),type='n',ylim=c(min(na.omit(c(ts[[i]]))),max(na.omit(c(ts[[i]])))),col='darkblue',bty='l',ylab=expression('log'[10]*' (Counts per survey)'),xlab='Year',main=paste(sp,GZ,sep=', '))
  lines(ts[1,]~c(seq(1993,2018)),col='navy')
  points(ts[1,]~c(seq(1993,2018)),col='white',pch=21,bg='navy',cex=1.2)
  lines(ts[2,]~c(seq(1993,2018)),col='darkred')
  points(ts[2,]~c(seq(1993,2018)),col='white',pch=21,bg='darkred',cex=1.2)
  legend(2016,c(max(na.omit(c(ts[[i]])))*1.05),c('RVC','REEF'),text.col=c('navy','darkred'),bty='n')
  
}

TS_plot.pdf<- function(ts,sp,GZ){
  pdf(paste(paste(sp,GZ,sep='_'),'.pdf',sep=''),width=8,height=6)
  plot(ts[1,]~c(seq(1993,2018)),type='n',ylim=c(min(na.omit(c(ts))),max(na.omit(c(ts)))),col='darkblue',bty='l',ylab=expression('log'[10]*' (Counts per survey)'),xlab='Year',main=paste(sp,GZ,sep=', '))
  lines(ts[1,]~c(seq(1993,2018)),col='navy')
  points(ts[1,]~c(seq(1993,2018)),col='white',pch=21,bg='navy',cex=1.2)
  lines(ts[2,]~c(seq(1993,2018)),col='darkred')
  points(ts[2,]~c(seq(1993,2018)),col='white',pch=21,bg='darkred',cex=1.2)
  legend(2016,c(max(na.omit(c(ts)))*1.05),c('RVC','REEF'),text.col=c('navy','darkred'),bty='n')
  dev.off(paste(paste(sp,GZ,sep='_'),'.pdf',sep=''))
}

TS__stan_plot.pdf<- function(ts,sp,GZ){
  ##extract hidden state
  x_mod1<- extract(test_1)$x
  
  
  pdf(paste(paste(sp,GZ,sep='_'),'.pdf',sep=''),width=8,height=6)
  plot(ts[1,]~c(seq(1993,2018)),type='n',ylim=c(min(na.omit(c(ts))),max(na.omit(c(ts)))),col='darkblue',bty='l',ylab=expression('log'[10]*' (Counts per survey)'),xlab='Year',main=paste(sp,GZ,sep=', '))
  lines(ts[1,]~c(seq(1993,2018)),col='navy')
  points(ts[1,]~c(seq(1993,2018)),col='white',pch=21,bg='navy',cex=1.2)
  lines(ts[2,]~c(seq(1993,2018)),col='darkred')
  points(ts[2,]~c(seq(1993,2018)),col='white',pch=21,bg='darkred',cex=1.2)
  legend(2016,c(max(na.omit(c(ts)))*1.05),c('RVC','REEF'),text.col=c('navy','darkred'),bty='n')
  dev.off(paste(paste(sp,GZ,sep='_'),'.pdf',sep=''))
}

setwd("C:/Users/14388/Desktop/reef_florida_keys_data/Key Largo time-series")
for(i in 1:nrow(mars_3403_stan)){

  TS_plot(ts_comp[[i]],sp=rvc.green.3403.sp[i],GZ='Key Largo')
  TS_plot.pdf(ts=ts_comp[[i]],sp=rvc.green.3403.sp[i],GZ='Key Largo')
  dev.off()
  }
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
                                                                                          "x0",'log_lik'),control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 1000, chains = 4, iter = 2000, thin = 1)
test_1.1 <- rstan::stan(model_code = mod_1.1, data = list(y = ypos, 
                                                      TT = ncol(Y), N = nrow(Y), n_pos = n_pos, ts_id=ts_pos, col_indx_pos = col_indx_pos, 
                                                      row_indx_pos = row_indx_pos), pars = c("sd_q", "x", "sd_r",'a',"u",
                                                                                             "x0",'log_lik'),control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 1000, chains = 4, iter = 2000, thin = 1)



test_2 <- rstan::stan(model_code = mod_2, 
                      data = list(y = ypos, TT = ncol(Y), N = nrow(Y), n_pos = n_pos, col_indx_pos = col_indx_pos, 
                      row_indx_pos = row_indx_pos), pars = c("sd_q", "x", "sd_r","x0",'log_lik'),
                      control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 1000, chains = 4, iter = 2000, thin = 1)

test_3 <- rstan::stan(model_code = mod_3.1, 
                      data = list(y = ypos, TT = ncol(Y), N = nrow(Y), n_pos = n_pos, col_indx_pos = col_indx_pos, 
                                  row_indx_pos = row_indx_pos), pars = c("sd_q", "x", "sd_r","x0",'log_lik'),
                      control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 1000, chains = 4, iter = 2000, thin = 1)

                                                                                             
                                                                                             
                                                                                                                                                                                   
test_2.1 <- rstan::stan(model_code = mod_2.1, data = list(y = ypos, 
                                                      TT = ncol(Y), N = nrow(Y), n_pos = n_pos, col_indx_pos = col_indx_pos, 
                                                      row_indx_pos = row_indx_pos), pars = c("sd_q", "x", "sd_r", "u", 
                                                                                             "x0"),control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 1000, chains = 4, iter = 2000, thin = 1)

test_3<- rstan::stan(model_code = mod_3, data = list(y = ypos, 
                                                      TT = ncol(Y), N = nrow(Y), n_pos = n_pos, col_indx_pos = col_indx_pos, 
                                                      row_indx_pos = row_indx_pos), pars = c("sd_q", "x", "sd_r", 



library(shinystan)
launch_shinystan(test_2)

m2<- extract(test_2)




mars_3403_stan<- data.frame(SP=NA,loo1=NA,loo2=NA,loo3=NA,Q1=NA,Q2=NA,Q2.reef=NA,R1.rvc=NA,R1.reef=NA,R2.rvc=NA,R2.reef=NA,mAbund.rvc=NA,mAbund.reef=NA,sdAbund.rvc=NA,sdAbund.reef=NA,dAIC1=NA,dAIC2=NA,years.rvc=NA,years.reef=NA)
for(i in 1:length(rvc.green.3403.sp)){
  spp<- filter(fish_reef,commonname==rvc.green.3403.sp[i])
  ts_comp[[i]]<- t(log10(rvc_3403_green[[i]]$mean_ssu_abundance)) #extract logged time-series
  m<- match(spp$commonname,REEF.green.3403.sp)
  if(is.na(m)==T){next}
  ts_comp[[i]]<- rbind(ts_comp[[i]],t(log10(REEF_3403_green[[m]]$mean_site_abund[1:26])))
  rownames(ts_comp[[i]])<- c(paste(spp$scientificname,'rvc',sep=":"),paste(spp$scientificname,'REEF',sep=":"))
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
  
  loo_comp<- loo::loo_compare(loo1,loo2,loo3)
  
  mars_3403_stan[i,1]=spp$commonname
  mars_3403_stan[i,2]=loo_comp["model1",1]
  mars_3403_stan[i,3]=loo_comp["model1",2]  
  mars_3403_stan[i,4]=loo_comp["model2",1]
  mars_3403_stan[i,5]=loo_comp["model2",2]
  mars_3403_stan[i,6]=loo_comp["model3",1]
  mars_3403_stan[i,7]=median(params_1$sd_q) 
  mars_3403_stan[i,8]=sort(params_1$sd_q)[length(params_1$sd_q)*0.975]
  mars_3403_stan[i,9]=median(params_2$sd_q) 
  mars_3403_stan[i,10]=sort(params_2$sd_q)[length(params_1$sd_q)*0.975]
  mars_3403_stan[i,11]=median(params_3$sd_q[,1]) 
  mars_3403_stan[i,12]=sort(params_3$sd_q[,1])[length(params_3$sd_q[,1])*0.975]
  mars_3403_stan[i,13]=median(params_3$sd_q[,2]) 
  mars_3403_stan[i,14]=sort(params_3$sd_q[,2])[length(params_3$sd_q[,2])*0.975]
  mars_3403_stan[i,6]=params.1$parMean[4]
  mars_3403_stan[i,7]=params.3$parMean[3]
  mars_3403_stan[i,8]=params.3$parMean[4]
  mars_3403_stan[i,9]=params.1$parMean[2]
  mars_3403_stan[i,10]=params.1$parMean[3] 
  mars_3403_stan[i,11]=params.3$parMean[1]
  mars_3403_stan[i,12]=params.3$parMean[2]
  mars_3403_stan[i,13]=10^(mean(na.omit(ts_comp[[i]][1,])))
  mars_3403_stan[i,14]=10^(mean(na.omit(ts_comp[[i]][2,])))
  mars_3403_stan[i,15]=sd(10^(na.omit(ts_comp[[i]][1,])))
  mars_3403_stan[i,16]=sd(10^(na.omit(ts_comp[[i]][2,])))
  mars_3403_stan[i,17]=mars_3403_stan[i,4]-min(mars_3403_stan[i,4:5])
  mars_3403_stan[i,18]=mars_3403_stan[i,5]-min(mars_3403_stan[i,4:5])
  if(mars_3403_stan[i,17]==0){mars_3403_stan[i,19]=1}
  if(mars_3403_stan[i,18]==0){mars_3403_stan[i,19]=2}
  mars_3403_stan[i,20]=length(na.omit(ts_comp[[i]][1,]))
  mars_3403_stan[i,21]=length(na.omit(ts_comp[[i]][2,]))
  
  
  #TS.plot.pdf(x=ts_comp[[i]],m=mars_3403$mod[i],sp=mars_3403$cName[i],GZ='Key Largo')
  # dev.off()
  print(i)
}

for(i in 1:nrow())
  
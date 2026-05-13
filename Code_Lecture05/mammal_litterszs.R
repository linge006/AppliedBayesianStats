
setwd("...")

library(rstan)
#library(gridExtra)
rstan_options (auto_write = TRUE)
options (mc.cores = parallel::detectCores ())

###

poisson_lin_disp <- "
data {
  int<lower=1> N; int<lower=1> K;
  int<lower=0> y[N];
  int<lower=1,upper=K> x[N];
  real<lower=0> pr_b_mn; real pr_b_sd; real pr_c_shape;
}

parameters {
  vector[K] beta;
  vector[N] e_raw;
  real<lower=0> sigma_e;
}

transformed parameters {
  // linear predictor and residual and latent variable
  vector[N] e_res; vector[N] latent; 
  
  e_res = sigma_e * e_raw;
  latent = beta[x] + e_res;
}

model {
  // priors
  beta ~ normal(pr_b_mn, pr_b_sd);
  e_raw ~ normal(0, 1);
  //sigma_e ~ exponential(1);
  sigma_e ~ cauchy(0, pr_c_shape);  // half-normal
  
  // log link
  y ~ poisson_log(latent);
}

generated quantities {
  vector[N] log_lik;
  
  int ypred[N]; real pred_res;

  //beta_rs = beta[2];
  
  for (n in 1:N){
    pred_res = normal_rng(0, sigma_e);
    ypred[n] = poisson_log_rng(fmin(latent[n],20) + pred_res);
    log_lik[n] = poisson_log_lpmf(y[n] | latent[n]);
  } // End for loop n
}
"

poisson_lin_disp_mod <- stan_model(model_code=poisson_lin_disp, model_name="poisson_lin_disp")

#

poisson_lin <- "
data {
  int<lower=1> N; int<lower=1> K;
  int y[N];
  int<lower=1,upper=K> x[N];
  real<lower=0> pr_b_sd;
}

parameters {
  real beta[K];
}

model {
  // priors
  //alpha ~ normal(0, pr_a_sd);
  beta ~ normal(10, pr_b_sd);
  
  // log link
  for (i in 1:N) {
    y[i] ~ poisson(beta[x[i]]); 
  }
}

generated quantities {
  real delta12; real delta13; real delta14; //real delta15; 
  real delta23; real delta24; //real delta25; 
  real delta34; //real delta35;
  //real delta45;
  
   delta12 = beta[1]-beta[2];  delta13 = beta[1]-beta[3];  delta14 = beta[1]-beta[4];  //delta15 = beta[1]-beta[5]; 
   delta23 = beta[2]-beta[3];  delta24 = beta[2]-beta[4];  //delta25 = beta[2]-beta[5]; 
   delta34 = beta[3]-beta[4];  //delta35 = beta[3]-beta[5];
   //delta45 = beta[4]-beta[5];
  
  real<lower=0> ypred[N]; real log_lik[N];
  for (n in 1:N){
    ypred[n] = poisson_rng(beta[x[n]]);
    log_lik[n] = poisson_lpmf(y[n] | beta[x[n]]);
  } // End for loop n
}
"

poisson_lin_mod <- stan_model(model_code=poisson_lin, model_name="poisson_lin")


# Read data
pig_lpb_dat <- read.csv("litters_per_breed.csv")
pig_lpp_dat <- read.csv("litters_per_parity.csv")
dog_lpb_dat <- read.csv("litterszs_dogbreeds.csv")

# Prepare list object for STAN
pig_lpb_dat_list <- list(N=nrow(pig_lpb_dat), K=length(unique(pig_lpb_dat$breed)), y=pig_lpb_dat$litter_size, 
                     x=as.integer(factor(pig_lpb_dat$breed)), pr_b_sd=1)


### Run model
poisson_lin_pig.rs1 <- sampling(poisson_lin_mod, data=pig_lpb_dat_list, iter=15e3, warmup=1e3, chains=4, thin=10)
print(summary(poisson_lin_pig.rs1, pars=c("beta"))$summary)
print(summary(poisson_lin_pig.rs1, pars=c("delta12","delta13","delta14","delta23","delta24","delta34"))$summary)

# Adjust prior to sd=10 and rerun model
pig_lpb_dat_list$pr_b_sd <- 10 # prior SD for alpha and beta
poisson_lin_pig.rs2 <- sampling(poisson_lin_mod, data=pig_lpb_dat_list, iter=15e3, warmup=1e3, chains=4, thin=10)
print(summary(poisson_lin_pig.rs2, pars=c("beta"))$summary)

# Adjust prior to sd=20 and rerun model
pig_lpb_dat_list$pr_b_sd <- 20 # prior SD for alpha and beta
poisson_lin_pig.rs3 <- sampling(poisson_lin_mod, data=pig_lpb_dat_list, iter=15e3, warmup=1e3, chains=4, thin=10)
print(summary(poisson_lin_pig.rs3, pars=c("beta"))$summary)

# Write simulation output to csv files
write.csv(summary(poisson_lin_pig.rs1, 
                  pars=c("beta","delta12","delta13","delta14",
                         "delta23","delta24","delta34",
                         "log_lik","ypred"))$summary, "pig_lpb_lin_p1.csv")
write.csv(summary(poisson_lin_pig.rs2, 
                  pars=c("beta","delta12","delta13","delta14",
                         "delta23","delta24","delta34",
                         "log_lik","ypred"))$summary, "pig_lpb_lin_p2.csv")
write.csv(summary(poisson_lin_pig.rs3, 
                  pars=c("beta","delta12","delta13","delta14",
                         "delta23","delta24","delta34",
                         "log_lik","ypred"))$summary, "pig_lpb_lin_p3.csv")

# Write traceplots to pdf files
pdf("lpb_lin.pdf", width=11, height=8.5)
traceplot(poisson_lin_pig.rs1, pars=c("beta"))
traceplot(poisson_lin_pig.rs2, pars=c("beta"))
traceplot(poisson_lin_pig.rs3, pars=c("beta"))
dev.off()

# load ggmcmc and gridExtra for plotting
library(ggmcmc)
library(gridExtra)

# Generate ggs (ggmcmc) objects per fitted model
P_ggs1 <- ggs(poisson_lin_pig.rs1)
P_ggs2 <- ggs(poisson_lin_pig.rs2)
P_ggs3 <- ggs(poisson_lin_pig.rs3)

# Set axis sizes, etc
p_theme <- theme(text=element_text(size=11), axis.text=element_text(size=12), strip.text=element_text(size=16), 
                 legend.title=element_text(size=9))

# Generate density plots per parameter
p12 <- ggs_density(subset(P_ggs1,Parameter=="delta12")) + p_theme
p13 <- ggs_density(subset(P_ggs1,Parameter=="delta13")) + p_theme
p14 <- ggs_density(subset(P_ggs1,Parameter=="delta14")) + p_theme
p23 <- ggs_density(subset(P_ggs1,Parameter=="delta23")) + p_theme
p24 <- ggs_density(subset(P_ggs1,Parameter=="delta24")) + p_theme
p34 <- ggs_density(subset(P_ggs1,Parameter=="delta34")) + p_theme

p212 <- ggs_density(subset(P_ggs2,Parameter=="delta12")) + p_theme
p213 <- ggs_density(subset(P_ggs2,Parameter=="delta13")) + p_theme
p214 <- ggs_density(subset(P_ggs2,Parameter=="delta14")) + p_theme
p223 <- ggs_density(subset(P_ggs2,Parameter=="delta23")) + p_theme
p224 <- ggs_density(subset(P_ggs2,Parameter=="delta24")) + p_theme
p234 <- ggs_density(subset(P_ggs2,Parameter=="delta34")) + p_theme

p312 <- ggs_density(subset(P_ggs3,Parameter=="delta12")) + p_theme
p313 <- ggs_density(subset(P_ggs3,Parameter=="delta13")) + p_theme
p314 <- ggs_density(subset(P_ggs3,Parameter=="delta14")) + p_theme
p323 <- ggs_density(subset(P_ggs3,Parameter=="delta23")) + p_theme
p324 <- ggs_density(subset(P_ggs3,Parameter=="delta24")) + p_theme
p334 <- ggs_density(subset(P_ggs3,Parameter=="delta34")) + p_theme

# Save multipanel parameter plots
pdf("Pig_litters_breed.pdf", width=11, height=6.5)
grid.arrange(p12,p13,p14,p23,p24,p34, ncol=3)
grid.arrange(p212,p213,p214,p223,p224,p234, ncol=3)
grid.arrange(p312,p313,p314,p323,p324,p334, ncol=3)
dev.off()


### Over-dispersed model runs ###

# Set priors for over-dispersed model
pig_lpb_dat_list$pr_b_mn <- 12
pig_lpb_dat_list$pr_b_sd <- 1.05
pig_lpb_dat_list$pr_c_shape <- 0.25

# Run model, print output, write to csv
poisson_lin_disp_pig.rs11 <- sampling(poisson_lin_disp_mod, data=pig_lpb_dat_list, iter=10e3, warmup=1e3, chains=3, thin=10)
print(summary(poisson_lin_disp_pig.rs11, pars=c("beta","sigma_e"))$summary)

write.csv(summary(poisson_lin_disp_pig.rs11, 
                  pars=c("beta","sigma_e","lp__","log_lik","latent","e_raw","e_res","ypred"))$summary, 
          "pig_lpb_lin_disp.csv")


# Save traceplot as pdf
pdf("pig_lpb_lin_disp.pdf", width=11, height=8.5)
traceplot(poisson_lin_disp_pig.rs11, pars=c("beta","sigma_e"))
dev.off()

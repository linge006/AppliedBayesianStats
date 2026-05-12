
setwd("C:/Users/hvanling/OneDrive - University of Guelph/Documents/CourseInBayesian/R_code_lecture02")

sheep_FN_dat <- read.csv("sheep_FN_sim.csv")

### Generate data input for STAN model as a list
sheep_FN_dat <- list(N=nrow(sheep_FN_dat), y=sheep_FN_dat$y, 
                     J=max(sheep_FN_dat$Study), EXP=sheep_FN_dat$Study,
                     x_DMI=sheep_FN_dat$x_DMI, x_CP=sheep_FN_dat$x_CP, x_NDF=sheep_FN_dat$x_NDF)


library(rstan)
## Loading required package: StanHeaders
## rstan (Version 2.18.2, GitRev: 2e1f913d3ca3)
## For execution on a local, multicore CPU with excess RAM we recommend calling
## options(mc.cores = parallel::detectCores()).
## To avoid recompilation of unchanged Stan programs, we recommend calling
## rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

FN_HierLin <- "
  data {
    int<lower=1> N;               // number of data points
    real y[N];                    // the FN data
    real x_DMI[N];                // the covariate 1 data
    real x_CP[N];                 // the covariate 2 data
    real x_NDF[N];                // the covariate 3 data
    int<lower=1> J;               // number of Studies
    int<lower=1, upper=J> EXP[N]; // the Studies
  } // end of data

parameters{
  real beta_0;            // regression coefficients
  real beta[3];           // regression coefficients
  real u[J];              // experiment intercepts
  real<lower=0> sigma_u;  // experiment std dev
  real<lower=0> sigma_e;  // residual std dev
} // end of parameters

model {
  for(i in 1:N){
    y[i] ~ normal(beta_0 + beta[1]*x_DMI[i] + beta[2]*x_CP[i] + beta[3]*x_NDF[i] + u[EXP[i]], sigma_e);
  }

  beta_0 ~ normal(0,300);
  beta ~ normal(0,300);
  u ~ normal(0,sigma_u);
  sigma_u ~ cauchy(0,15); 
  sigma_e ~ cauchy(0,15);
} // end of model

"

FN_HierLinMod <- stan_model(model_code=FN_HierLin,model_name="FN_HierLin")


### Run the MCMC simulation of the STAN model
for_fitattime <- Sys.time()
sheep.rs2 <- sampling(FN_HierLinMod, data=sheep_FN_dat, 
                      iter=15e3, warmup=1e3, thin=25, chains=2,
                      pars=c("beta_0","beta","u","sigma_u","sigma_e")); times_for_mod <- list(for_fitattime,Sys.time())
print(times_for_mod)

print(summary(sheep.rs2, pars=c("beta_0","beta","u","sigma_u","sigma_e"))$summary)
### Write output to csv file
write.csv(summary(sheep.rs2, pars=c("beta_0","beta","sigma_u","sigma_e"))$summary, "sheep_FN_Hier_out.csv")

### Save traceplot in pdf format
#pdf("ModPars_vs_iterations.pdf", height=7, width=11)
traceplot(sheep.rs2, pars=c("beta_0","beta","sigma_u","sigma_e"))
#dev.off()

### Vectorized mixed-model

FN_HierLin_X <- "
  data {
    int<lower=0> N;               // number of data points
    real y[N];                    // the CH4 data
    int<lower=1> K;               // number of predictors
    matrix[N,K+1] x;                  // the N covariate dataframe
    int<lower=1> J;               //number of experiments
    int<lower=1, upper=J> EXP[N]; // the experiments
  } // end of data

parameters{
  vector[K+1] beta;           // slope rcs
  real<lower=0> sigma_e;    // residual std dev
  vector[J] u;              // experiment intercepts
  real<lower=0> sigma_u;    // experiment std dev
  } // end of parameters

model {
  y ~ normal(x*beta + u[EXP], sigma_e);
  
  beta ~ normal(0,300);
  u ~ normal(0,sigma_u);
  sigma_u ~ cauchy(0,15); 
  sigma_e ~ cauchy(0,15);
}

"

FN_HierLin_X_mod <- stan_model(model_code=FN_HierLin_X,model_name="FN_HierLin_X")

### Adjust data list for vectorized STAN model
sheep_FN_dat_X <- list(N=sheep_FN_dat$N, y=sheep_FN_dat$y, K=3,
                       x=data.frame(x_0=1, DMI=sheep_FN_dat$x_DMI, CP=sheep_FN_dat$x_CP, NDF=sheep_FN_dat$x_NDF),
                       J=sheep_FN_dat$J, EXP=sheep_FN_dat$EXP)

# Note first column of data.frame object x
print(head(sheep_FN_dat_X$x))

### Run the HMC simulation of the STAN model
X_fitattime <- Sys.time()
sheep.rs3 <- sampling(FN_HierLin_X_mod, data=sheep_FN_dat_X, 
                      iter=15e3, warmup=1e3, thin=25, chains=2,
                      pars=c("beta","u","sigma_u","sigma_e")); times_X_mod <- list(X_fitattime,Sys.time())
print(list(for_mod=times_for_mod,X_mod=times_X_mod))

print(summary(sheep.rs3, pars=c("beta","u","sigma_u","sigma_e"))$summary)

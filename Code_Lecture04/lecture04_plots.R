setwd("C:/Users/hvanling/OneDrive - University of Guelph/Documents/CourseInBayesian/R_code_lecture04")

worm_Zn_dat <- read.table("Efetida_Zn.txt", header=T, sep=",") # Read dataset
print(head(worm_Zn_dat)) # print head of dataset

# Add a column with 1's and 2's for the uptake and elimination phases
t_n <- 14 # Set transfer time
worm_Zn_dat <- data.frame(worm_Zn_dat, Phase=1+as.integer(worm_Zn_dat$Time > t_n)) 

# Plot the data 
pdf("Earthworms_Zn_vs_Time.pdf", height=4.75, width=8)
par(mar=c(4.25,4.5,0.5,0.5))
plot(toxicant~Time, data=worm_Zn_dat[worm_Zn_dat$Time <= t_n,], pch=16, ylab=bquote("Zinc content in earthworms (" * mu * "g/g)"), xlab="Time (days)",
     ylim=c(0,350), xlim=c(0,28), xaxt="n", cex.axis=1.35, cex.lab=1.4, col="red") 
axis(1, at=seq(0,28,7), labels=seq(0,28,7), cex.axis=1.35)
points(toxicant~Time, data=worm_Zn_dat[worm_Zn_dat$Time > t_n,], pch=16, col="dodgerblue") 
legend("topright", legend=c("Uptake phase","Elimination phase"), pch=16, col=c("red","dodgerblue"), cex=1.25)
dev.off()

#####

library(ggplot2)
library(ggmcmc)
library(rstan)
## Loading required package: StanHeaders
## rstan (Version 2.18.2, GitRev: 2e1f913d3ca3)
## For execution on a local, multicore CPU with excess RAM we recommend calling
## options(mc.cores = parallel::detectCores()).
## To avoid recompilation of unchanged Stan programs, we recommend calling
## rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

library(loo)

mod_NLin_ToxKin <- "
// The data block contains all data that is read in from R
data {
  int<lower=0> N; // Number of observations
  real x[N];  // time
  real y[N];  // Zn concentration in earthworms
  real z; int t_d; // Zn concentration in soil and soil transfer time
  int grp[N]; int no_sigma; // Sigma stratification variable and number of sigma's
  real a_C0; real b_C0; // Prior parameters for C_0
  real mean_k1; real sd_k1; // Prior parameters for k_1
  real mean_k2; real sd_k2; // Prior parameters for k_2
  
} // End of data block


// The parameters block contains all model parameters that need to be estimated
parameters {
  real<lower=0> C_0; 
  real<lower=0> k[2]; // Contains k1 and k2. 
  real<lower=0> sigma_e[no_sigma]; // residual sd
} // End of parameters block


// The transformed parameters block contains the one-compartmental function for toxicokinetic modeling
transformed parameters {
  real m[N];
  
  for (i in 1:N){
    m[i] = (x[i] <= t_d ? k[1]/k[2]*z*(1-exp(-k[2]*x[i])) : k[1]/k[2]*z*(exp(-k[2]*(x[i]-t_d))-exp(-k[2]*x[i])));
  } // End of for loop i
} // End of transformed parameters block


// The model block contains the priors for all parameters and the likelihood of the models 
model {
  // priors
  C_0 ~ normal(a_C0, b_C0); // Assign a gamma prior to C_0
  k[1] ~ normal(mean_k1, sd_k1); k[2] ~ normal(mean_k2, sd_k2); // Assign normal priors to k1 and k2 
  sigma_e ~ cauchy(0,5); // Assign a half-cauchy prior to residual variance parameter sigma
  
  // likelihoods
  for (j in 1:N) y[j] ~ normal(C_0 + m[j], sigma_e[grp[j]]);// Toxkin model 
} // End of model block


// The generated quantities block contains all quantities to be returned in addition to the estimated parameters
generated quantities {

  // Initialize (additional) quantities to be generated
  vector[N] log_lik; // to use loo package and compute waic - the logLikelihood MUST be called log_lik!
  real pred[N]; // Quantities for fitted values (not including noise/error estimate)
  real res[N];  // Residuals
  
  // Assign values to the initialized quantities to be returned
  for (n in 1:N){
    pred[n] = C_0 + (x[n] <= t_d ? k[1]/k[2]*z*(1-exp(-k[2]*x[n])) : k[1]/k[2]*z*(exp(-k[2]*(x[n]-t_d))-exp(-k[2]*x[n])));
    res[n] = y[n] - pred[n];
    
    log_lik[n] = normal_lpdf(y[n] | pred[n], sigma_e[grp[n]]);   // Toxkin model
    
    } // End of for loop n
    
} // End of generated quantities block

"

#

print("Compiling one-compartmental stan model. Note that this may take a couple of minutes...!")
stan_mod_NLin_ToxKin <- stan_model(model_code=mod_NLin_ToxKin, model_name="mod_NLin_ToxKin") # Compile stan model

#####

worm_dat.list <- list(N=nrow(worm_Zn_dat), #Number of observations
                      x=worm_Zn_dat$Time,  # time
                      y=worm_Zn_dat$toxicant,  # Zn concentration in earthworms
                      z=mean(worm_Zn_dat$C_exp), # Zn conc in soil
                      t_d=t_n, # soil transfer time
                      grp=rep(1,nrow(worm_Zn_dat)), # Sigma stratification variable
                      no_sigma=1, # number of sigma's
                      a_C0=75, b_C0=5, # Prior parameters for C_0
                      mean_k1=1, sd_k1=1, # Prior parameters for k_1
                      mean_k2=1, sd_k2=1) # Prior parameters for k_2

worm_Zn_homosc.rs <- sampling(stan_mod_NLin_ToxKin, data=worm_dat.list, iter=15e3,
         pars=c("C_0","k","sigma_e","res","log_lik"), warmup=1e3, thin=20, chains=4)

# Update data object for heteroscedastic model
worm_dat.list$grp <- worm_Zn_dat$Phase
worm_dat.list$no_sigma <- 2
worm_Zn_heterosc.rs <- sampling(stan_mod_NLin_ToxKin, data=worm_dat.list, iter=15e3,
                           pars=c("C_0","k","sigma_e","res","log_lik"), warmup=1e3, thin=20, chains=4)

# Print/Inspect output
print(summary(worm_Zn_homosc.rs, pars=c("C_0","k","sigma_e"))$summary)
print(summary(worm_Zn_heterosc.rs, pars=c("C_0","k","sigma_e"))$summary)

# Generate traceplots
traceplot(worm_Zn_homosc.rs, pars=c("C_0","k","sigma_e"))
traceplot(worm_Zn_heterosc.rs, pars=c("C_0","k","sigma_e"))

### ggmcmc 

# Generate ggs objects for plotting
S1 <- ggs(worm_Zn_homosc.rs)
S2 <- ggs(worm_Zn_heterosc.rs)

# Update plotting settings (e.g. text sizes on plots, etc)
p_theme <- theme(text=element_text(size=13), axis.text=element_text(size=14), strip.text=element_text(size=16), 
                 legend.title=element_text(size=10))

# Generate traceplots and density plots - assign to objects g1 and g2
g1 <- ggs_traceplot(S1,family="sigma_e") + p_theme
g2 <- ggs_density(S1,family="sigma_e") + p_theme

# Generate density plots
p_0 <- ggs_density(S1,family="C_0") + p_theme
p_1 <- ggs_density(subset(S1,Parameter=="k[1]"),greek=T) + p_theme
p_2 <- ggs_density(subset(S1,Parameter=="k[2]"),greek=T) + p_theme
p_s <- ggs_density(S1,family="sigma_e") + p_theme

# Load gridExtra library for multiple panel plotting
library(gridExtra)

# Write the density plots to pdf files
pdf("Density_vs_Chain_worm_Zn.pdf", height=5, width=11)
grid.arrange(g1,g2,ncol=2)
dev.off()
pdf("Densities_Zn_worms.pdf", height=6, width=11)
grid.arrange(p_0,p_1,p_2,p_s,ncol=2)
dev.off()


### Alternative plotting option using bayesplot package

library(bayesplot)

posterior_homo <- as.array(worm_Zn_homosc.rs)
mcmc_dens_overlay(posterior_homo, pars = c("C_0", "k[1]", "k[2]","sigma_e[1]"))


### Function to compute loo

compute_loo <- function(stan_out.rs){
  loglik_mod.rs <- extract_log_lik(stan_out.rs, parameter_name="log_lik", merge_chains=T) # Extract log likelihood; merge_chains=T indicates that log-likelihoods should be merged across different chains if applicable
  
  LL <- as.array(stan_out.rs, pars="log_lik") # Convert log Likelihood into an array
  r_eff <- relative_eff(exp(LL), cores=1) # compute likelihood rather than log-likelihood and the relative efficiency of the MCMC chains based on these likelihoods. Calculation should be done using a single core.
  ltest <- loo(loglik_mod.rs, r_eff=r_eff, cores=1) # Compute the LOOIC using extracted log-likelihood and r_eff using a single core (cores=1).
  return(ltest) # Return LOOIC
} # End of compute_loo()

# Compute WAIC/looic for the two models
homosc.loo <- compute_loo(worm_Zn_homosc.rs)
heterosc.loo <- compute_loo(worm_Zn_heterosc.rs)

# Print looic for the two models
print(homosc.loo$estimates)
print(heterosc.loo$estimates)


### STAN code for model that computes posterior predictive distribution (PDD)

mod_NLin_ToxKin02 <- "
// The data block contains all data that is read in from R
data {
  int<lower=0> N; // Number of observations
  real x[N];  // time
  real y[N];  // Zn concentration in earthworms
  real z; int t_d; // Zn concentration in soil and soil transfer time
  int grp[N]; int no_sigma; // Sigma stratification variable and number of sigma's
  int<lower=N> P; real pred_t[P]; int s_grp[P]; // ...
  real a_C0; real b_C0; // Prior parameters for C_0
  real mean_k1; real sd_k1; // Prior parameters for k_1
  real mean_k2; real sd_k2; // Prior parameters for k_2
  
} // End of data block


// The parameters block contains all model parameters that need to be estimated
parameters {
  real<lower=0> C_0; 
  real<lower=0> k[2]; // Contains k1 and k2. 
  real<lower=0> sigma_e[no_sigma]; // residual sd
} // End of parameters block


// The transformed parameters block contains the one-compartmental function for toxicokinetic modeling
transformed parameters {
  real m[N];
  
  for (i in 1:N){
    m[i] = (x[i] <= t_d ? k[1]/k[2]*z*(1-exp(-k[2]*x[i])) : k[1]/k[2]*z*(exp(-k[2]*(x[i]-t_d))-exp(-k[2]*x[i])));
  } // End of for loop i
} // End of transformed parameters block


// The model block contains the priors for all parameters and the likelihood of the models 
model {
  // priors
  C_0 ~ normal(a_C0, b_C0); // Assign a gamma prior to C_0
  k[1] ~ normal(mean_k1, sd_k1); k[2] ~ normal(mean_k2, sd_k2); // Assign normal priors to k1 and k2 
  sigma_e ~ cauchy(0,5); // Assign a half-cauchy prior to residual variance parameter sigma
  
  // likelihoods
  for (j in 1:N) y[j] ~ normal(C_0 + m[j], sigma_e[grp[j]]);   // Toxkin model 
} // End of model block


generated quantities {
  // Initialize (additional) quantities to be generated
  real log_lik[N]; // to use loo package and compute waic - the logLikelihood MUST be called log_lik!
  real pred[N]; // Quantities for fitted values (not including noise/error estimate)
  real res[N];  // Residuals
  real cpred[P]; real ypred[P]; // Predicted values (including noise/error estimate)
  
  // Assign values to the initialized quantities to be returned
  for (n in 1:N){
    pred[n] = C_0 + (x[n] <= t_d ? k[1]/k[2]*z*(1-exp(-k[2]*x[n])) : k[1]/k[2]*z*(exp(-k[2]*(x[n]-t_d))-exp(-k[2]*x[n])));
    res[n] = y[n] - pred[n];
    
    log_lik[n] = normal_lpdf(y[n] | pred[n], sigma_e[grp[n]]); // Toxkin model
    
    } // End of for loop n
    
    for (p in 1:P){
    cpred[p] = (pred_t[p] <= t_d ? k[1]/k[2]*z*(1-exp(-k[2]*pred_t[p])) : k[1]/k[2]*z*(exp(-k[2]*(pred_t[p]-t_d))-exp(-k[2]*pred_t[p])));

    ypred[p] = normal_rng(C_0 + cpred[p], sigma_e[s_grp[p]]); // Toxkin model
    
  } // End of for loop p
} // End of generated quantities block

"

#

print("Compiling one-compartmental stan model. Note that this may take a couple of minutes...!")
stan_mod_NLin_ToxKin02 <- stan_model(model_code=mod_NLin_ToxKin02, model_name="mod_NLin_ToxKin02") # Compile stan model

# Set times at which PDD is computed
t_rng <- seq(1,t_n,0.1) # From 1 to t_n (=14 is this case) 

# rng_grp contains choice for sigma_e per data point of the PDD
rng_grp <- rep(1,length(t_rng)); t_rng <- c(t_rng,t_n+t_rng); rng_grp <- c(rng_grp,rng_grp+1)

worm_dat.list02 <- list(N=nrow(worm_Zn_dat), #Number of observations
                        x=worm_Zn_dat$Time,  # time
                        y=worm_Zn_dat$toxicant,  # Zn concentration in earthworms
                        z=mean(worm_Zn_dat$C_exp), # Zn conc in soil
                        t_d=t_n, # soil transfer time
                        grp=rep(1,nrow(worm_Zn_dat)), # Sigma stratification variable
                        no_sigma=1, # number of sigma's
                        P=length(t_rng), pred_t=t_rng, s_grp=rep(1,length(t_rng)),
                        a_C0=75, b_C0=5, # Prior parameters for C_0
                        mean_k1=1, sd_k1=1, # Prior parameters for k_1
                        mean_k2=1, sd_k2=1) # Prior parameters for k_2

# Run homoscedastic model with PDD
worm_Zn_homosc.rs02 <- sampling(stan_mod_NLin_ToxKin02, data=worm_dat.list02, iter=10e3,
                                pars=c("C_0","k","sigma_e","log_lik","ypred"), warmup=1e3, thin=25, chains=4)

worm_dat.list02$grp <- worm_Zn_dat$Phase # Update phase variable (sigma = {1,2}) for uptake and elimination phases (length=N) 
worm_dat.list02$s_grp <- rng_grp # Update phase variable (sigma = {1,2}) for uptake and elimination phases (length=P)
worm_dat.list02$no_sigma <- 2 # Number of sigma's (now 2, was 1 for the homoscedastic case)
worm_Zn_heterosc.rs02 <- sampling(stan_mod_NLin_ToxKin02, data=worm_dat.list02, iter=10e3,
                                  pars=c("C_0","k","sigma_e","log_lik","ypred"), warmup=1e3, thin=25, chains=4)

# Extract data summary
summ_homosc02 <- summary(worm_Zn_homosc.rs02, pars=c("C_0","k","sigma_e","ypred"))$summary
print(summary(worm_Zn_heterosc.rs02, pars=c("C_0","k","sigma_e"))$summary)

# Generate traceplots
traceplot(worm_Zn_homosc.rs02, pars=c("C_0","k","sigma_e"))
traceplot(worm_Zn_heterosc.rs02, pars=c("C_0","k","sigma_e"))

### Plot Posterior predictive distribution

source("plot_PostPredDist.R") # load function to plot PDD

pdf("PostPredDistr_worms_Zn.pdf", width=8, height=4)
par(mfrow=c(1,2), mar=c(1.5,3,1.5,0.5), oma=c(3,3,1,1))
Plot_PostPredDistr(summary(worm_Zn_homosc.rs02, pars=c("C_0","k","sigma_e","ypred"))$summary, 
                   summary(worm_Zn_heterosc.rs02, pars=c("C_0","k","sigma_e","ypred"))$summary, 
                   raw_data=worm_Zn_dat, "Time","toxicant", "", y_lim=c(-15,355), main02="Heteroscedastic")
mtext("Time (days)",1, outer=T, cex=1.35, line=1)
mtext(bquote("Zinc content (" * mu * "g/g)"), 2, outer=T, cex=1.35, line=0.5)
dev.off()

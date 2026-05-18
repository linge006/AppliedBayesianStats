
setwd("...")

library(rstan)
#library(gridExtra)
rstan_options (auto_write = TRUE)
options (mc.cores = parallel::detectCores ())


### Poisson+over-dispersion model

poisson_linpred_disp <- "
data {
  int<lower=1> N; int<lower=1> K;
  int<lower=0> y[N];
  matrix[N,K] x;
  real<lower=0> prior_sd; real prior_alpha; real prior_beta;
}

parameters {
  real alpha;
  vector[K] beta;
  vector[N] e_raw;
  real<lower=0> sigma_e;
}

transformed parameters {
  // linear predictor and residual and latent variable
  vector[N] e_res; vector[N] eta; vector[N] latent; 
  
  e_res = sigma_e * e_raw;
  eta = alpha + x * beta;
  latent = eta + e_res;
}

model {
  // priors
  alpha ~ normal(prior_alpha, prior_sd);
  beta ~ normal(prior_beta, prior_sd);
  e_raw ~ normal(0, 1);
  //sigma_e ~ exponential(1);
  sigma_e ~ normal(0, 0.1);  // half-normal
  
  // log link
  y ~ poisson_log(latent);
}

generated quantities {
  vector[N] log_lik;
  
  int ypred[N]; real pred_res; vector[N] eta_rs;

  eta_rs = alpha + x[,1] * beta[1];
  
  for (n in 1:N){
    pred_res = normal_rng(0, sigma_e);
    ypred[n] = poisson_log_rng(fmin(eta_rs[n],20) + pred_res);
    log_lik[n] = poisson_log_lpmf(y[n] | eta[n]+e_res[n]);
  } // End for loop n
}
"

poisson_linpred_disp_mod <- stan_model(model_code=poisson_linpred_disp, model_name="poisson_linpred_disp")


### Negative-Binomial model

NegBin_LinPred <- "
data {
  int<lower=1> N;
  int<lower=1> K;

  int<lower=0> y[N];
  matrix[N,K] x;

  real<lower=0> prior_sd;
  real prior_alpha;
  real prior_beta;
}

parameters {
  real alpha;
  vector[K] beta;

  // NB overdispersion parameter
  real<lower=0> phi;
}

transformed parameters {
  vector[N] eta;

  eta = alpha + x * beta;
}

model {
  // priors
  alpha ~ normal(prior_alpha, prior_sd);
  beta ~ normal(prior_beta, prior_sd);

  // dispersion prior
  phi ~ gamma(2, 0.1); // phi ~ exponential(1);

  // likelihood
  y ~ neg_binomial_2_log(eta, phi);
}

generated quantities {
  vector[N] log_lik;
  int ypred[N]; vector[N] eta_rs;

  eta_rs = alpha + x[,1] * beta[1];

  for (n in 1:N) {
    log_lik[n] = neg_binomial_2_log_lpmf(y[n] | eta[n], phi);

    ypred[n] = neg_binomial_2_log_rng(fmin(eta_rs[n], 16), phi);
  }
}
"

NegBin_LinPred_mod <- stan_model(model_code=NegBin_LinPred, model_name="NegBin_LinPred")

###

# Load Swedish traffic incidents dataset

data(Traffic, package = "MASS")

Traffic <- Traffic[order(Traffic$day),]
print(head(Traffic))

# Create data list object
Traffic_data_list <- list(N=nrow(Traffic), K=2, y=Traffic[,4], x=Traffic[,2:3])

# Convert no/yes into 0's and 1's, respectively
Traffic_data_list$x$limit <- as.integer(Traffic_data_list$x$limit)-1


# Set priors
Traffic_data_list$prior_sd <- 1 # prior SD for alpha and beta
Traffic_data_list$prior_alpha <- 3.0
Traffic_data_list$prior_beta <- 0.0

# Run simulation
poisson_linpred_disp.rs1 <- sampling(poisson_linpred_disp_mod, data=Traffic_data_list, iter=10e3, warmup=1e3, chains=4, thin=10)

write.csv(summary(poisson_linpred_disp.rs1, 
                  pars=c("alpha","beta","lp__", "eta","eta_rs","e_raw","e_res","log_lik","ypred"))$summary, 
          "poisson_linpred_disp.csv")

print(summary(poisson_linpred_disp.rs1, pars=c("alpha","beta","sigma_e","lp__"))$summary)
print(summary(poisson_linpred_disp.rs1, pars=c("eta_rs","e_raw","e_res","log_lik"))$summary)
print(summary(poisson_linpred_disp.rs1, pars=c("ypred"))$summary)

traceplot(poisson_linpred_disp.rs1, pars=c("alpha","beta","sigma_e","lp__"))


### NegBin Fit

Traffic_data_list$prior_sd <- 1.0
Traffic_data_list$prior_alpha <- 3.0
Traffic_data_list$prior_beta <- 0.0

NegBin_LinPred.rs1 <- sampling(NegBin_LinPred_mod, data=Traffic_data_list, iter=10e3, warmup=1e3, chains=4, thin=10)

write.csv(summary(NegBin_LinPred.rs1, pars=c("alpha","beta","phi","lp__","eta","y_rep"))$summary, 
          "NegBin_LinPred.csv")

print(summary(NegBin_LinPred.rs1, pars=c("alpha","beta","phi","lp__"))$summary)
print(summary(NegBin_LinPred.rs1, pars=c("log_lik","eta"))$summary)
print(summary(NegBin_LinPred.rs1, pars=c("y_rep"))$summary)

traceplot(NegBin_LinPred.rs1, pars=c("alpha","beta","phi","lp__"))


### Plotting using ggmcmc

library(ggmcmc)

Pois.ggs <- ggs(poisson_linpred_disp.rs1)
NegB.ggs <- ggs(NegBin_LinPred.rs1)

ggs_dens_Pois <- ggs_density(Pois.ggs,family="alpha|beta|sigma", greek=T)
ggs_dens_NegB <- ggs_density(NegB.ggs,family="alpha|beta|phi", greek=T)


library(gridExtra) # For multipanel plotting

pdf("ChainDens_Pois_NegBin.pdf", height=8, width=11)
grid.arrange(ggs_dens_Pois,ggs_dens_NegB,ncol=2)
dev.off()


### Plot_PostPredDistr() function to Plot Posterior predictive distribution ###

Plot_PostPredDistr <- function(summ_mod01, summ_mod02, raw_data, var_x, var_y, title, y_lim=c(-15,345), main01="Homoscedastic", main02="Heteroscedastic"){
  pars <- as.numeric(summ_mod01[1:3,1]); pars_hsk <- as.numeric(summ_mod02[1:3,1]) # Extract C_0, k_1 and k_2 for both models
  x_max <- max(raw_data[,var_x]) # Extract maximum value of raw data toxicant concentration in the organism
  #C_exp <- mean(raw_data$C_exp[raw_data[var_x,] <= 0.5*x_max]) # Compute mean of the toxicant concentration in the exposure medium
  print(x_max)
  t_rng <- raw_data[,var_x]#; t_rng <- c(t_rng,rev(t_rng)) # Generate time points at which error margin will be computed
  
  pred_rows <- grepl("ypred",rownames(summ_mod01)) # Select Predicted values (including noise/error estimate) for homoscedastic model
  pred_hsk_rows <- grepl("ypred",rownames(summ_mod02)) # Select Predicted values (including noise/error estimate) from heteroscedastic model
  
  c_y2 <- summ_mod01[pred_rows,"97.5%"]; c_y1 <- summ_mod01[pred_rows,"2.5%"] # Extract 2.5 and 97.5 percentiles of ypred for homoscedastic model
  c_y4 <- summ_mod02[pred_hsk_rows,"97.5%"]; c_y3 <- summ_mod02[pred_hsk_rows,"2.5%"] # Extract 2.5 and 97.5 percentiles of ypred for heteroscedastic model
  print('dog')
  # Model 1 plotting
  plot(NULL, xlim=c(0,x_max), ylim=y_lim, las=1,
       cex.lab=1.35, cex.axis=1.25, cex.main=1.35, ylab=title, xlab="", main=paste0(main01,""), xaxt='n') # Plot model's predicted values as a curve
  axis(1,seq(0,x_max,0.25*x_max),cex.axis=1.25) # Plot x-axis tick marks
  print('peer')
  polygon(c(t_rng, rev(t_rng)), c(c_y2, rev(c_y1)),
          col="lightgrey", lty=0) # Plot 95% credible interval of posterior predictive distribution as a grey shaded area
  points(raw_data[,var_x], raw_data[,var_y], pch=20) # Plot raw data points
  lines(raw_data[,var_x], summ_mod01[pred_rows,"mean"], col="cyan", lwd=1.5)
  abline(h=c(20,25,40,45,50), lty=2, col="darkgrey")
  
  # Model 2 plotting
  plot(NULL, xlim=c(0,x_max), ylim=y_lim, 
       cex.lab=1.35, cex.axis=1.25, cex.main=1.35, ylab="", xlab="", main=paste0(main02,""), xaxt='n', yaxt='n') # Plot model's predicted values as a curve
  axis(1,seq(0,x_max,0.25*x_max),cex.axis=1.25) # Plot x-axis tick marks
  polygon(c(t_rng, rev(t_rng)), c(c_y4, rev(c_y3)),
          col="lightgrey", lty=0) # Plot 95% credible interval of posterior predictive distribution as a grey shaded area
  points(raw_data[,var_x], raw_data[,var_y], pch=20) # Plot raw data points
  lines(raw_data[,var_x], summ_mod02[pred_rows,"mean"], col="magenta", lwd=1.5)
  abline(h=c(20,25,40,45,50), lty=2, col="darkgrey")
  
} # End Plot_PostPredDistr()

getwd()
pdf("PPD_Traffic_Pois_vs_NB.pdf", width=8, height=4)
par(mfrow=c(1,2), mar=c(1.0,2,1.5,0.25), oma=c(2.5,2,0.25,0.25))
Plot_PostPredDistr(summary(poisson_linpred_disp.rs1, pars=c("beta","sigma_e","ypred"))$summary, 
                   summary(NegBin_LinPred.rs1, pars=c("beta","phi","ypred"))$summary, 
                   raw_data=Traffic, "day","y", "", y_lim=c(0,60), main01="Poisson", main02="Negative-Binomial")
mtext("Time (days)",1, outer=T, cex=1.35, line=1.25)
mtext("Traffic incidents", 2, outer=T, cex=1.35, line=0.5)
dev.off()

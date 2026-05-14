
setwd("...")

### Calibration FME HIV virus model using STAN ###

# K. Soetaert and T. Petzoldt, 2010
# Inverse Modelling, Sensitivity and Monte Carlo Analysis in R Using Package FME
# 10.18637/jss.v033.i03


library(FME) # load deSolve package for solving ODEs

# Derivative specification
# Need to create a function including derivatives, name of function is HIV_FME
# t = times, initial = initial conditions, Parms = parameter values

HIV_FME <- function(Parms, State, Times, solver='lsoda'){
  HIV_FME2 <- function (times, initial, Parms){
    
    with(as.list(c(initial, Parms)),
         {
           if ((times %% 1) <= 1e-3 & times > 0.5) {cat('Simulation at',round(times,1),'hour...\n')} # Print marked times
           
           # Differentials
           dC <- lam - rho * C - bet * C * V
           dI <- bet * C * V - delt * I
           dV <- n * delt * I - c * V - bet * C * V
           
           # Return derivatives and other output
           return(list(
             c(dC, dI, dV),
             c(logV=log(V),dC=dC,dI=dI,dV=dV)))
         } # End of with(.)
    )
  } # End of derivative specification/HIV_FME2(.)
  
  cat('Model parameters are:\n'); print(Parms)
  
  #Solve model using ode for original solution
  cat('Simulation started at:'); print(Sys.time())
  out <- ode(y = State, times = Times, func = HIV_FME2, parms = Parms, method=solver)
  cat('Simulation phinished at:'); print(Sys.time())
  
  return(as.data.frame(out))
}  # End of HIV_FME model


# Initial parameter values
Parms_HIV <- c(bet=0.00002, rho=0.15, delt=0.55, c=5.5, lam=80, n=900)

# Initial values for state variables
c0 <- 100
v0 <- 50e3
i0 <- (-200750 + 5.5*v0) / (900*0.55)
initial_HIV <- c(C=c0, I=i0, V=v0)


# Set simulation time/step size
times_HIV <- c(seq(0, 0.8, 0.1), seq(2, 60, 2))


# Run model and assign output to object 
out_HIV <- HIV_FME(Parms=Parms_HIV, State=initial_HIV, Times=times_HIV, solver = 'lsodes')


# Data generation/simulation - error % of with respect to mean parameter value
error <- 0.10 # adjust this value to increase/decrease the noise level e.g. 0.05, 0.15, 0.25, etc


# build/simulate data
dat <- as.data.frame(out_HIV[out_HIV$time %in% c(2,4,6,seq(10,60,10)),1:4])
dat[,2] <- rnorm(length(dat[,2]),mean=dat[,2],sd=error*mean(dat[,2]))
dat[,3] <- rnorm(length(dat[,3]),mean=dat[,3],sd=error*mean(dat[,3]))
dat[,4] <- rnorm(length(dat[,4]),mean=dat[,4],sd=error*mean(dat[,4]))


# Plot model output and simulated data
par(mfrow=c(1,3), mar=c(3,4.5,1.5,0.5), oma=c(1,0,0,0))
plot(out_HIV$time, out_HIV$C, ylab="T", xlab="", type="l", ylim=c(100,375), cex.lab=1.35, cex.axis=1.35, cex.main=1.35,
     main="Model calibration to data")
points(C~time,data=dat, col='dodgerblue', pch=16, cex=1.35)
plot(out_HIV$time, out_HIV$I, ylab = "", xlab = "", type = "l", ylim=c(0,200), cex.axis=1.35, 
     cex.main=1.35, main="Model sensitivity to parameters")
points(I~time,data=dat, col='dodgerblue', pch=16, cex=1.35)
plot(out_HIV$time, out_HIV$V, ylab = "", xlab = "", type = "l", ylim=c(2e3,4e4), cex.axis=1.35, 
     cex.main=1.35, main="Model sensitivity to parameters")
points(V~time,data=dat, col='dodgerblue', pch=16, cex=1.35)
mtext("Time", 1, outer=T, cex=1.35)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### RSTAN coding ## ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

## Prepare data as a list for Stan using informative and weakly-informative priors
data_HIV_inf <-   list(n_obs=nrow(dat), n_diffeq=3, y0=initial_HIV, t0=0, ts=c(2,4,6,seq(10,60,10)), 
                       y=dat[,-1], sd_scale=0.25, cauchy_scale=5) # Informative prior list
data_HIV_wkinf <- list(n_obs=nrow(dat), n_diffeq=3, y0=initial_HIV, t0=0, ts=c(2,4,6,seq(10,60,10)), 
                       y=dat[,-1], sd_scale=..., cauchy_scale=...) # Set weakly-informative (or non-informative?) prior

### RSTAN HIV model

HIV_normal_mv <- "

functions {
  real[] HIV(real t, real[] y, real[] theta,
             real[] x_r, int[] x_i){
    
    real T = y[1];
    real I = y[2];
    real V = y[3];
    
    real bet = theta[1];
    real rho = theta[2];
    real delt = theta[3];
    real c = theta[4];
    real lam = theta[5];
    real n = 900.0;
    
    real dT_dt = lam - rho * T - bet * T * V;
    real dI_dt = bet * T * V - delt * I;
    real dV_dt = n * delt * I - c * V - bet * T * V;
    
    return{dT_dt, dI_dt, dV_dt};
  } // End of HIV function 
} // End functions block


data{
  int<lower=1> n_obs; // number of times sampled
  int<lower=1> n_diffeq; // number of differentials of system 
  real y0[n_diffeq]; // Initial values of state variables
  real t0; // Initial time point (zero)
  real ts[n_obs]; // time samples of observations
  row_vector[n_diffeq] y[n_obs];
  real sd_scale; real cauchy_scale;
} // End data block


transformed data {
  real x_r[0];
  int x_i[0];
} // End tf data block


parameters{
  real<lower=0> theta[5];
  vector<lower=0>[n_diffeq] sigma_e;  // init residual SD matrix
} // End parameters block


transformed parameters{
  real y_hat[n_obs, n_diffeq]; // Output from ODE solve
  row_vector[n_diffeq] mu[n_obs];      // init matrix for 'fixed-effects' regression coefficients
  
  y_hat = integrate_ode_rk45(HIV, y0, t0, ts, theta, x_r, x_i); // solve ODE part of model
  for (i in 1:n_obs) mu[i] = to_matrix(y_hat)[i,]; // transform ODE model output from real to row_vector
  
} // End tf parameters block


model {
  // priors
  theta[1] ~ normal(4.0e-5, sd_scale*4.0e-5); // truncated at 0 uniform(0,4e-5); // 
  theta[2] ~ normal(0.3, sd_scale*0.3); // truncated at 0 uniform(0,0.3); // 
  theta[3] ~ normal(1.0, sd_scale*1.0); // truncated at 0 uniform(0,1.1); // 
  theta[4] ~ normal(11.0, sd_scale*11.0); // truncated at 0 uniform(0,11); // 
  theta[5] ~ normal(100.0, sd_scale*100.0); // truncated at 0 uniform(0,160); // 
  sigma_e ~ cauchy(0, cauchy_scale);
  
  // Likelihood/sampling distribution
  y ~ multi_normal(mu, diag_matrix(sigma_e) ); // likelihood
} // End model block


generated quantities {
  vector[n_diffeq] pred_y[n_obs];
  pred_y = multi_normal_rng(mu, diag_matrix(sigma_e));
} // End generated quantities

"

###
library(rstan)
rstan_options (auto_write = TRUE)
options (mc.cores = parallel::detectCores ())

HIV_model_mv <- stan_model(model_code=HIV_normal_mv, model_name="HIV_normal_mv")


### Run model with informative prior ...
fit_HIV_mv_inf <- sampling(HIV_model_mv,
                       data=data_HIV_inf,
                       iter=10e3, warmup=1e3, thin=15, # 12, 6, 2
                       chains=4)

Post_mv_STAN_inf <- summary(fit_HIV_mv_inf, pars=c("theta","sigma_e","pred_y","lp__"))$summary
print(Post_mv_STAN_inf)

write.csv(Post_mv_STAN_inf,paste0("post_HIV_mv_inf_",error,".csv"))

pdf(paste0("traces_HIV_mv_inf_",error,".pdf"), height=8, width=11)
traceplot(fit_HIV_mv_inf, pars=c("theta","sigma_e")) + 
  theme(axis.text  = element_text(size = 12),
        axis.title = element_text(size = 13),
        strip.text = element_text(size = 13) )
dev.off()

### ggmcmc

library(ggmcmc); library(gridExtra)

HIV_mv_inf.ggs <- ggs(fit_HIV_mv_inf) # Generate ggs object for plotting

p_theme <- theme(text=element_text(size=14), axis.text=element_text(size=12), strip.text=element_text(size=16), 
                 legend.title=element_text(size=11)) # Set axis/label sizes

# Density plots
d01 <- ggs_density(subset(HIV_mv_inf.ggs,Parameter=="theta[1]")) + p_theme
d02 <- ggs_density(subset(HIV_mv_inf.ggs,Parameter=="theta[2]")) + p_theme
d03 <- ggs_density(subset(HIV_mv_inf.ggs,Parameter=="theta[3]")) + p_theme
d04 <- ggs_density(subset(HIV_mv_inf.ggs,Parameter=="theta[4]")) + p_theme
d05 <- ggs_density(subset(HIV_mv_inf.ggs,Parameter=="theta[5]")) + p_theme

pdf(paste0("HIV_densplot_inf_theta.pdf"), width=14, height=8)
grid.arrange(d01, d02, d03, d04, d05, ncol=3) # multipanel plot 
dev.off()

d11 <- ggs_density(subset(HIV_mv_inf.ggs,Parameter=="sigma_e[1]")) + p_theme
d12 <- ggs_density(subset(HIV_mv_inf.ggs,Parameter=="sigma_e[2]")) + p_theme
d13 <- ggs_density(subset(HIV_mv_inf.ggs,Parameter=="sigma_e[3]")) + p_theme

pdf(paste0("HIV_densplot_inf_sigma.pdf"), width=14, height=6)
grid.arrange(d11, d12, d13, ncol=3) # multipanel plot
dev.off()


### Posterior predictive distribution plotting

Parms_HIV[1:5] <- Post_mv_STAN_inf[1:5,c("mean")]
out_HIV_inf <- HIV_FME(Parms=Parms_HIV, State=initial_HIV, Times=times_HIV, solver = 'lsodes')


# Function to plot PDD

plot_post_pred <- function(Y_var, y_lim, x_lab, posterior=Post_mv_STAN, daytah=dat, out_HIV_mod=out_HIV, t_xyz=c(2,4,6,seq(10,60,10))){
  
  plot_mod <- paste0(Y_var,"~time"); plot(as.formula(plot_mod), out_HIV_mod, type="l", ylim=y_lim, xlab=x_lab, cex.lab=1.5, cex.axis=1.5)
  
  if (Y_var=="C") {DE_no <- ",1]"} else if (Y_var=="I") {DE_no <- ",2]"} else if (Y_var=="V" | Y_var=="logV") {DE_no <- ",3]"}
  
  pred_X <- grep(DE_no,rownames(posterior))
  c_y2 <- posterior[pred_X,"97.5%"]; c_y1 <- posterior[pred_X,"2.5%"]
  
  polygon(c(t_xyz, rev(t_xyz)), c(c_y2, rev(c_y1)),
          col="grey", lty=0)
  
  points(as.formula(plot_mod), data=daytah, type="p", pch=16)
  lines(as.formula(plot_mod), out_HIV_mod)
} # End plot_post_pred()


# Set y-axis bounds for C, I and V plotting
C_bounds <- c(100,400); I_bounds <- c(30,140); V_bounds <- c(3e3,10e3)

pdf(paste0("HIV_STAN_PostPred_inf_",error,".pdf"), height=6, width=10)
par(mfrow=c(1,3), oma=rep(0,4), mar=c(4.5,4.5,0.5,0.5))
plot_post_pred("C", y_lim=C_bounds, x_lab="", posterior=Post_mv_STAN_inf, out_HIV_mod=out_HIV_inf)
plot_post_pred("I", y_lim=I_bounds, x_lab="Time (h)", posterior=Post_mv_STAN_inf, out_HIV_mod=out_HIV_inf)
plot_post_pred("V", y_lim=V_bounds, x_lab="", posterior=Post_mv_STAN_inf, out_HIV_mod=out_HIV_inf)
legend('topright', c('Data points','50% quantile','95% Credible interval'), pch=c(16,-1,-1), lty=c(0,1,1), col=c('black','black','grey'), lwd=c(1,1,5))
dev.off()

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

### Run model with weakly informative prior ...

fit_HIV_mv_wkinf <- sampling(HIV_model_mv,
                           data=data_HIV_wkinf,
                           iter=2e3, warmup=1e3, thin=15, # 12, 6, 2
                           chains=1)

Post_mv_STAN_wkinf <- summary(fit_HIV_mv_wkinf, pars=c("theta","sigma_e","pred_y"))$summary
print(Post_mv_STAN_wkinf)

write.csv(Post_mv_STAN_wkinf,paste0("post_HIV_mv_wkinf_",error,".csv"))

pdf(paste0("traces_HIV_mv_wkinf_",error,".pdf"), height=8.5, width=11)
traceplot(fit_HIV_mv_wkinf, pars=c("theta","sigma_e"))
dev.off()

pdf(paste0("HIV_densplot_wkinf_",error,".pdf"), width=8.5, height=8.5)
stan_dens(fit_HIV_mv_wkinf, pars=c("theta","sigma_e"), fill="grey")
dev.off()


### Posterior predictive distribution plotting

Parms_HIV[1:5] <- Post_mv_STAN_wkinf[1:5,c("mean")]
out_HIV_wkinf <- HIV_FME(Parms=Parms_HIV, State=initial_HIV, Times=times_HIV, solver = 'lsodes')

pdf(paste0("HIV_STAN_PostPred_wkinf_",error,".pdf"), height=6, width=10)
par(mfrow=c(1,3), oma=rep(0,4), mar=c(4.5,4.5,0.5,0.5))
plot_post_pred("C", y_lim=C_bounds, x_lab="", posterior=Post_mv_STAN_wkinf, out_HIV_mod=out_HIV_wkinf)
plot_post_pred("I", y_lim=I_bounds, x_lab="Time (h)", posterior=Post_mv_STAN_wkinf, out_HIV_mod=out_HIV_wkinf)
plot_post_pred("V", y_lim=V_bounds, x_lab="", posterior=Post_mv_STAN_wkinf, out_HIV_mod=out_HIV_wkinf)
legend('topright', c('Data points','50% quantile','95% Credible interval'), pch=c(16,-1,-1), lty=c(0,1,1), col=c('black','black','grey'), lwd=c(1,1,5))
dev.off()

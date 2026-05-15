
### Calibration ODE model using STAN ###

# L?o Grinsztajn, Elizaveta Semenova, Charles C. Margossian, Julien Riou (2019). 
# Bayesian workflow for disease transmission modeling in Stan
# arXiv
# https://arxiv.org/pdf/2006.02985.pdf

setwd("...")
library(deSolve) # load deSolve package for solving ODEs


# Derivative specification
# Need to create a function including derivatives, name of function is COViD
# t = times, initial = initial conditions, Parms = parameter values

COViD <- function(Parms, State, Times, solver='lsoda')
{
  COViD2 <- function (times, initial, Parms){
    
    with(as.list(c(initial, Parms)),
         {
           if ((times %% 60) <= 1e-3 & times > 9) {cat('Simulation at',round(times,1),'min...\n')} # Print marked times
           
           # state variables
           
           f1 <- beta*S*I/N
           f2 <- gamma*I
           dS <- -f1
           dI <- f1 - f2
           dR <- f2

           # Return derivatives and other output
           return(list(
             c(dS, dI, dR),
             c(f1=f1, f2=f2)))
         } # End of with(.)
    )
  } # End of derivative specification/COViD2(.)
  
  cat('Model parameters are:\n'); print(Parms)
  
  #Solve model using ode for original solution
  cat('Simulation started at:'); print(Sys.time())
  out <- ode(y = State, times = Times, func = COViD2, parms = Parms, method=solver)
  cat('Simulation phinished at:'); print(Sys.time())
  
  return(as.data.frame(out))
}  # End of COViD model

###


Parms_C <- c(beta=1.72, gamma=0.53)

# Initial values
N <- 763
i0 <- 1
s0 <- N - i0
r0 <- 0
initial_C <- c(S=s0, I=i0, R=r0)

# Set simulation time/step size
end <- 20
times_C <- seq(0,end,0.1)


# Run model and assign output to object 
out_C <- COViD(Parms=Parms_C, State=initial_C, Times=times_C, solver='lsodes')
tail(out_C)


# Plotting
LogicSeq <- which(!grepl("time",colnames(out_C)))
#pdf('COViD_01.pdf', width=29.7/2.54, height=21/2.54)
par(mfrow=c(2,3), mar=c(4.5,4.5,2,0.5), oma=c(1,1,1,1)) # set 3 rows and 4 columns per plot 
sapply(LogicSeq, function(x) plot(out_C[,x]~out_C$time,type="l",xlab="Time (days)",ylab="",
                                  cex.lab=1.25,cex.axis=1.25,las=1, main=colnames(out_C[x]))) # plotting
#dev.off()


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

SIR_distribution <- "

functions {
  real[] sir(real t, real[] y, real[] theta,
             real[] x_r, int[] x_i){
    
    ...
    
    return{dS_dt, dI_dt, dR_dt};
  } // End of sir function
} // End functions block


data{
  int<lower=1> n_days; // Number of observations days
  real y0[3]; // Initial values of state variables
  real t0; // t at t_zero
  real ts[n_days]; // Time points at which observations were done
  int N;
  int cases[n_days];
} // End data block


transformed data {
  real x_r[0];
  int x_i[1] = {N};
} // End tf data block


parameters{
  real<lower=0> beta;
  real<lower=0> gamma;
  real<lower=0> theta_inv;
} // End parameters block


transformed parameters{
  real y[n_days, 3];
  real theta = .1 / theta_inv;
  {
    real theta[2];
    theta[1] = beta;
    theta[2] = gamma;
    
    y = integrate_ode_rk45(sir, y0, t0, ts, theta, x_r, x_i);
  }
} // End tf parameters block


model {
  // priors
  beta ~ distribution(1, 1); // truncated at 0
  gamma ~ distribution(1, 1); // truncated at 0
  theta_inv ~ distribution(1);
  
  // Likelihood/sampling distribution
  // col(matrix x, int n)
  cases ~ distribution(col(to_matrix(y), 2), theta);
} // End model block


generated quantities {
  real R0 = beta/gamma;
  real recovery_time = 1/gamma;
  real pred_cases[n_days];
  pred_cases = distribution_rng(col(to_matrix(y), 2) + 1e-5, theta);
} // End generated quantities

"


##### ##### 

library(tidyverse)
library(outbreaks)

# time series of cases
cases <- influenza_england_1978_school$in_bed # Number of students in bed

# times
n_days <- length(cases)
t <- seq(0, n_days, by = 1)
t0 <- 0
t <- t[-1]

# initial conditions
y0 = c(S = s0, I = i0, R = r0)

# data for Stan
data_SIR <- list(n_days = n_days, y0 = y0, t0 = t0, ts = t, N = N, cases = cases)


# Compile STAN model
library(rstan)
#library(gridExtra)
rstan_options (auto_write = TRUE)
options (mc.cores = parallel::detectCores ())
set.seed(3) # for reproducibility

SIR_model <- stan_model(model_code=SIR_distribution, model_name="SIR_distribution")

fit_SIR_distribution <- sampling(SIR_model,
                           data = data_SIR,
                           warmup = 1e3, iter = 2e3,
                           chains = 2)
ODE_out_rs <- summary(fit_SIR_distribution, pars=c("beta","gamma","theta_inv","theta","R0","recovery_time","pred_cases","lp__"))$summary

write.table(summary(fit_SIR_distribution)$summary, "SIR_ODE.csv", sep=",")
pdf("SIR_ODE_chns.pdf", width=11, height=8.5)
traceplot(fit_SIR_distribution, pars=c("beta","gamma","theta_inv","theta","R0","recovery_time","pred_cases","lp__"))
dev.off()

#

pairs(fit_SIR_distribution, pars=c("beta","gamma","theta_inv","theta","R0","recovery_time"))


### Plot posterior predictive distribution

c_y1 <- ODE_out_rs[grepl("pred_cases", rownames(ODE_out_rs)),"2.5%"]
c_y2 <- ODE_out_rs[grepl("pred_cases", rownames(ODE_out_rs)),"97.5%"]

#png("SIR_PDD.png", width=7, height=4, units="in", res=300)
tiff("SIR_PDD.tiff", width=7, height=4.5, units="in", res=300, compression="lzw")
par(mar=c(4.5,4.5,1,1))
plot(NULL, ylim=c(0,500), xlim=c(0,14), xlab="Time (days)", ylab="Number of infected people",
     cex.lab=1.25, cex.axis=1.25)
polygon(c(data_SIR$ts, rev(data_SIR$ts)), c(c_y2, rev(c_y1)),
        col="lightgreen", lty=0) # Plot 95% credible interval of posterior predictive distribution as a grey shaded area
lines(1:14,ODE_out_rs[grepl("pred_cases", rownames(ODE_out_rs)),"mean"], lwd=1.35,)
points(data_SIR$ts, data_SIR$cases, pch=16, col="dodgerblue")
dev.off()
  

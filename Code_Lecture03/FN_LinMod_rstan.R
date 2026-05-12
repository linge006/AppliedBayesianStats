
setwd("C:/Users/hvanling/OneDrive - University of Guelph/Documents/CourseInBayesian/R_code_lecture02")

#

FN_FixedEf <- "
data {
  int<lower=0> N;               // number of data points
  real y[N];                    // the Fecal nitrogen data
  real x_DMI[N];                  // the DMI data
  real x_CP[N];                  // the CP data
  real x_NDF[N];                  // the NDF data
} // End of data

parameters{
  real<lower=0> sigma_e;    // residual std dev
  real beta_0; real beta[3];             // regression coefficients
} // End of parameters

model {
  for(i in 1:N){
    y[i] ~ normal(beta_0 + beta[1]*x_DMI[i] + beta[2]*x_CP[i] + beta[3]*x_NDF[i], sigma_e);
  } // End of for loop i
  beta ~ normal(0,200);
  sigma_e ~ cauchy(0,15);
} // End of model

"

### RSTAN code ###

library(rstan)
## Loading required package: StanHeaders
## rstan (Version 2.18.2, GitRev: 2e1f913d3ca3)
## For execution on a local, multicore CPU with excess RAM we recommend calling
## options(mc.cores = parallel::detectCores()).
## To avoid recompilation of unchanged Stan programs, we recommend calling
## rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

### Initialize STAN code in R
FN_FixedEf_mod <- stan_model(model_code=FN_FixedEf,model_name="FN_FixedEf")

### Generate data input for STAN model as a list 
n_sample <- 160; n_studies <- 10

set.seed(6708)

# Simulate fecal Nitrogen data
sheep_FN_dat <- data.frame(x_DMI = rnorm(n_sample, 0.977, 0.209), # DMI
                           x_CP = rnorm(n_sample, 17.7, 4.33),   # CP
                           x_NDF = rnorm(n_sample, 41.6, 13.3),   # NDF
                           s_s = rep(rnorm(n_studies, 0, sqrt(1.68)), each=n_sample/n_studies),
                           e_res = rnorm(n_sample, 0, sqrt(2.25)) )


sheep_FN_dat <- data.frame(sheep_FN_dat, Study=rep(1:n_studies, each=n_sample/n_studies),
                           FN = -4.36 + 6.03*sheep_FN_dat$x_DMI + 0.109*sheep_FN_dat$x_CP + 
                             0.115*sheep_FN_dat$x_NDF + sheep_FN_dat$s_s + sheep_FN_dat$e_res)

# Simulate CH4 data
set.seed(6708)
sheep_CH4_dat <- data.frame(x_DMI = rnorm(n_sample, 0.977, 0.209), # DMI
                           x_NDF = rnorm(n_sample, 41.6, 13.3),   # NDF
                           x_BW = rnorm(n_sample, 41.9, 7.16),   # CP
                           s_s = rep(rnorm(n_studies, 0, sqrt(5.81)), each=n_sample/n_studies),
                           e_res = rnorm(n_sample, 0, sqrt(7.45)) )


sheep_CH4_dat <- data.frame(sheep_CH4_dat, Study=rep(1:n_studies, each=n_sample/n_studies),
                           CH4 = -4.00 + 12.9*sheep_CH4_dat$x_DMI + 0.121*sheep_CH4_dat$x_NDF + 
                             0.132*sheep_CH4_dat$x_BW + sheep_CH4_dat$s_s + sheep_CH4_dat$e_res)
library(nlme)
lme(CH4~x_DMI+x_NDF+x_BW, random=~1|s_s, data=sheep_CH4_dat)
lme(CH4~x_DMI+x_NDF+x_CP, random=~1|s_s, data=sheep_CH4_dat)

write.csv(sheep_FN_dat[,c("FN","Study","x_DMI","x_NDF","x_CP")], "sheep_FN_sim.csv", row.names=F)
write.csv(sheep_CH4_dat[,c("CH4","Study","x_DMI","x_NDF","x_BW")], "sheep_CH4_sim.csv", row.names=F)


sheep_FN_dat$Study <- as.factor(sheep_FN_dat$Study)

# Prep input data for STAN as a list
sheep_FN_dat <- list(N=n_sample, y=sheep_FN_dat$FN, 
                     x_DMI=sheep_FN_dat$x_DMI, x_CP=sheep_FN_dat$x_CP, x_NDF=sheep_FN_dat$x_NDF)

### Run the MCMC simulation of the STAN model
sheep.rs1 <- sampling(FN_FixedEf_mod, data=sheep_FN_dat, pars=c("beta_0","beta","sigma_e","lp__"),
                      iter=15e3, warmup=1e3, thin=25, chains=2)

### Write output to csv file
print(summary(sheep.rs1, pars=c("beta_0","beta","sigma_e","lp__"))$summary)
write.csv(summary(sheep.rs1, pars=c("beta_0","beta","sigma_e","lp__"))$summary, "sheep_FN_FixedEf_out.csv")

### Save traceplot in pdf format
pdf("ModPars_vs_iterations.pdf", height=6, width=9.5)
traceplot(sheep.rs1, pars=c("beta_0","beta","sigma_e"))+ theme(
  axis.text  = element_text(size = 12),
  axis.title = element_text(size = 12),
  strip.text = element_text(size = 12))
dev.off()

### Save traceplot in pdf format including warmup/burnin
pdf("sigma_vs_iterations.pdf", height=5, width=9.5)
traceplot(sheep.rs1, pars=c("sigma_e"), inc_warmup = TRUE) + theme(
            axis.text  = element_text(size = 14),
            axis.title = element_text(size = 14),
            strip.text = element_text(size = 14) )
dev.off()


### Plot parameter densities

library(ggmcmc)
library(gridExtra)

### Run the MCMC simulation of the STAN model again using 4 chains
sheep.rs2 <- sampling(FN_FixedEf_mod, data=sheep_FN_dat, pars=c("beta_0","beta","sigma_e","lp__"),
                      iter=15e3, warmup=1e3, thin=25, chains=4)
print(summary(sheep.rs2, pars=c("beta_0","beta","sigma_e","lp__"))$summary)

### Generate ggs object
S_sheep_02 <- ggs(sheep.rs2)

# Set axis/label sizes
shp_theme <- theme(text=element_text(size=12), axis.text=element_text(size=13), strip.text=element_text(size=16), 
                 legend.title=element_text(size=12))

# Plot density per parameter
s02_0 <- ggs_density(subset(S_sheep_02,Parameter=="beta_0"), greek=T) + shp_theme
s02_1 <- ggs_density(subset(S_sheep_02,Parameter=="beta[1]") ) + shp_theme
s02_2 <- ggs_density(subset(S_sheep_02,Parameter=="beta[2]") ) + shp_theme
s02_3 <- ggs_density(subset(S_sheep_02,Parameter=="beta[3]") ) + shp_theme
s02_5 <- ggs_density(S_sheep_02, family="sigma_e", greek=T) + shp_theme

# Plot a density panel per parameter on one pdf file
pdf("LinMod_ParDensities.pdf", height=6, width=11)
grid.arrange(s02_0, s02_1, s02_2, s02_3, 
             s02_5, ncol=3)
dev.off()

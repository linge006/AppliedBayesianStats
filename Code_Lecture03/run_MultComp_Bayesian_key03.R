library(rstan)
## Loading required package: StanHeaders
## rstan (Version 2.18.2, GitRev: 2e1f913d3ca3)
## For execution on a local, multicore CPU with excess RAM we recommend calling
## options(mc.cores = parallel::detectCores()).
## To avoid recompilation of unchanged Stan programs, we recommend calling
## rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
#options(mc.cores = NUM_CORES)
rstan_options(auto_write = TRUE)

library(loo)

model_rbd <- "
data {
  int<lower=1> N; // number of data points
  int<lower=1> no_trt;
  int trt[N]; // the treatments
  int<lower=1> no_block; // 
  int block[N]; // the variance groups
  vector[N] y; // the data
  real mean_y; real beta_sd; real sigma_b_sd; real sigma_e_sd;
  } // end of data

parameters{
  vector[no_trt] beta; // treatment means
  vector[no_block] blck;
  real<lower=1> nu_0; real<lower=1> nu_b;
  real<lower=0> sigma_b; // block std dev
  real<lower=0> sigma_e; // res std dev
  } // end of parameters

model {
  beta ~ normal(mean_y,beta_sd);
  blck ~ normal(0,sigma_b);
  nu_0 ~ exponential(1); nu_b ~ exponential(1);
  sigma_b ~ student_t(nu_b,0,sigma_b_sd);
  sigma_e ~ student_t(3,0,sigma_e_sd);
  
  //sigma_b ~ normal(0, sigma_b_sd) T[0, ]; sigma_e ~ normal(0, sigma_e_sd) T[0, ]; // adjust 5 to match your data scale
  //sigma_e ~ cauchy(0,15); //sigma_b ~ cauchy(0,15);

  y ~ student_t(nu_0, beta[trt] + blck[block], sigma_e);
} // End of model

generated quantities {
  real log_lik[N];
  
  vector[N] mu = beta[trt] + blck[block];
  
  for (n in 1:N){
    log_lik[n] = student_t_lpdf(y[n] | nu_0, mu[n], sigma_e);
  } // End for loop n
} // End of generated quantities
"

comp_mod_rbd <- stan_model(model_code=model_rbd, model_name="model_rbd")

###

setwd("C:/Users/hvanling/OneDrive - University of Guelph/Documents/CourseInBayesian/R_code_lecture03")

source("Group_assgn.R"); #source("Coef_comp.R"); source("Composite_fun.R"); source("Prob_comp.R")

MilkY_dat <- read.csv("MilkY_CattleBreeds_rbd.csv")


# To help see why we might consider a heteroscedastic model, let's get a box plot.
library(ggplot2)
library(dplyr)


p03 <- ggplot(MilkY_dat, aes(x=Breed, y=MY_FirstLact, group=Breed, fill=Breed) ) + 
  geom_boxplot() + ylab("Milk yield (kg/305 d)") +
  theme(axis.text=element_text(size=17), axis.title=element_text(size=17))
p03

pdf("Boxplot_MY_FirstLact.pdf", height=5, width=9); print(p03); dev.off()


##### #####

library(partitions)

max.grps <- length(unique(MilkY_dat$Breed))

# Generate treatment designs
parts <- setparts(max.grps)

# Generate treatment design matrix
trt_mat03 <- apply(parts, 2, grp.asn, npg=table(MilkY_dat$Breed))

mod.coefs01 <- matrix(-ceiling(sum(MilkY_dat$MY_FirstLact)), nrow=nrow(parts), ncol=ncol(parts))
cat("mod.coefs: \n"); print(mod.coefs01)
mod.sd01 <- matrix(-ceiling(sum(MilkY_dat$MY_FirstLact)), nrow=nrow(parts), ncol=ncol(parts))
cat("mod.sd: \n"); print(mod.sd01)

mod.coefs02 <- matrix(-ceiling(sum(MilkY_dat$MY_FirstLact)), nrow=nrow(parts), ncol=ncol(parts))
cat("mod.coefs: \n"); print(mod.coefs02)
mod.sd02 <- matrix(-ceiling(sum(MilkY_dat$MY_FirstLact)), nrow=nrow(parts), ncol=ncol(parts))
cat("mod.sd: \n"); print(mod.sd02)


# Run STAN model using Non-informative priors
Time.rs <- Sys.time(); fit.rs01 <- NULL
for (i in 1:ncol(trt_mat03)){
  reg.dat01 <- list(N=nrow(MilkY_dat), trt=trt_mat03[,i], no_trt=max(trt_mat03[,i]), block=MilkY_dat$Period, 
                    no_block=length(unique(MilkY_dat$Period)), y=MilkY_dat$MY_FirstLact,
                    mean_y=mean(MilkY_dat$MY_FirstLact), beta_sd=1e5, sigma_b_sd=1e5, sigma_e_sd=1e5)
  fit.rs01[[length(fit.rs01)+1]] <- sampling(comp_mod_rbd, data=reg.dat01, iter=102e3, warmup=2e3, thin=20, chains=3,
                                             pars=c("beta","nu_0","nu_b","blck","sigma_b","sigma_e","log_lik"))
  print(fit.rs01[[i]], pars=c("beta","nu_0","nu_b","sigma_b","sigma_e"))
  Time.rs <- list(Time.rs,Sys.time())
  mod_emm01 <- summary(fit.rs01[[i]],pars=c("beta"))$summary[,"mean"]
  mod_emm_sd01 <- summary(fit.rs01[[i]],pars=c("beta"))$summary[,"sd"]
  betas01 <- paste0("beta[",parts[,i],"]")
  if(i==1) {mod.coefs01[,1] <- mod_emm01; mod.sd01[,1] <- mod_emm_sd01
  } else {mod.coefs01[,i] <- mod_emm01[betas01]; mod.sd01[,i] <- mod_emm_sd01[betas01]}
} # End of for loop i


library(ggmcmc)

pdf("Pairs_NonInf_Fits01.pdf", height=7, width=7)
ggs_pairs(ggs(fit.rs01[[1]]), family = c("beta|nu_0|nu_b|sigma_e|sigma_b"))
ggs_pairs(ggs(fit.rs01[[2]]), family = c("beta|nu_0|nu_b|sigma_e|sigma_b"))
ggs_pairs(ggs(fit.rs01[[3]]), family = c("beta|nu_0|nu_b|sigma_e|sigma_b"))
ggs_pairs(ggs(fit.rs01[[4]]), family = c("beta|nu_0|nu_b|sigma_e|sigma_b"))
ggs_pairs(ggs(fit.rs01[[5]]), family = c("beta|nu_0|nu_b|sigma_e|sigma_b"))
dev.off()

ggs_traceplot(ggs(fit.rs01[[5]]), family = c("beta|nu_0|nu_b|sigma_e|sigma_b"))

pdf("Traces_logLik_fit01_MY_FirstLact.pdf", width=11, height=8.5)
traceplot(fit.rs01[[5]], pars="beta")
dev.off()


library(gridExtra)

pdf("Densities_logLik_fit01_MY_FirstLact.pdf", width=11, height=8.5)
grid.arrange(
  ggs_density(subset(ggs(fit.rs01[[5]]),Parameter==c("log_lik[1]","log_lik[2]","log_lik[3]","log_lik[4]"))),
  ggs_density(subset(ggs(fit.rs01[[5]]),Parameter==c("log_lik[5]","log_lik[6]","log_lik[7]","log_lik[8]"))),
  ggs_density(subset(ggs(fit.rs01[[5]]),Parameter==c("log_lik[9]","log_lik[10]","log_lik[11]","log_lik[12]"))), ncol=3)
grid.arrange(
  ggs_density(subset(ggs(fit.rs01[[5]]),Parameter==c("log_lik[13]","log_lik[14]","log_lik[15]","log_lik[16]"))),
  ggs_density(subset(ggs(fit.rs01[[5]]),Parameter==c("log_lik[17]","log_lik[18]","log_lik[19]","log_lik[20]"))),
  ggs_density(subset(ggs(fit.rs01[[5]]),Parameter==c("log_lik[21]","log_lik[22]","log_lik[23]","log_lik[24]"))), ncol=3)
dev.off()


# Run STAN model using Informative priors
Time.rs <- Sys.time(); fit.rs02 <- NULL
for (j in 1:ncol(trt_mat03)){
  reg.dat02 <- list(N=nrow(MilkY_dat), trt=trt_mat03[,j], no_trt=max(trt_mat03[,j]), block=MilkY_dat$Period, 
                  no_block=length(unique(MilkY_dat$Period)), y=MilkY_dat$MY_FirstLact,
                  mean_y=mean(MilkY_dat$MY_FirstLact), beta_sd=250, sigma_b_sd=250, sigma_e_sd=250)
  fit.rs02[[length(fit.rs02)+1]] <- sampling(comp_mod_rbd, data=reg.dat02, iter=102e3, warmup=2e3, thin=20, chains=3,
                                             pars=c("beta","nu_0","nu_b","blck","sigma_b","sigma_e","log_lik"))
  print(fit.rs02[[j]], pars=c("beta","nu_0","nu_b","sigma_b","sigma_e"))
  Time.rs <- list(Time.rs,Sys.time())
  mod_emm02 <- summary(fit.rs02[[j]],pars=c("beta"))$summary[,"mean"]
  mod_emm_sd02 <- summary(fit.rs02[[j]],pars=c("beta"))$summary[,"sd"]
  betas02 <- paste0("beta[",parts[,j],"]")
  if(j==1) {mod.coefs02[,1] <- mod_emm02; mod.sd02[,1] <- mod_emm_sd02
  } else {mod.coefs02[,j] <- mod_emm02[betas02]; mod.sd02[,j] <- mod_emm_sd02[betas02]}
} # End of for loop j


pdf("Pairs_Inf_Fits02.pdf", height=7, width=7)
ggs_pairs(ggs(fit.rs02[[1]]), family = c("beta|nu_0|nu_b|sigma_e|sigma_b"))
ggs_pairs(ggs(fit.rs02[[2]]), family = c("beta|nu_0|nu_b|sigma_e|sigma_b"))
ggs_pairs(ggs(fit.rs02[[3]]), family = c("beta|nu_0|nu_b|sigma_e|sigma_b"))
ggs_pairs(ggs(fit.rs02[[4]]), family = c("beta|nu_0|nu_b|sigma_e|sigma_b"))
ggs_pairs(ggs(fit.rs02[[5]]), family = c("beta|nu_0|nu_b|sigma_e|sigma_b"))
dev.off()

ggs_traceplot(ggs(fit.rs02[[2]]), family = c("beta|sigma_e|sigma_b"))


### Extract LOOIC values ###
compute_LOO <- function(fit100){
  LL <- extract_log_lik(fit100, parameter_name="log_lik", merge_chains=F) # Extract log-likelihood
  r_eff <- relative_eff(LL) # r_eff <- relative_eff(exp(LL))
  ltest <- loo(LL, r_eff=r_eff, moment_match = TRUE)
  print(Sys.time())
  ltest
} # End compute_LOOIC()

fit.rs_LOOIC01 <- lapply(fit.rs01, compute_LOO)
fit.rs_LOOIC02 <- lapply(fit.rs02, compute_LOO)

# Plot k value diagnostic plots
pdf("Pareto_k_rs12.pdf", height=5, width=7)
lapply(fit.rs_LOOIC01, plot)
lapply(fit.rs_LOOIC02, plot)
dev.off()

# Check if any k^hat exceeds 0.7
which(fit.rs_LOOIC01[[1]]$diagnostics$pareto_k > 0.7)
which(fit.rs_LOOIC02[[1]]$diagnostics$pareto_k > 0.7)

### Extract loo estimates
fit.rs_ex_LOO01 <- lapply(fit.rs_LOOIC01, function(x) x$estimates)
fit.rs_ex_LOO02 <- lapply(fit.rs_LOOIC02, function(x) x$estimates)


#### Extract model posterior probabilities ###

multcomp_prob <- function(fit.rs_LOOIC00, trt_mat00, parts00){
  mod.looics <- sapply(fit.rs_LOOIC00, function(x) x["looic","Estimate"]) # vector with just LOOIC values
  
  delta_mod.looics <- mod.looics - min(mod.looics)
  pre.probs <- exp(-0.5*delta_mod.looics)
  probs <- pre.probs/sum(pre.probs)
  names(probs) <- paste0("pi_",1:ncol(trt_mat00))
  
  char.mods <- apply(parts00, 2, paste, collapse="")
  
  #mod.looic <- mod.looics
  return(list(Pairing=char.mods, LOOICs=mod.looics, delta_B_k=delta_mod.looics, Pi_k=probs))
} # End of multocmp_prob()


Prob_B_k.rs01 <- multcomp_prob(fit.rs_ex_LOO01, trt_mat03, parts)
Prob_B_k.rs02 <- multcomp_prob(fit.rs_ex_LOO02, trt_mat03, parts)

write.csv(list2DF(Prob_B_k.rs01, nrow=5), "MY_FirstLact_NonInf.csv")
write.csv(list2DF(Prob_B_k.rs02, nrow=5), "MY_FirstLact_Inf.csv")


### Bayesian model averaging

BMA_comp <- function(mod.coefs00, mod.sd00, probs00, n_Treat=3){
  bmest.pre <- matrix(NA, nrow=nrow(setparts(n_Treat)), ncol=ncol(setparts(n_Treat)))
  
  for(p in 1:ncol(bmest.pre) ) {
    bmest.pre[,p] <- mod.coefs00[,p]*probs00[p]
  } # End for loop p
  
  bma.coefs <- rowSums(bmest.pre); 
  bma.coefs <- data.frame(mean=mod.coefs00 %*% as.matrix(probs00),
                          sd=mod.sd00 %*% as.matrix(probs00)); row.names(bma.coefs) <- paste0("Diet_",LETTERS[1:n_Treat])
  cat("Bayesian model averaging coefficients: \n"); print(bma.coefs)
  #write.csv(bma.coefs, "bma.coefs_MY_FirstLact.csv")
  return(bma.coefs)
} # BMA_comp()

BMA_pars01 <- BMA_comp(mod.coefs01, mod.sd01, Prob_B_k.rs01$Pi_k)
BMA_pars02 <- BMA_comp(mod.coefs02, mod.sd02, Prob_B_k.rs02$Pi_k)

# Write all (BMA) parameter outputs to files
write.table(rbind(mod.coefs01, mod.sd01), "MY_FirstLact_NonInf.csv", append=T, sep=",")
write.table(rbind(mod.coefs02, mod.sd02), "MY_FirstLact_Inf.csv", append=T, sep=",")

write.table(BMA_pars01, "MY_FirstLact_NonInf.csv", append=T, sep=",")
write.table(BMA_pars02, "MY_FirstLact_Inf.csv", append=T, sep=",")


##### Compute P-values per comparison ##### 

rs.pairs <- matrix(c(1,2, 1,3, 2,3), nrow=3, ncol=2, byrow=TRUE)

probs01 <- Prob_B_k.rs01$Pi_k
probs02 <- Prob_B_k.rs02$Pi_k

compute_Pr <- function(rs.pairs00,probs00,parts00){
  if(!all(is.na(rs.pairs00))){
    pw.bma.prob <- NULL
    plab <- NULL
    for (j in 1:dim(rs.pairs00)[1]){
      pairvec <- NULL
      for (i in 1:dim(parts00)[2]){
        pairvec[i] <- (parts00[rs.pairs00[j,1], i] == parts00[rs.pairs00[j,2], i])
      }
      pw.bma.prob[j] <- sum(probs00[pairvec])
      plab[j] <- paste(rs.pairs00[j,], collapse = "")
    }
    
    bma.df <- data.frame(Pairing = plab,
                         Prob    = pw.bma.prob)
    bma.df <- bma.df[order(bma.df[,2], decreasing = TRUE),]
    cat("Bayesian model averaging df \n"); print(bma.df)
  }
} # End of compute_Pr()


# P-values per comparison
Pr_pairs01 <- compute_Pr(rs.pairs,probs01,parts)
Pr_pairs02 <- compute_Pr(rs.pairs,probs02,parts)


### one_compartment() function ###

one_compartment <- function(C_0,k_1,k_2,C_exposure,time,t_e=t_d){
  
  ifelse(time <= t_e, 
         C_0+k_1/k_2*C_exposure*(1-exp(-k_2*time)), 
         C_0+k_1/k_2*C_exposure*(exp(-k_2*(time-t_e))-exp(-k_2*time)))
  
} # End one_compartment()

### Plot_PostPredDistr() function ###

Plot_PostPredDistr <- function(summ_mod01, summ_mod02, raw_data, var_x, var_y, title, y_lim=c(-15,345), main01="Homoscedastic", main02="Heteroscedastic"){
  pars <- as.numeric(summ_mod01[1:3,1]); pars_hsk <- as.numeric(summ_mod02[1:3,1]) # Extract C_0, k_1 and k_2 for both models
  x_max <- max(raw_data[,var_x]) # Extract maximum value of raw data toxicant concentration in the organism
  C_exp <- mean(raw_data$C_exp[raw_data$Time <= 0.5*x_max]) # Compute mean of the toxicant concentration in the exposure medium
  
  t_rng <- seq(min(raw_data$Time),0.5*x_max,0.1); t_rng <- c(t_rng,0.5*x_max+t_rng) # Generate time points at which error margin will be computed
  
  pred_rows <- grepl("ypred",rownames(summ_mod01)) # Select Predicted values (including noise/error estimate) for homoscedastic model
  pred_hsk_rows <- grepl("ypred",rownames(summ_mod02)) # Select Predicted values (including noise/error estimate) from heteroscedastic model
  
  c_y2 <- summ_mod01[pred_rows,"97.5%"]; c_y1 <- summ_mod01[pred_rows,"2.5%"] # Extract 2.5 and 97.5 percentiles of ypred for homoscedastic model
  c_y4 <- summ_mod02[pred_hsk_rows,"97.5%"]; c_y3 <- summ_mod02[pred_hsk_rows,"2.5%"] # Extract 2.5 and 97.5 percentiles of ypred for heteroscedastic model
  
  # Model 1 plotting
  plot(NULL, xlim=c(0,x_max), ylim=y_lim, las=1,
       cex.lab=1.35, cex.axis=1.25, cex.main=1.35, ylab=title, xlab="", main=paste0(main01,""), xaxt='n') # Plot model's predicted values as a curve
  axis(1,seq(0,x_max,0.25*x_max),cex.axis=1.25) # Plot x-axis tick marks
  polygon(c(t_rng, rev(t_rng)), c(c_y2, rev(c_y1)),
          col="lightgrey", lty=0) # Plot 95% credible interval of posterior predictive distribution as a grey shaded area
  points(raw_data[,var_x], raw_data[,var_y], pch=20) # Plot raw data points
  curve(one_compartment(pars[1],pars[2],pars[3],C_exp,x,0.5*x_max), add=T) # Plot model's predicted values as a curve again
  
  # Model 2 plotting
  plot(NULL, xlim=c(0,x_max), ylim=y_lim, 
       cex.lab=1.35, cex.axis=1.25, cex.main=1.35, ylab="", xlab="", main=paste0(main02,""), xaxt='n', yaxt='n') # Plot model's predicted values as a curve
  axis(1,seq(0,x_max,0.25*x_max),cex.axis=1.25) # Plot x-axis tick marks
  polygon(c(t_rng, rev(t_rng)), c(c_y4, rev(c_y3)),
          col="lightgrey", lty=0) # Plot 95% credible interval of posterior predictive distribution as a grey shaded area
  points(raw_data[,var_x], raw_data[,var_y], pch=20) # Plot raw data points
  curve(one_compartment(pars_hsk[1],pars_hsk[2],pars_hsk[3],C_exp,x,0.5*x_max), add=T) # Plot model's predicted values as a curve again
  
} # End Plot_PostPredDistr()

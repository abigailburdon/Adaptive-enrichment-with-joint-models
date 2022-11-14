library(rstudioapi)
library(parallel)
library(truncnorm)
library(JM)
library(numDeriv)
library(spatstat)
library(mvtnorm)
library(nlme)
setwd(dirname(getActiveDocumentContext()$path ))
source("RMST_functions.R")
source("numerical_integration.R")
source("error_spending_threshold.R")
source("generate_patients.R")
source("conditional_score.R")
source("all_models.R")
source("enrichment_trial.R")

N <- 10 #Number of simulations
gamma <- 0.8
sigma_sq <- 1
var2 <- 5
alpha <- 0.025
beta <- 0.1
K <- 2
eta1 <- -0.5
eta2 <- 0
b2 <- c(-0.25,0)
lambda <- 2/3
b_mean = c(4.23, 1.81)
n <- 1e3
obs_times <- c(seq(0,0.25,by = 2/52), seq(12/52, 5, by = 1/12))
b_sigma <- matrix(c(2.5,1.7,1.7,var2), 2,2)
knot_point <- 1
t_star <- 5
p1 <- 0.6

set.seed(as.integer(12345))

cl <- makeCluster(16)  
invisible(clusterEvalQ(cl,{
  library(JM)
  library(parallel)
  library(truncnorm)
  library(JM)
  library(numDeriv)
  library(spatstat)
  library(mvtnorm)
  library(nlme)
  source("RMST_functions.R")
  source("numerical_integration.R")
  source("error_spending_threshold.R")
  source("generate_patients.R")
  source("conditional_score.R")
  source("all_models.R")
  source("enrichment_trial.R")
}))
clusterExport(cl, c("gamma","eta1", "eta2", "sigma_sq","lambda",
                    "t_star", "knot_point", "b_sigma", "b2", "alpha", "beta"))

#########################
### Design parameters ###
#########################
#Number of events at each analysis and maximum information
start_time <- Sys.time()
{
#Threshold selection criteria - required information at the first interim
#analysis in subgroup S1
root_find <- optim(par = c(4, 0.5), function(x){
  I1 <- x[1]
  C_val <- x[2]
  P1 <- (1-pnorm(C_val+eta1*sqrt(I1)))*pnorm(C_val)
  P2 <- pnorm(C_val+eta1*sqrt(I1))*(1-pnorm(C_val))
  PF <- (1-pnorm(C_val+eta1*sqrt(I1)))*(1-pnorm(C_val))
  P_empty <- 1-(P1+P2+PF)
  val1 <- PF-P_empty
  val2 <- P1-p1
  return(val1^2+val2^2)
})
C_val <- root_find$par[2]
I1_required <- root_find$par[1]

#Vector of true parameters in the joint model - give the same names as the MLE
#from the package JM
theta_true <- c("betas.(Intercept)"=b_mean[1],
               "betas.obs_time"=b_mean[2],
               "betas.obs_time:treatment"=b2[1],
               "sigma"=log(sqrt(sigma_sq)),
               "gammas.treatment"=eta1,
               "alpha.alpha"=gamma,
               "D1"=log(sqrt(b_sigma[1,1])),
               "D2"=b_sigma[2,1]/sqrt(b_sigma[1,1]),
               "D4"=log(sqrt(b_sigma[2,2])))

#Find baseline hazard as the value of h0 such that median survival is 1 year
#for patients in control group (alternatively under H0)
h0_val <- uniroot(function(h0){
  h0_val <- h0
  theta_true["xi.xi.1"] <- log(1/h0_val)
  theta_true["xi.xi.2"] <- log(1/(h0_val*0.6))
  theta_true["gammas.treatment"] <- 0
  theta_true["betas.obs_time:treatment"] <- 0
  return(S.b0b1(theta_true, knot_point, 1)-0.5)},
  lower = 0.001, upper = 1e5)$root
theta_true["xi.xi.1"] <- log(1/h0_val)
theta_true["xi.xi.2"] <- log(1/(h0_val*0.6))

##Relationship between number of events and information
#Generate data
all.data <- generate_ps(1e3, gamma_v = gamma, eta_f = c(eta1, eta2),
                        sigma_sq = sigma_sq, h0_val=h0_val,
                        recruit_end =2, b2 = b2, b_sigma=b_sigma,
                        lambda = lambda, knot_point = knot_point)   
surv.data <- all.data[[1]]
long.data <- all.data[[2]]
#Calendar times of events in subgroup S1
t_vals1 <- surv.data[surv.data$status==1&surv.data$group==1,"surv_times"] + 
  surv.data[surv.data$status==1&surv.data$group==1,"arrival_times"] 
t_vals1 <- t_vals1[order(t_vals1)]
#Calendar times of events in subgroup S2
t_vals2 <- surv.data[surv.data$status==1&surv.data$group==2,"surv_times"] + 
  surv.data[surv.data$status==1&surv.data$group==2,"arrival_times"] 
t_vals2 <- t_vals2[order(t_vals2)]
#Start from 10th event to avoid complications with noise information
event.min <- 10


### Specific to conditional score method
#Information at each event time in subgroup 1
clusterExport(cl, c("surv.data", "long.data"))
I_vals1 <- parSapply(cl, t_vals1[seq(event.min, 100, length=20)], function(t_val){
  
  data_t <- data_cs(surv.data, long.data, t_val)
  surv_1 <- subset(data_t[[1]], group ==1)
  long_1 <- subset(data_t[[2]], group ==1)
  I1 <- treat_effect1(surv_1, long_1, eta1, t_star, knot_point)[2]
  return(I1)
  
})
I_vals2 <- parSapply(cl, t_vals2[seq(event.min, 100, length=20)], function(t_val){
  
  data_t <- data_cs(surv.data, long.data, t_val)
  surv_2 <- subset(data_t[[1]], group ==2)
  long_2 <- subset(data_t[[2]], group ==2)
  I2 <- treat_effect1(surv_2, long_2, eta2, t_star, knot_point)[2]
  return(I2)
  
})

#Constants m_j describing relationship between number of events and information
m_vals <- rep(0, 2)
lm1 <- lm(I~n_event-1,data = data.frame(
  I=I_vals1,n_event = seq(event.min, 100, length=20)))
lm2 <- lm(I~n_event-1,data = data.frame(
  I=I_vals2,n_event = seq(event.min, 100, length=20)))
m_vals[1] <- 1/lm1$coefficients
m_vals[2] <- 1/lm2$coefficients

#Number of events required at each analysis
d11 <- I1_required*m_vals[1]
events1 <- ceiling(d11)
d2 <- n2_calc(K,alpha,beta,delta01=eta1,delta02=eta2,C_val,
              lambda,1/m_vals[1],1/m_vals[2],I1_val=I1_required)
events2 <- ceiling(d2)

#Maximum information needed at final analysis
I_max <- d2/(lambda/m_vals[1]+(1-lambda)*m_vals[2])
  
###Repeat using RMST analysis - this has different I_max because the 
#methods have different endpoints
#For RMST, must have 2 events after knot point. Find all event times for this subset of patients
t_vals1 <- surv.data[surv.data$status==1&surv.data$group==1&surv.data$surv_times>1,"arrival_times"] +
  surv.data[surv.data$status==1&surv.data$group==1&surv.data$surv_times>1, "surv_times"]
t_vals2 <- surv.data[surv.data$status==1&surv.data$group==2&surv.data$surv_times>1,"arrival_times"] +
  surv.data[surv.data$status==1&surv.data$group==2&surv.data$surv_times>1, "surv_times"]
n_events1 <- sapply(t_vals1[1:20], function(t_val){
  sum(surv.data[surv.data$group==1&surv.data$status==1, "surv_times"] + 
        surv.data[surv.data$group==1&surv.data$status==1, "arrival_times"] <= t_val)
})
n_events2 <- sapply(t_vals2[1:20], function(t_val){
  sum(surv.data[surv.data$group==2&surv.data$status==1, "surv_times"] + 
        surv.data[surv.data$group==2&surv.data$status==1, "arrival_times"] <= t_val)
})

#Information at each event time in subgroup 1
I_vals1 <- parSapply(cl, t_vals1[1:20], function(t_val){
  
  data_t <- data_cs(surv.data, long.data, t_val)
  surv_1 <- subset(data_t[[1]], group ==1)
  long_1 <- subset(data_t[[2]], group ==1)
  I1 <- treat_effect2(surv_1, long_1, eta1, t_star, knot_point)[2]
  return(I1)
  
})
I_vals2 <- parSapply(cl, t_vals2[1:20], function(t_val){
  
  data_t <- data_cs(surv.data, long.data, t_val)
  surv_2 <- subset(data_t[[1]], group ==2)
  long_2 <- subset(data_t[[2]], group ==2)
  I2 <- treat_effect1(surv_2, long_2, eta2, t_star, knot_point)[2]
  return(I2)
  
})

#Constants m_j describing relationship between number of events and information
m_vals <- rep(0, 2)
lm1 <- lm(I~n_event-1,data = data.frame(
  I=I_vals1,n_event = n_events1))
lm2 <- lm(I~n_event-1,data = data.frame(
  I=I_vals2,n_event = n_events2))
m_vals[1] <- 1/lm1$coefficients
m_vals[2] <- 1/lm2$coefficients

#RMST objects at true parameter values
RMST_true1 <- RMST_treat(theta_true, t_star, knot_point)
RMST_true0 <- RMST_control(theta_true, t_star, knot_point)
RMST_true <- RMST_true1-RMST_true0

#Threshold selection criteria - required information at the first interim
#analysis in subgroup S1
root_find <- optim(par = c(4, 0.5), function(x){
  I1 <- x[1]
  C_val <- x[2]
  P1 <- (1-pnorm(C_val-RMST_true*sqrt(I1)))*pnorm(C_val)
  P2 <- pnorm(C_val-RMST_true*sqrt(I1))*(1-pnorm(C_val))
  PF <- (1-pnorm(C_val-RMST_true*sqrt(I1)))*(1-pnorm(C_val))
  P_empty <- 1-(P1+P2+PF)
  val1 <- PF-P_empty
  val2 <- P1-p1
  return(val1^2+val2^2)
})
C_val.RMST <- root_find$par[2]
I1_required.RMST <- root_find$par[1]

#Maximum information needed at final analysis
d2.RMST <- n2_calc(K,alpha,beta,delta01=RMST_true,delta02=0,C_val.RMST,
              lambda,1/m_vals[1],1/m_vals[2],I1_val=I1_required.RMST)
I_max.RMST <- d2.RMST/(lambda/m_vals[1]+(1-lambda)*m_vals[2])
}
end_time <- Sys.time()
paste("Calculation of design parameters takes", end_time-start_time)

########################
### Simulation study ###
########################
all.seeds <- 1e9*runif(N, -1, 1)
clusterExport(cl, c("all.seeds", "h0_val", "events1", "events2",
                    "RMST_true", "C_val", "C_val.RMST", "I_max", "I_max.RMST"))
start_time <- Sys.time()
vals <- parSapply(cl, 1:N, function(sim.rep){
  
  #For each method, set the seed so that the same dataset is generated
  
  #Conditional score
  set.seed(as.integer(all.seeds[sim.rep]))
  trial.1 <- run_trial(treat_effect1, t_star, knot_point, 1e3, gamma, eta1, eta2,
                       -eta1, 0, sigma_sq, h0_val, b2, b_sigma, lambda,
                       events1, events2, C_val, I_max)
  
  #RMST
  set.seed(as.integer(all.seeds[sim.rep]))
  trial.2 <- run_trial(treat_effect2, t_star, knot_point, 1e3, gamma, eta1, eta2,
                       RMST_true, 0, sigma_sq, h0_val, b2, b_sigma, lambda,
                       events1, events2, C_val.RMST, I_max.RMST)
  
  #Cox model withough longitudinal data
  set.seed(as.integer(all.seeds[sim.rep]))
  trial.3 <- run_trial(treat_effect3, t_star, knot_point, 1e3, gamma, eta1, eta2,
                       -eta1, 0, sigma_sq, h0_val, b2, b_sigma, lambda,
                       events1, events2, C_val, I_max)
  
  #Cox model using longitudinal data as time-varying covariate
  set.seed(as.integer(all.seeds[sim.rep]))
  trial.4 <- run_trial(treat_effect4, t_star, knot_point, 1e3, gamma, eta1, eta2,
                       -eta1, 0, sigma_sq, h0_val, b2, b_sigma, lambda,
                       events1, events2, C_val, I_max)
  
  return(c(trial.1, trial.2, trial.3, trial.4))
  
})
end_time <- Sys.time()
paste("Simulation stuy of", N, "replicates takes", end_time-start_time)


### Power
rowMeans(as.matrix(vals[c(5, 31, 57, 83),]))

         
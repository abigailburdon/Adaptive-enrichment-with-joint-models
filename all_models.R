###Functions that take a dataset and return treatment effect estimate and information

#Although inputs t_star and knot_point are only relevant to method 2 (RMST),
#include them so that a generic treatment effect function can be called.

#Conditional score method
treat_effect1 <- function(surv.data, long.data, eta, t_star, knot_point){
  
  #Relabel patients
  id.data <- data.frame(old.id = surv.data$patient, new.id = 1:nrow(surv.data))
  surv.data$patient <- id.data$new.id
  long.data$patient <- sapply(long.data$patient, function(x){
    id.data[id.data$old.id == x, "new.id"]
  })
  
  
  # Conditional score items
  cs_lists <- cs_indep(surv.data, long.data)
  t_vals <- cs_lists[[1]]
  rs <- cs_lists[[2]]
  theta_vals <- cs_lists[[3]]
  X_vals <- cs_lists[[4]]
  sigma_hat <- sigma_guess(surv.data, long.data, nrow(surv.data),
                           t_vals, rs, theta_vals, X_vals)
  
  #Boundaries of optimisation function
  if(gamma==0){
    min_vals <- c(-2.5,-2,0)
  }else{
    min_vals <- c(0,-2,0)
  }
  max_vals <- c(2, 2, Inf)
  
  #Find root until within boundaries
  count <- 1
  root_find <- optim(
    par = c(gamma,eta,sigma_sq),
    fn = function(x) {
      out <- cond_score(x[1],x[2],x[3],surv.data,long.data,
                        t_vals,rs,theta_vals,X_vals)
      out[1]^2+out[2]^2},
    method = "L-BFGS-B", lower = min_vals, upper = max_vals)
  
  while((root_find$par[1] == min_vals[1] | root_find$par[1] == max_vals[1] |
         root_find$par[2] == min_vals[2] | root_find$par[2] == 2) &
        count < 20){
    root_find <- optim(
      par = runif(3, min_vals, c(max_vals[1:2], 20)),
      fn = function(x) {
        out <- cond_score(x[1],x[2],x[3],surv.data,long.data,
                          t_vals,rs,theta_vals,X_vals)
        out[1]^2+out[2]^2},
      method = "L-BFGS-B", lower = min_vals, upper = max_vals)
    
    count <- count +1
  }
  
  # Treatment effect estimate
  theta_hat <- c(root_find$par[1],root_find$par[2], root_find$par[3])
  
  
  # Sandwich variance 
  AB <- sand_var(theta_hat[1], theta_hat[2], theta_hat[3],
                 surv.data, long.data,
                 t_vals, rs, theta_vals, X_vals)
  A <- AB[[1]]
  B <- AB[[2]]
  #A_inv <- matrix(c(A[2,2],-A[1,2],-A[2,1],A[1,1]),2,2)/(A[1,1]*A[2,2]-A[2,1]*A[1,2])
  V_mat <- B[2,2]/(A[2,2]^2)
  I <- nrow(surv.data)^2/(nrow(t_vals)*V_mat)
  
  #return(c(theta_hat[2],  I*exp(-0.27*gamma*sigma_sq)))
  return(c(theta_hat[2],  I))
  
}
#Restricted mean survival time
treat_effect2 <- function(surv.data, long.data, eta, t_star, knot_point){
  
  
  
  #Relabel patients
  id.data <- data.frame(old.id = surv.data$patient, new.id = 1:nrow(surv.data))
  surv.data$patient <- id.data$new.id
  long.data$patient <- sapply(long.data$patient, function(x){
    id.data[id.data$old.id == x, "new.id"]
  })
  
  long.model <- lme(long_obs~ obs_time*treatment - treatment,
                    random =~ obs_time|patient,
                    data = long.data,
                    control = lmeControl(opt = "optim"))
  surv.model <- coxph(Surv(surv_times, status)~treatment,
                      data = surv.data, x=T)
  joint.model <- jointModel(long.model, surv.model, "obs_time",
                            method = "piecewise-PH-aGH",
                            control = list(knots=1))
  
  #Parameter estimates
  theta_hat <- unlist(joint.model$coefficients)
  theta_hat <- theta_hat[c(1:10,12)]
  theta_hat[c(4,7,8)] <- log(theta_hat[c(4,7,8)])
  theta_hat[10] <- theta_hat[10]/sqrt(theta_hat[9])
  theta_hat[c(9,11)] <- log(sqrt(theta_hat[c(9,11)]))
  
  #Calculate RMST
  RMST_hat1 <- RMST_treat(theta_hat, t_star, knot_point)
  RMST_hat0 <- RMST_control(theta_hat, t_star, knot_point)
  RMST_hat <- RMST_hat0-RMST_hat1
  
  #Variance of RMST
  d1 <- grad(function(theta_hat){RMST_treat(theta_hat, t_star, knot_point)}, theta_hat)
  d0 <- grad(function(theta_hat){RMST_control(theta_hat, t_star, knot_point)}, theta_hat)    
  RMST_var <- t(d1-d0)%*%vcov(joint.model)%*%(d1-d0)
  
  I <- 1/RMST_var
  if(is.na(I)|I<0){
    I <- 0
  }
  
  return(c(RMST_hat,  I))
}
##Cox model without longitudinal data
treat_effect3 <- function(surv.data, long.data, eta, t_star, knot_point){
  
  max_t <- max(surv.data[surv.data$status==1,"surv_times"])
  if(max_t <= 2){
    RMST_hat <- 0
    I <- 0
  }else{
    time   = surv.data$surv_times
    status = surv.data$status
    arm    = surv.data$treatment
    m.RMST <- rmst2(time, status, arm, tau=2)
    
    RMST_hat <- m.RMST$RMST.arm0$rmst[1]-m.RMST$RMST.arm1$rmst[1]
    I <- 1/(m.RMST$RMST.arm1$rmst[2]^2+m.RMST$RMST.arm0$rmst[2]^2)
  }
  
  return(c(RMST_hat,  I))
}
#Cox model with longitudinal data as time-varying covariate
treat_effect4 <- function(surv.data, long.data, eta, t_star, knot_point){
  
  #Relabel patients
  id.data <- data.frame(old.id = surv.data$patient, new.id = 1:nrow(surv.data))
  surv.data$patient <- id.data$new.id
  long.data$patient <- sapply(long.data$patient, function(x){
    id.data[id.data$old.id == x, "new.id"]
  })
  
  long.model <- lme(long_obs~ obs_time*treatment - treatment,
                    random =~ obs_time|patient,
                    data = long.data,
                    control = lmeControl(opt = "optim"))
  surv.model <- coxph(Surv(surv_times, status)~treatment,
                      data = surv.data, x=T)
  joint.model <- jointModel(long.model, surv.model, "obs_time",
                            method = "piecewise-PH-aGH",
                            control = list(knots=1))
  
  #Parameter estimates
  theta_hat <- unlist(joint.model$coefficients)
  theta_hat <- theta_hat[c(1:10,12)]
  theta_hat[c(4,7,8)] <- log(theta_hat[c(4,7,8)])
  theta_hat[10] <- theta_hat[10]/sqrt(theta_hat[9])
  theta_hat[c(9,11)] <- log(sqrt(theta_hat[c(9,11)]))
  
  #Calculate RMST
  RMST_hat1 <- RMST_treat(theta_hat, 2, 1)
  RMST_hat0 <- RMST_control(theta_hat, 2, 1)
  RMST_hat <- RMST_hat1-RMST_hat0
  
  #Variance of RMST
  #d1 <- grad(function(theta_hat){RMST_treat(theta_hat, 2, 1)}, theta_hat)
  #d0 <- grad(function(theta_hat){RMST_control(theta_hat, 2, 1)}, theta_hat)    
  RMST_var <- t(d1.2-d0.2)%*%vcov(joint.model)%*%(d1.2-d0.2)
  
  I <- 1/RMST_var
  if(is.na(I)|I<0){
    I <- 0
  }
  
  return(c(RMST_hat,  I))
}
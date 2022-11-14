#Survival function
S <- function(theta_hat, knot_point, u, y0, y1){
  
  #Longitudinal coefficient
  gamma <- theta_hat["alpha.alpha"]
  #Treatment effect
  eta <- theta_hat["gammas.treatment"]
  #Baseline hazard
  h0_vals  <- exp(-theta_hat[c("xi.xi.1","xi.xi.2")])
  #Random effects
  b0_mean <- theta_hat["betas.(Intercept)"]
  b1_mean <- theta_hat["betas.obs_time"]
  b2 <- theta_hat["betas.obs_time:treatment"]
  #Variance of random effects
  b0_sigma <- exp(theta_hat["D1"])
  b1_sigma <- exp(theta_hat["D4"])
  rho <- theta_hat["D2"]*b0_sigma
  
  
  #For Gauss-Hermit integration in 2D, transform correlated 
  #bivariate normal to get independent normals
  b0 <- y0
  b1 <- -rho*(y1-y0)/(b0_sigma^2)
  
  #Cummulative hazard
  if(gamma==0&u<=knot_point){
    H_val <- u*exp(eta)/h0_vals[1]
  }
  if(gamma==0&u>knot_point){
    int1 <- knot_point*exp(eta)/h0_vals[1]
    int2 <- (u-knot_point)*exp(eta)/h0_vals[2]
    H_val <- int1+int2
  }
  if(gamma!=0&u<=knot_point){
    v1 <- exp(gamma*b0+eta)
    v2 <- exp(gamma*u*(b1+b2))-1
    v3 <- h0_vals[1]*gamma*(b1+b2)
    H_val <- v1*v2/v3 
  }
  if(gamma!=0&u>knot_point){
    v1 <- exp(gamma*b0+eta)
    v2 <- exp(gamma*knot_point*(b1+b2))-1
    v3 <- h0_vals[1]*gamma*(b1+b2)
    int1 <- v1*v2/v3 
    v1 <- exp(gamma*b0+eta)
    v2 <- exp(gamma*u*(b1+b2))-exp(gamma*knot_point*(b1+b2))
    v3 <- h0_vals[2]*gamma*(b1+b2)
    int2 <- v1*v2/v3 
    H_val <- int1+int2
  }
  
  #Survival function
  S_val <- exp(-H_val)
  
  return(S_val)
}
S_Vec <- Vectorize(S, vectorize.args = "u")

#Survival function integrated over random effects
S.b0b1 <- function(theta_hat, knot_point, u){
  
  #Longitudinal coefficient
  gamma <- theta_hat["alpha.alpha"]
  #Treatment effect
  eta <- theta_hat["gammas.treatment"]
  #Baseline hazard
  h0_vals  <- exp(-theta_hat[c("xi.xi.1","xi.xi.2")])
  #Random effects
  b0_mean <- theta_hat["betas.(Intercept)"]
  b1_mean <- theta_hat["betas.obs_time"]
  b2 <- theta_hat["betas.obs_time:treatment"]
  #Variance of random effects
  b0_sigma <- exp(theta_hat["D1"])
  b1_sigma <- exp(theta_hat["D4"])
  rho <- theta_hat["D2"]*b0_sigma
  
  #Integrate over b0
  int.b0 <- function(theta_hat,knot_point,u,y1){
    gauss.hermite(function(y0){S(theta_hat,knot_point,u,y0,y1)},
                  b_mean[1],sqrt(b_sigma[1,1]),order=8)
  }
  int.b0b1 <- function(theta_hat,knot_point,u){
    gauss.hermite(function(y1){int.b0(theta_hat,knot_point,u,y1)},
                  b_mean[1]-b_sigma[1,1]*b_mean[2]/b_sigma[1,2],
                  sqrt(b_sigma[1,1])*sqrt((sqrt(b_sigma[1,1])*sqrt(b_sigma[2,2])/b_sigma[1,2])^2-1),order=8)
  }
  return(int.b0b1(theta_hat,knot_point,u))
}

#RMST for treatment group
RMST_treat <- function(theta_hat, t_star, knot_point){
  
  
  #Longitudinal coefficient
  gamma <- theta_hat["alpha.alpha"]
  #Treatment effect
  eta <- theta_hat["gammas.treatment"]
  #Baseline hazard
  h0_vals  <- exp(-theta_hat[c("xi.xi.1","xi.xi.2")])
  #Random effects
  b0_mean <- theta_hat["betas.(Intercept)"]
  b1_mean <- theta_hat["betas.obs_time"]
  b2 <- theta_hat["betas.obs_time:treatment"]
  #Variance of random effects
  b0_sigma <- exp(theta_hat["D1"])
  b1_sigma <- exp(theta_hat["D4"])
  rho <- theta_hat["D2"]*b0_sigma
  
  
  #RMST with b0 and b1 fixed
  RMST <- function(theta_hat, knot_point, y0, y1){
    integrate(function(u){S_Vec(theta_hat, knot_point, u, y0, y1)}, lower = 0, upper = t_star)$value
  }
  
  #Integrate over b0
  int.b0 <- function(theta_hat, knot_point, y1){
    gauss.hermite(function(y0){RMST(theta_hat,knot_point,y0,y1)},
                  b0_mean,b0_sigma,order=8)
  }
  
  RMST_val <- gauss.hermite(function(y1){int.b0(theta_hat, knot_point, y1)},
                            b0_mean-b0_sigma^2*b1_mean/rho,
                            b0_sigma*sqrt((b0_sigma*b1_sigma/rho)^2-1),order=8)
  
  return(RMST_val)
}
#RMST for control group
RMST_control <- function(theta_hat, t_star, knot_point){
  

  #Set parameters describing treatment effect to zero
  theta_hat["gammas.treatment"] <- 0
  theta_hat["betas.obs_time:treatment"] <- 0
  
  #Longitudinal coefficient
  gamma <- theta_hat["alpha.alpha"]
  #Treatment effect
  eta <- theta_hat["gammas.treatment"]
  #Baseline hazard
  h0_vals  <- exp(-theta_hat[c("xi.xi.1","xi.xi.2")])
  #Random effects
  b0_mean <- theta_hat["betas.(Intercept)"]
  b1_mean <- theta_hat["betas.obs_time"]
  b2 <- theta_hat["betas.obs_time:treatment"]
  #Variance of random effects
  b0_sigma <- exp(theta_hat["D1"])
  b1_sigma <- exp(theta_hat["D4"])
  rho <- theta_hat["D2"]*b0_sigma
  
  
  #RMST with b0 and b1 fixed
  RMST <- function(theta_hat, knot_point, y0, y1){
    integrate(function(u){S_Vec(theta_hat, knot_point, u, y0, y1)}, lower = 0, upper = t_star)$value
  }
  
  #Integrate over b0
  int.b0 <- function(theta_hat, knot_point, y1){
    gauss.hermite(function(y0){RMST(theta_hat,knot_point,y0,y1)},
                  b0_mean,b0_sigma,order=8)
  }
  
  RMST_val <- gauss.hermite(function(y1){int.b0(theta_hat, knot_point, y1)},
                            b0_mean-b0_sigma^2*b1_mean/rho,
                            b0_sigma*sqrt((b0_sigma*b1_sigma/rho)^2-1),order=8)
  
  return(RMST_val)
}

#Simulate entire data set of patients
generate_ps <- function(n,
                        gamma_v=0.8,
                        eta_f=c(-0.25, 0),
                        lambda_censor=5e-5,
                        sigma_sq=3.4,
                        h0_val = 118,
                        obs_times = c(seq(0,0.23,by = 2/52), seq(12/52, 2, by = 1/12)),
                        recruit_end = 2,
                        b2 = c(-0.5,0),
                        b_mean = c(4.23, 1.81),
                        b_sigma = matrix(c(2.5,1.7,1.7,5),2,2),
                        lambda = 2/3,
                        knot_point)


{
  
  #Group 
  S <- rbinom(n, 1, lambda)
  S[S==0] <- 2
  
  #fixed covariates
  Z <- rbinom(n, 1, prob = 0.5)
  
  ###Variables to fix
  b_mat <- rmvnorm(n, b_mean, b_sigma)
  b_mat[,2] <- b_mat[,2]+b2[S]*Z
  
  #Generate survival times for each patient
  knot_point <- 1
  h0_val <- h0_val*c(1, 0.6)
  T_val <- sapply(1:n, function(i){
    
    H <- function(t){
      if(gamma_v==0&t<=knot_point){
        H_val <- t*exp(eta_f[S[i]]*Z[i])/h0_val[1]
      }
      if(gamma_v==0&t>knot_point){
        H_val <- exp(eta_f[S[i]]*Z[i])*(knot_point*(1/h0_val[1]-1/h0_val[2])+t/h0_val[2])
      }
      if(gamma_v!=0&t<=knot_point){
        v1 <- exp(gamma_v*b_mat[i,1]+eta_f[S[i]]*Z[i])
        v2 <- exp(gamma_v*t*b_mat[i,2])-1
        v3 <- h0_val[1]*gamma_v*b_mat[i,2]
        H_val <- v1*v2/v3 
      }
      if(gamma_v!=0&t>knot_point){
        v1 <- exp(gamma_v*b_mat[i,1]+eta_f[S[i]]*Z[i])/(gamma_v*b_mat[i,2])
        v2 <- (1/h0_val[1]-1/h0_val[2])*exp(gamma_v*knot_point*b_mat[i,2])
        v3 <- exp(gamma_v*t*b_mat[i,2])/h0_val[2]-1/h0_val[1]
        H_val <- v1*(v2+v3)
      }
      return(H_val)
    }
    
    #### inversion method e.g T = H^-1(u) ##########
    
    #generate uniform random variable for use in the inversion method
    u <- runif(1,0,1)
    
    if(H(10*obs_times[length(obs_times)]) < -log(1-u)){
      Ts <- Inf
    }else{
      Ts <- uniroot(function(t) H(t) + log(1-u), lower=0,
                    upper = 2*obs_times[length(obs_times)],
                    extendInt = "upX")$root
    }
    
    return(Ts)})
  
  ####Censoring times ######
  C_val <- rexp(n, lambda_censor)
  status <-  as.numeric(as.numeric(T_val) < as.numeric(C_val))
  T_val[status==0] <- C_val[status==0] #replace event times with censoring times for
  T_val <- as.numeric(T_val)           #censored observations
  
  #Generate arrival times
  arrival_times <- runif(n, 0, recruit_end)
  
  #Order to help with large samples
  gen_ordered <- order(T_val+arrival_times)
  
  #All survival data
  surv.data <- data.frame(patient=c(1:n),
                          treatment=Z[gen_ordered],
                          surv_times = T_val[gen_ordered],
                          status = status[gen_ordered],
                          arrival_times = arrival_times[gen_ordered],
                          group = S[gen_ordered])
  
  
  
  #Generate longitudinal data for each patient
  long.data <- data.frame(patient = NULL, time_obs = NULL, long_obs = NULL)
  for(i in 1:n){
    
    #generate a random number of longitudinal measurement times
    obs_time <- obs_times[which(obs_times<=surv.data$surv_times[i])]
    num_obs <- length(obs_time)
    
    
    #at each measurement time find true value of longitudinal observation, then 
    #add measurement error
    long_obs <- c(b_mat[gen_ordered[i],]%*%t(matrix(c(rep(1, num_obs),obs_time), ncol = 2)))+
      rnorm(num_obs, 0, sqrt(sigma_sq))
    
    long.data <- rbind(long.data, data.frame(patient = rep(surv.data$patient[i], num_obs),
                                             obs_time,long_obs,
                                             group = rep(surv.data$group[i], num_obs),
                                             treatment = rep(surv.data$treatment[i], num_obs)))}
  
  return(list(surv.data, long.data))}

#Function to truncate both survival and longitudinal data sets at time time_k
data_cs <- function(surv.data, long.data, time_k){
  
  
  # Survival observations
  surv.data_k <- surv.data[surv.data$arrival_times <= time_k,]
  #Censor unobserved observations and change time
  un_obs_k <- surv.data_k$surv_times+surv.data_k$arrival_times > time_k
  surv.data_k$surv_times[un_obs_k] <- time_k-surv.data_k$arrival_times[un_obs_k]
  surv.data_k$status[un_obs_k] <- 0
  
  #Longitudinal observations
  ld <- sapply(1:nrow(surv.data_k), function(i){
    t_arrive <- surv.data_k$arrival_times[i]
    long_all_i <- long.data[long.data$patient == surv.data_k$patient[i],]
    long_val <- long_all_i[long_all_i$obs_time+t_arrive<=time_k,]
    long_val$patient <- rep(surv.data_k$patient[i], nrow(long_val))
    long_val
  }, simplify = F)
  long.data_k <- do.call(rbind, ld)
  
  return(list(surv.data_k, long.data_k))
}
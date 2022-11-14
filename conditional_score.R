#Functions needed to calculate treatment effect estimate and information levels
#when using the conditional score analysis method

cs_indep <- function(surv.data, long.data){
  #Some necessary parts of the conditional score function do not involve gamma, eta
  #and sigma - for root finding, these functions only need to be evaluated once, not
  #for varying values of gamma, eta, sigma.
  n <- nrow(surv.data) #number of patients
  
  
  
  #Find all event times which will contribute to the conditional score sum
  t_vals <- sapply(1:n, function(i){
    #Event times only contribute if the observation is exact and there are enough
    #longitudinal measurements before t to do regression
    no_censoring <- surv.data$status[i]==1
    all_long_i <- long.data[long.data$patient==i,]
    enough_measures <- nrow(all_long_i)>=2
    return(no_censoring&enough_measures)})
  
  t_vals <- data.frame(patient=which(t_vals),
                       time = surv.data[which(t_vals),"surv_times"])
  
  if(nrow(t_vals)==0){
    rs <- list(NULL)
    theta_vals <- list(NULL)
    X_vals <- list(NULL)
  }else{
    
    
    #A function returning all patients at risk at a particular time
    #A patient is at risk at time t if they are still alive and also if they have
    #enough longitudinal measurements before t to do regression
    risk_set <- function(t, surv.data, long.data){
      ds <- surv.data[surv.data$surv_times>=t,"patient"]
      ds_obs <- table(factor(long.data[long.data$obs_time <=t,"patient"], lev = ds))
      return(ds[ds_obs>=2])
    }
    rs <- sapply(1:nrow(t_vals), function(j){
      risk_set(t_vals[j, "time"], surv.data, long.data)
    }, simplify = F)
    
    
    #variance of the predicted observation
    theta <- function(i, t, surv.data, long.data){
      #function takes as inputs: patient number and time,
      #returns variance of the predicted longitudinal measurement
      ss <- long.data[long.data$patient == i, ]
      ts <- ss[ss$obs_time <=t , "obs_time"]
      1/length(ts) + ((t-mean(ts))^2) / sum((ts-mean(ts))^2)
    }
    theta_Vec <- Vectorize(theta, vectorize.args = "i")
    theta_vals <- sapply(1:nrow(t_vals), function(j){
      theta_Vec(rs[[j]], t_vals[j, "time"], surv.data, long.data)
    }, simplify = F)
    
    X_hat <- function(i, t, long.data){
      #function takes as inputs: patient number and time,
      #returns value of predicted longitudinal  measurement
      set1 <- long.data[long.data$patient==i,]
      set2 <- set1[set1$obs_time<=t,]
      x_obs <- set2$obs_time
      y_obs <- set2$long_obs
      
      x_bar <- mean(x_obs)
      y_bar <- mean(y_obs)
      
      #ordinary least squares
      b1_hat <- sum((x_obs-x_bar)*(y_obs-y_bar))/sum(((x_obs-x_bar))^2)
      b0_hat <- y_bar-b1_hat*x_bar
      
      return(b0_hat+b1_hat*t)}
    
    X_vals <- list(NULL)
    for(j in 1:nrow(t_vals)){
      X_vals[[j]] <- sapply(rs[[j]], function(i) X_hat(i, t_vals[j, "time"], long.data))
    }
  }
  
  
  
  #return all values that will be used again in the part dependent on gamma, eta
  return(list(t_vals, rs, theta_vals, X_vals))
}

sigma_guess <- function(surv.data, long.data, n,
                        t_vals, rs, theta_vals, X_vals){
  
  #Residual sum of squares needed for estimating sigma_sq
  RSS <- function(i, long.data){
    #Function resturning the resisdual sum of squares for patient i
    set <- subset(long.data, patient==i)
    x_obs <- set$obs_time  #some cases result in error when all observations are zero
    
    if(nrow(set)>2 & sum(x_obs)!=0){
      #Patient i only contributes if they have more than two longitudinal measurement
      #this is for regression
      y_obs <- set$long_obs
      
      x_bar <- mean(x_obs)
      y_bar <- mean(y_obs)
      
      b1_hat <- sum((x_obs-x_bar)*(y_obs-y_bar))/sum(((x_obs-x_bar))^2)
      b0_hat <- y_bar-b1_hat*x_bar
      
      num <- sum((y_obs-b0_hat-b1_hat*x_obs)^2)
      denom <- length(y_obs)-2
    }else{
      num <- 0
      denom <- 0
    }
    return(c(num, denom))
  }
  rss <- rowSums(sapply(1:n, function(i) RSS(i, long.data)))
  
  return(rss[1]/rss[2])
  
}

cond_score <- function(gamma, eta, sigma_sq,
                       surv.data, long.data,
                       t_vals, rs, theta_vals, X_vals){
  #Define empty lists
  S <- list(NULL)
  S_i <- rep(0, nrow(t_vals))
  E0 <- list(NULL)
  Z <- list(NULL)
  Z_i <- rep(0, nrow(t_vals))
  
  if(nrow(t_vals) ==0){
    vals <- matrix(0, 1, 2)
  }else{
    for(j in 1:nrow(t_vals)){
      
      #keep a marker for the patient who's event time we're interested in
      #where they lie in the risk set
      patient_marker <- which(rs[[j]]== t_vals[j,"patient"])
      
      
      
      #sufficient statistic
      S[[j]] <- X_vals[[j]]
      S_i[j] <- S[[j]][patient_marker]+gamma*sigma_sq*theta_vals[[j]][patient_marker]
      S[[j]][patient_marker] <- S_i[j]
      
      
      #covariates
      Z[[j]] <- surv.data[surv.data$patient%in%rs[[j]],"treatment"]
      Z_i[j] <- Z[[j]][patient_marker]
      
      
      #condtitional density
      E0[[j]] <- exp(gamma*S[[j]]-gamma^2*sigma_sq*theta_vals[[j]]/2 +eta*Z[[j]])
      
      
    }
    
    
    #Conditional score function - create vector where each element is the contribution to
    #the conditional score function at an event time
    
    #Vector E1 at each event time
    E1_1 <- mapply(function(x,y) x*y, S, E0)
    E1_2 <- mapply(function(x,y) x*y, Z, E0)
    
    #contribution at event time T_i
    vals <- matrix(0, nrow = nrow(t_vals), ncol = 2)
    for(i in 1:nrow(t_vals)){
      vals[i,1] <- S_i[i]-sum(E1_1[[i]])/sum(E0[[i]])
      vals[i,2] <- Z_i[i]-sum(E1_2[[i]])/sum(E0[[i]])
    }
    
  }
  
  return(colSums(vals[!is.na(vals[,1]),]))
}

sand_var <- function(gamma, eta, sigma_sq,
                     surv.data, long.data,
                     t_vals, rs, theta_vals, X_vals){
  
  #Define empty lists
  S <- list(NULL)
  S_i <- rep(0, nrow(t_vals))
  E0 <- list(NULL)
  Z <- list(NULL)
  Z_i <- rep(0, nrow(t_vals))
  for(j in 1:nrow(t_vals)){
    
    #keep a marker for the patient who's event time we're interested in
    #where they lie in the risk set
    patient_marker <- which(rs[[j]]== t_vals[j,"patient"])
    
    
    
    #sufficient statistic
    S[[j]] <- X_vals[[j]]
    S_i[j] <- S[[j]][patient_marker]+gamma*sigma_sq*theta_vals[[j]][patient_marker]
    S[[j]][patient_marker] <- S_i[j]
    
    
    #covariates
    Z[[j]] <- surv.data[surv.data$patient%in%rs[[j]],"treatment"]
    Z_i[j] <- Z[[j]][patient_marker]
    
    
    #condtitional density
    E0[[j]] <- exp(gamma*S[[j]]-gamma^2*sigma_sq*theta_vals[[j]]/2 +eta*Z[[j]])
    
    
  }
  
  #Vector E1 at each event time
  E1_1 <- mapply(function(x,y) x*y, S, E0)
  E1_2 <- mapply(function(x,y) x*y, Z, E0)
  
  #contribution at event time T_i
  B_vals <- matrix(0, nrow = nrow(t_vals), ncol = 2)
  for(i in 1:nrow(t_vals)){
    B_vals[i,1] <- S_i[i]-sum(E1_1[[i]])/sum(E0[[i]])
    B_vals[i,2] <- Z_i[i]-sum(E1_2[[i]])/sum(E0[[i]])
  }
  B <- cov(B_vals)
  #B <- cov(rbind(B_vals, matrix(0, n-nrow(t_vals), 2)))
  
  
  A <- jacobian(function(x){
    cond_score(x[1], x[2], sigma_sq, surv.data, long.data,
               t_vals, rs, theta_vals, X_vals)
  }, x=c(gamma, eta))/nrow(surv.data)
  
  return(list(A,B))
  
}

cs_grad <- function(gamma, eta, sigma_sq,
                    surv.data, long.data,
                    t_vals, rs, theta_vals, X_vals){
  
  
  A <- sand_var(gamma, eta, sigma_sq,
                surv.data, long.data,
                t_vals, rs, theta_vals, X_vals)[[1]]
  
  U <- cond_score(gamma, eta, sigma_sq,
                  surv.data, long.data,
                  t_vals, rs, theta_vals, X_vals)
  
  d_gamma <- 2*U[1]*A[1,1] + 2*U[2]*A[2,1]
  d_eta <- 2*U[1]*A[1,2] + 2*U[2]*A[2,2]
  
  return(c(d_gamma, d_eta))
}


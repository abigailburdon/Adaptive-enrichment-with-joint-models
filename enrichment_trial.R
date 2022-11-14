#Function which performs an enrichment trial and outputs treatment effect
#estimates and information levels.
#Input takes a generic function "treat_effect" which can be changed depending on
#analysis method
run_trial <- function(treat_effect, t_star, knot_point, n, gamma, eta1, eta2, delta1, delta2, sigma_sq,
                      h0_val, b2, b_sigma, lambda, events1, events2, C_val, I_max){
  
  #Generate data for first analysis
  all.data <- generate_ps(n, gamma_v = gamma, eta_f = c(eta1, eta2),
                          sigma_sq = sigma_sq, h0_val=h0_val,
                          recruit_end = n/(2*52), b2 = b2, b_sigma=b_sigma,
                          lambda = lambda, knot_point = knot_point)   
  surv.data <- all.data[[1]]
  long.data <- all.data[[2]]
  
  #Events based analysis - stop when n1 events reached
  t_vals <- surv.data[surv.data$status==1&surv.data$group==1,"surv_times"] + 
    surv.data[surv.data$status==1&surv.data$group==1,"arrival_times"] 
  analysis_times <- rep(0,2)
  analysis_times[1] <- t_vals[order(t_vals)[events1]]
  
  #remove patients who haven't arrived yet
  ps <- surv.data[surv.data$arrival_times <= analysis_times[1],"patient"]
  surv.data <- surv.data[ps,]
  long.data <- long.data[long.data$patient%in%ps,]
  #Relabel patients
  id.data <- data.frame(old.id = surv.data$patient, new.id = 1:nrow(surv.data))
  surv.data$patient <- id.data$new.id
  long.data$patient <- sapply(long.data$patient, function(x){
    id.data[id.data$old.id == x, "new.id"]
  })
  n1 <- nrow(surv.data)
  n2 <- n - n1
  
  #Data at first analysis
  all.data_1 <- data_cs(surv.data, long.data, analysis_times[1])
  surv.data_1 <- all.data_1[[1]]
  long.data_1 <- all.data_1[[2]]
  
  #Z_statistic for each subgroup
  TE1 <- try(treat_effect(subset(surv.data_1, group==1),
                          subset(long.data_1, group==1), eta1,
                          t_star, knot_point), silent=T)
  if(class(TE1)=="try-error"){
    TE1 <- rep(0,2)
  }
  theta_hat1 <- TE1[1]
  I1 <- TE1[2]
  Z1 <- -theta_hat1*sqrt(I1)
  TE2 <- try(treat_effect(subset(surv.data_1, group==2),
                          subset(long.data_1, group==2), eta2,
                          t_star, knot_point), silent=T)
  if(class(TE2)=="try-error"){
    TE2 <- rep(0,2)
  }
  theta_hat2 <- TE2[1]
  I2 <- TE2[2]
  Z2 <- -theta_hat2*sqrt(I2)
  #Z-statistic for full population
  theta_hatf <- lambda*theta_hat1+(1-lambda)*theta_hat2
  If <- 1/(lambda^2/I1 + (1-lambda)^2/I2)
  Zf <- -theta_hatf*sqrt(If)
  
  
  ###Threshold selection rule
  enrich <- c(Z1, Z2) > C_val
  no_inf <- (Z1==0)|(Z2==0)
  if(no_inf){
    TEf <- try(treat_effect(surv.data_1, long.data_1,
                            lambda*eta1+(1-lambda)*eta2,
                            t_star, knot_point), silent=T)
    if(class(TEf)=="try-error"){
      TEf <- rep(0,2)
    }
    theta_hatf <- TEf[1]
    If <- TEf[2]
    Zf <- -theta_hatf*sqrt(If)
    enrich <- c(T,T)
  }
  
  #Subgroup 1 selected
  if(enrich[1]&(!enrich[2])){
    
    new.data <- generate_ps(n2, gamma_v = gamma, eta_f = c(eta1, eta1),
                            sigma_sq = sigma_sq, h0_val=h0_val,
                            recruit_end = n2/(2*52), b2 = c(b2[1], b2[1]),
                            b_sigma=b_sigma, lambda = 1, knot_point = knot_point)
    new.data[[1]]$patient <- new.data[[1]]$patient+n1
    new.data[[1]]$arrival_times <- new.data[[1]]$arrival_times + analysis_times[1]
    new.data[[2]]$patient <- new.data[[2]]$patient+n1
    surv.data <- rbind(surv.data,new.data[[1]])
    long.data <- rbind(long.data,new.data[[2]])
    
    #Events based analysis - stop when n2 events in subgroup S1 reached
    t_vals <- surv.data[surv.data$status==1&surv.data$group==1,"surv_times"] + 
      surv.data[surv.data$status==1&surv.data$group==1,"arrival_times"] 
    analysis_times[2] <- t_vals[order(t_vals)[events2]]
    
    
    ### Data at analysis k=2 (can be extended for general K)
    all.data_k <- data_cs(surv.data, long.data, analysis_times[2])
    surv.data_k <- all.data_k[[1]]
    long.data_k <- all.data_k[[2]]
    n2 <- nrow(surv.data_k)
    
    #Treatment effect and information
    TE <- try(treat_effect(subset(surv.data_k, group ==1),
                           subset(long.data_k, group ==1), eta1,
                           t_star, knot_point), silent=T)
    if(class(TE)=="try-error"){
      TE <- rep(0,2)
    }
    I1 <- c(I1, TE[2])
    theta_hat1 <- c(theta_hat1, TE[1])
    
    #Estimates of information in non-selected subgroup
    obs_events <- sum(surv.data_1$group==2&surv.data_1$status==1)
    I2 <- c(I2, events2*I2/obs_events)
    obs_events <- sum(surv.data_1$status==1)
    If <- c(If, events2*If/obs_events)
    
    #Treatment effect and information from groups not selected
    TE2_follow <- try(treat_effect(subanalysis_times[1], set(surv.data_k, group ==2),
                                   subset(long.data_k, group ==2), eta2,
                                   t_star, knot_point), silent=T)
    if(class(TE2_follow)=="try-error"){
      TE2_follow <- rep(0,2)
    }
    TEf_follow <- c(lambda*theta_hat1[2]+(1-lambda)*TE2_follow[1],
                    1/(lambda^2/I1[2] + (1-lambda)^2/TE2_follow[2]))
      
    #Theta has is the one from subgroup S1
    theta_hat <- theta_hat1
    I <- I1
    
    #Treatment effect and information from groups not selected but followed-up
    theta_hat1_follow <- theta_hat1
    theta_hat2_follow <- c(theta_hat2, TE2_follow[1])
    theta_hatf_follow <- c(theta_hatf, TEf_follow[1])
    I1_follow <- I1
    I2_follow <- c(I2[1], TE2_follow[2])
    If_follow <- c(If[1], TEf_follow[2])
    
  }
  #Subgroup 2 selected
  if((!enrich[1])&enrich[2]){
    
    new.data <- generate_ps(n2, gamma_v = gamma, eta_f = c(eta2, eta2),
                            sigma_sq = sigma_sq, h0_val=h0_val,
                            recruit_end = n2/(2*52), b2 = c(b2[2], b2[2]),
                            b_sigma=b_sigma, lambda = 0, knot_point = knot_point)
    new.data[[1]]$patient <- new.data[[1]]$patient+n1
    new.data[[1]]$arrival_times <- new.data[[1]]$arrival_times + analysis_times[1]
    new.data[[2]]$patient <- new.data[[2]]$patient+n1
    surv.data <- rbind(surv.data,new.data[[1]])
    long.data <- rbind(long.data,new.data[[2]])
    
    #Events based analysis - stop when n2 events in subgroup S1 reached
    t_vals <- surv.data[surv.data$status==1&surv.data$group==2,"surv_times"] + 
      surv.data[surv.data$status==1&surv.data$group==2,"arrival_times"] 
    analysis_times[2] <- t_vals[order(t_vals)[events2]]
    
    ### Data at analysis k=2 (can be extended for general K)
    all.data_k <- data_cs(surv.data, long.data, analysis_times[2])
    surv.data_k <- all.data_k[[1]]
    long.data_k <- all.data_k[[2]]
    n2 <- nrow(surv.data_k)
    
    #Treatment effect and information
    TE <- try(treat_effect(subset(surv.data_k, group ==2),
                           subset(long.data_k, group ==2), eta2,
                           t_star, knot_point), silent=T)
    if(class(TE)=="try-error"){
      TE <- rep(0,2)
    }
    I2 <- c(I2, TE[2])
    theta_hat2 <- c(theta_hat2, TE[1])
    
    #Estimates of information in non-selected subgroup
    obs_events <- sum(surv.data_1$group==1&surv.data_1$status==1)
    I1 <- c(I1, events2*I1/obs_events)
    obs_events <- sum(surv.data_1$status==1)
    If <- c(If, events2*If/obs_events)
    
    #Treatment effect and information from groups not selected
    TE1_follow <- try(treat_effect(subset(surv.data_k, group ==1),
                                   subset(long.data_k, group ==1), eta1,
                                   t_star, knot_point), silent=T)
    if(class(TE1_follow)=="try-error"){
      TE1_follow <- rep(0,2)
    }
    TEf_follow <- c(lambda*TE1_follow[1]+(1-lambda)*theta_hat2[2],
                    1/(lambda^2/TE1_follow[2] + (1-lambda)^2/I2[2]))
      
    #Theta_hat is the one from subgroup S2
    theta_hat <- theta_hat2
    I <- I2
    
    #Treatment effect and information from groups not selected but followed-up
    theta_hat1_follow <- c(theta_hat1, TE1_follow[1])
    theta_hat2_follow <- theta_hat2
    theta_hatf_follow <- c(theta_hatf, TEf_follow[1])
    I1_follow <- c(I1[1], TE1_follow[2])
    I2_follow <- I2
    If_follow <- c(If[1], TEf_follow[2])
    
  }
  #Continue in full population 
  if(enrich[1]&enrich[2]){
    
    new.data <- generate_ps(n2, gamma_v = gamma, eta_f = c(eta1, eta2),
                            sigma_sq = sigma_sq, h0_val=h0_val,
                            recruit_end = n2/(2*52), b2 = c(b2[1], b2[2]),
                            b_sigma=b_sigma, lambda = lambda, knot_point = knot_point)
    new.data[[1]]$patient <- new.data[[1]]$patient+n1
    new.data[[1]]$arrival_times <- new.data[[1]]$arrival_times + analysis_times[1]
    new.data[[2]]$patient <- new.data[[2]]$patient+n1
    surv.data <- rbind(surv.data,new.data[[1]])
    long.data <- rbind(long.data,new.data[[2]])
    
    #Events based analysis - stop when n2 events in subgroup S1 reached
    t_vals <- surv.data[surv.data$status==1,"surv_times"] + 
      surv.data[surv.data$status==1,"arrival_times"] 
    analysis_times[2] <- t_vals[order(t_vals)[events2]]
    
    ### Data at analysis k=2 (can be extended for general K)
    all.data_k <- data_cs(surv.data, long.data, analysis_times[2])
    surv.data_k <- all.data_k[[1]]
    long.data_k <- all.data_k[[2]]
    n2 <- nrow(surv.data_k)
    
    #Treatment effect and information
    TE1 <- try(treat_effect(subset(surv.data_k, group == 1),
                            subset(long.data_k, group ==1), eta1,
                            t_star, knot_point), silent=T)
    if(class(TE1)=="try-error"){
      TE1 <- rep(0,2)
    }
    TE2 <- try(treat_effect(subset(surv.data_k, group == 2),
                            subset(long.data_k, group ==2), eta2,
                            t_star, knot_point), silent=T)
    if(class(TE2)=="try-error"){
      TE2 <- rep(0,2)
    }
    If <- c(If, 1/(lambda^2/TE1[2] + (1-lambda)^2/TE2[2]))
    theta_hatf <- c(theta_hatf, lambda*TE1[1]+(1-lambda)*TE2[1])
    
    #Estimates of information in non-selected subgroup
    obs_events <- sum(surv.data_1$group==1&surv.data_1$status==1)
    I1 <- c(I1, events2*I1/obs_events)
    obs_events <- sum(surv.data_1$group==2&surv.data_1$status==1)
    I2 <- c(I2, events2*I2/obs_events)
    
    #Treatment effect and information from groups not selected
    TE1_follow <- TE1
    TE2_follow <- TE2
      
    #Theta_hat is the one from the full population, F
    theta_hat <- theta_hatf
    I <- If
    
    #Treatment effect and information from groups not selected but followed-up
    theta_hat1_follow <- c(theta_hat1, TE1_follow[1])
    theta_hat2_follow <- c(theta_hat2, TE2_follow[1])
    theta_hatf_follow <- theta_hatf
    I1_follow <- c(I1[1], TE1_follow[2])
    I2_follow <- c(I2[1], TE2_follow[2])
    If_follow <- If
    
  }
  #No subgroup selected
  if((!enrich[1])&(!enrich[2])){
    #Result of trial
    theta_hat <- rep(0, 2)
    I <- rep(If, 2)
    result <- 0
    which_k <- 1
    n2 <- 0
    
    theta_hat1_follow <- c(theta_hat1, 0)
    theta_hat2_follow <- c(theta_hat2, 0)
    theta_hatf_follow <- c(theta_hatf, 0)
    I1_follow <- rep(I1[1], 2)
    I2_follow <- rep(I2[1], 2)
    If_follow <- rep(If[1], 2)
  }
  
  #Boundary points for any subgroup selected
  if(enrich[1]|enrich[2]){
    
    #Analysis of trial
    Z <- -theta_hat*sqrt(I)
    bounds <- error_spend(2,alpha,beta,C_val,delta1,delta2,lambda,I1, I2, If,I_max)
    cross_lower <- (Z < bounds[[1]])
    cross_lower_which <- ifelse(sum(cross_lower)==0, 3, min(which(cross_lower)))
    cross_upper <- (Z > bounds[[2]])
    cross_upper_which <- ifelse(sum(cross_upper)==0, 3, min(which(cross_upper)))
    result <- cross_upper_which < cross_lower_which
    which_k <- min(cross_lower_which, cross_upper_which)
  }
  
  
  ### Follow-up time
  fut_tot <- apply(surv.data, 1, function(x){
    min(analysis_times[which_k],x[3]+x[5])
  })
  fut <- fut_tot-surv.data$arrival_times
  mean_fut_arrived <- mean(fut[fut>0])
  
  
  ### Number of hospital visits
  mean_visits <- mean(sapply(surv.data[surv.data$arrival_times<analysis_times[which_k], "patient"], function(i){
    sub.i <- long.data[long.data$patient==i,]
    sum(sub.i$obs_time+surv.data[surv.data$patient==i, "arrival_times"] <= analysis_times[which_k])
  }))
  
  
  return(c(theta_hat, I, result, which_k, mean_fut_arrived, mean_visits, enrich,
           theta_hat1_follow,theta_hat2_follow,theta_hatf_follow,
           I1_follow,I2_follow,If_follow,n1,n2,analysis_times))
}
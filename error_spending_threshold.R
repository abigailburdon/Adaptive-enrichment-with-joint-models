#Choice of form for the error spending functions - taking input t which is the
#obtained fraction of the maximum information.
error_spend_f = function(t, a){
  return(min(a*t^2,a))}
error_spend_g = function(t, b){
  return(min(b*t^2,b))}

### Function which takes information levels and threshold criteria and returns
#the boundaries of the enrichment trial
error_spend <- function(
    K,           #number of groups
    alpha,       #required significance,
    beta,        #power requirement
    C_val,       #threshold selection value
    delta01,     #value of theta01 under alternative
    delta02,     #value of theta02 under alternative
    lambda,      #breakdown of Z3
    I1,I2,If,    #observed information levels
    I_max){      #maximum information
  
  
  #Pass anaylses that go backwards (or very close) in information
  for(k in 1:(K-1)){
    if(I1[k+1] < I1[k]){I1[k+1] <- I1[k]}
    if(I2[k+1] < I2[k]){I2[k+1] <- I2[k]}
    if(If[k+1] < If[k]){If[k+1] <- If[k]}}
  
  
  
  #############################
  ### Numerical integration ###
  #############################
  #gridpoints
  r_ind <- c(sapply(2:K, function(k) (If[k]-If[k-1])/If[k]),abs((I_max-If)/I_max))
  r_ind[r_ind=="NaN"] <- Inf
  if(min(r_ind)>0.01){r <- 16}
  if(min(r_ind)<=0.01){r <- 64}
  if(min(r_ind)<=0.001){r <- 128}
  
  x_unbounded <- grid(r,0)
  
  #Create empty vectors to store information
  h1 <- h2 <- hf <- list()  #density of continuation region
  z1 <- z2 <- zf <- list()  #list of grid points
  w1 <- w2 <- wf <- list()  #list of weights
  a_root <- rep(0, K) #lower boundaries
  b_root <- rep(0, K) #upper boundaries
  
  
  #Desired amount of error to be spent  
  #Type 1 error
  pi_1k <- rep(0,K)
  pi_1k[1] <- error_spend_f(If[1]/I_max, alpha)
  pi_1k[2:K] <- sapply(2:K, function(k)
  {error_spend_f(If[k]/I_max, alpha)-error_spend_f(If[k-1]/I_max, alpha)})
  #Type 2 error
  pi_2k <- rep(0,K)
  pi_2k[1] <- error_spend_g(If[1]/I_max, beta)
  pi_2k[2:K] <- sapply(2:K, function(k)
  {error_spend_g(If[k]/I_max, beta)-error_spend_g(If[k-1]/I_max, beta)})
  
  
  #Avoid underrunning by spending all error at final anaylsis
  if(If[K] < I_max){
    pi_1k[K] <- alpha-error_spend_f(If[K-1]/I_max, alpha)} 
  
  #########################
  ### Density functions ###
  #########################
  norm1 <- function(delta01, delta02, C_val){
    (1-pnorm(C_val, mean=delta01*sqrt(I1[1]), sd=1))*
      pnorm(C_val, mean=delta02*sqrt(I2[1]), sd=1)
  }
  f_Z1 <- function(z1, delta01, delta02, C_val){
    #Returns joint density of selecting subgroup 1 and Z1
    h <- dtruncnorm(z1, a=C_val, b=Inf, mean = delta01*sqrt(I1[1]), sd = 1)
    return(h*norm1(delta01, delta02, C_val))
  }
  
  norm2 <- function(delta01, delta02, C_val){
    pnorm(C_val, mean=delta01*sqrt(I1[1]), sd=1)*
      (1-pnorm(C_val, mean=delta02*sqrt(I2[1]), sd=1))
  }
  f_Z2 <- function(z2, delta01, delta02, C_val){
    #Returns joint density of selecting subgroup 2 and Z2
    h <- dtruncnorm(z2, a=C_val, b=Inf, mean = delta02*sqrt(I2[1]), sd = 1)
    return(h*norm2(delta01, delta02, C_val))
  }
  
  norm3 <- function(delta01, delta02, C_val){
    (1-pnorm(C_val, mean=delta01*sqrt(I1[1]), sd=1))*
      (1-pnorm(C_val, mean= delta02*sqrt(I2[1]), sd=1))
  }
  f_Z3 <- function(z3, delta01, delta02, C_val){
    #Returns joint density of selecting full population and Z3
    const1 <- 1/sqrt(lambda)
    const2 <- 1/sqrt(1-lambda)
    h <- integrate(function(x){
      p_val1 <- dtruncnorm(x/sqrt(lambda), a=C_val, b=Inf, mean = delta01*sqrt(I1[1]), sd=1)
      p_val2 <- dtruncnorm((z3-x)/sqrt(1-lambda), a=C_val, b=Inf,
                           mean = delta02*sqrt(I2[1]),sd=1)
      return(p_val1*p_val2)
    }, lower = C_val*sqrt(lambda), upper = z3-C_val*sqrt(1-lambda))$value
    return(const1*const2*h*norm3(delta01, delta02, C_val))
  }
  f_Z3_Vec <- Vectorize(f_Z3, vectorize.args = "z3")
  
  f_Z2gZ1 <- function(z1,z2,I1,I2,theta){
    #Returns the probability of moving from z1 in analysis k-1 to z2 in analysis k.
    dI <- I2-I1 #Information difference
    f_val <- sqrt(I2)*dnorm((z2*sqrt(I2)-z1*sqrt(I1)-theta*dI)/sqrt(dI))/sqrt(dI)
    return(f_val)
  }
  
  spent_error <- function(b,I2,I1,h1,z1,theta){
    #Returns the probability of being below b
    #Gridpoints and weights
    x <- grid_bounded(x_unbounded, c(-Inf, b))
    z2 <- grid_z(x)
    w2 <- weights(z2)
    #Density 
    h_current <- sapply(1:length(z2), function(i){
      fs <- f_Z2gZ1(z1,z2[i],I1,I2,theta)
      return(sum(h1*fs*w2[i]))})
    return(sum(h_current))}
  
  ##################
  ### Analysis 1 ###
  ##################
  
  #P(R) at k=1
  b1 <- uniroot(function(b1){
    integrate(function(z){f_Z1(z,0,0,C_val)}, b1, Inf)$value+
      integrate(function(z){f_Z2(z,0,0,C_val)}, b1, Inf)$value+
      integrate(function(z){f_Z3_Vec(z,0,0,C_val)}, b1, Inf)$value-
      pi_1k[1]}, lower = 0, upper = 5, extendInt = "yes")$root
  b_root[1] <- b1
  
  
  #P(A) at k=1
  a1 <- uniroot(function(a1){
    val <- integrate(function(z){f_Z1(z,delta01,delta02,C_val)},C_val, a1)$value
    val/norm1(delta01, delta02, C_val)-pi_2k[1]
  }, lower = C_val, upper = 10, extendInt = "yes")$root
  
  
  #a1 <- uniroot(function(a1){
  #  val <- integrate(function(z){
  #    f_Z1(z,delta01,delta02,C_val)+f_Z3_Vec(z,delta01,delta02,C_val)},-Inf, a1)$value/
  #    (norm1(delta01, delta02, C_val)+norm3(delta01, delta02, C_val))
  #  val-pi_2k[1]}, lower = C_val, upper = 10, extendInt = "yes")$root
  
  #a1 <- uniroot(function(a1){
  #  val <- -10
  #  while(val==-10){
  #    val <- try(integrate(function(z){f_Z1(z,delta01,delta02,C_val)},
  #                         -Inf, a1)$value,silent=T)
  #    if(class(val)=="try-error"){
  #      val <- -10
  #      a1 <- a1+1e-6
  #    }
  #  }
  #  val-pi_2k[1]
  #}, lower = C_val, upper = 10, extendInt = "yes")$root
  a_root[1] <- a1
  
  #Rare case where boundaries cross at analysis 1
  k <- 2
  if(a1 > b1 | I1[k] <= I1[k-1] | I2[k] <= I2[k-1] | If[k] <= If[k-1]){
    K_final <- 1
  }else{
    K_final <- K
    
    ### Analysis 1
    x <- grid_bounded(x_unbounded, c(a_root[1],b_root[1]))
    z1[[1]] <- z2[[1]] <- zf[[1]] <- grid_z(x)
    w1[[1]] <- weights(z1[[1]])
    w2[[1]] <- weights(z2[[1]])
    wf[[1]] <- weights(zf[[1]])
    
    #Densities
    h1[[1]] <- w1[[1]]*f_Z1(z1[[1]], 0, 0, C_val)
    h2[[1]] <- w2[[1]]*f_Z2(z2[[1]], 0, 0, C_val)
    hf[[1]] <- wf[[1]]*f_Z3_Vec(zf[[1]], 0, 0, C_val)
    
    ############################
    ### Analyses k = 2,...,K ###
    ############################
    
    k <- 2
    #P(R) at k=2
    b_root[k] <- uniroot(function(b){
      sum(h1[[k-1]])+sum(h2[[k-1]])+sum(hf[[k-1]])-
        spent_error(b,I1[k],I1[k-1],h1[[k-1]],z1[[k-1]],0)-
        spent_error(b,I2[k],I2[k-1],h2[[k-1]],z2[[k-1]],0)-
        spent_error(b,If[k],If[k-1],hf[[k-1]],zf[[k-1]],0)-
        pi_1k[k]
    }, lower = 0, upper = 5, extendInt = "yes")$root
    
  }
  
  a_root[K_final] <- b_root[K_final]
  
  
  return(list(a_root, b_root, pi_1k, pi_2k))
}

#Function which calculates required number of total events 
n2_calc <- function(
    K,
    alpha,
    beta,
    delta01,
    delta02,
    C_val,
    lambda,
    sigma1,
    sigma2,
    I1_val
){
  
  #Find the inflation factor needed to get power=1-beta
  n2_val <- uniroot(function(n2){
    
    I1 <- c(I1_val, n2/sigma1)
    I2 <- c(I1_val*(1-lambda)/lambda, n2/sigma2)
    If <- c(I1_val/lambda, n2/(lambda*sigma1+(1-lambda)*sigma2))
    I_max <- If[2]
    
    bounds <- error_spend(K,alpha,beta,C_val,delta01,delta02,lambda,I1,I2,If,I_max)
    a_root <- bounds[[1]]
    b_root <- bounds[[2]]
    
    #gridpoints
    r_ind <- c(sapply(2:K, function(k) (If[k]-If[k-1])/If[k]), abs((I_max-If)/I_max))
    r_ind[r_ind=="NaN"] <- Inf
    if(min(r_ind)>0.01){r <- 16}
    if(min(r_ind)<=0.01){r <- 64}
    if(min(r_ind)<=0.001){r <- 128}
    x_unbounded <- grid(r,0)
    
    #Create empty vectors to store information
    h1 <- list()       #density of continuation region
    z <- list()        #list of grid points
    w <- list()        #list of weights
    x <- grid_bounded(x_unbounded, c(a_root[1],b_root[1]))
    z[[1]] <- grid_z(x)
    w[[1]] <- weights(z[[1]])
    
    
    #########################
    ### Density functions ###
    #########################
    norm1 <- function(delta01, delta02, C_val){
      (1-pnorm(C_val, mean=delta01*sqrt(I1[1]), sd=1))*
        pnorm(C_val, mean=delta02*sqrt(I2[1]), sd=1)
    }
    f_Z1 <- function(z1, delta01, delta02, C_val){
      #Returns joint density of selecting subgroup 1 and Z1
      h <- dtruncnorm(z1, a=C_val, b=Inf, mean = delta01*sqrt(I1[1]), sd = 1)
      return(h*norm1(delta01, delta02, C_val))
    }
    f_Z2gZ1 <- function(z1,z2,I1,I2,theta){
      #Returns the probability of moving from z1 in analysis k-1 to z2 in analysis k.
      dI <- I2-I1 #Information difference
      f_val <- sqrt(I2)*dnorm((z2*sqrt(I2)-z1*sqrt(I1)-theta*dI)/sqrt(dI))/sqrt(dI)
      return(f_val)
    }
    spent_error <- function(b,I2,I1,h1,z1,theta){
      #Returns the probability of being below b
      #Gridpoints and weights
      x <- grid_bounded(x_unbounded, c(-Inf, b))
      z2 <- grid_z(x)
      w2 <- weights(z2)
      #Density 
      h_current <- sapply(1:length(z2), function(i){
        fs <- f_Z2gZ1(z1,z2[i],I1,I2,theta)
        return(sum(h1*fs*w2[i]))})
      return(sum(h_current))}
    
    #P(R|S=1)
    h1[[1]] <- w[[1]]*f_Z1(z[[1]], delta01, delta02, C_val)
    
    k <- 2
    val <- (integrate(function(z){f_Z1(z,delta01,delta02,C_val)}, b_root[1], Inf)$value+
              sum(h1[[1]])-
              spent_error(b_root[2],I1[k],I1[k-1],h1[[k-1]],z[[k-1]],delta01))/
      norm1(delta01, delta02, C_val)
    
    return(val-1+beta)
  },
  lower=ceiling(I1_val*sigma1), upper=200, extendInt = "yes")$root
  return(n2_val)
}
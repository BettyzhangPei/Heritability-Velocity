Final.result.AI.ReML <- function(initial.par, parts, part.result.AI.ReML, Data)
{ 
  
  # parts: number of parts need to be divided for whole datasets
  # final.result.AI.ReML: function to estimate all variance components, two heritability metrics and their variances
  # Data: the full dataset used for analysis 
   
   G1 = Data$G
   N0 = dim(G1)[1]
   
   
   result_part = matrix(0, nrow= 16, ncol = parts)
   
   for (a in 1: parts)
   {
     
   # Randomly shuffle numbers: 
   shuffled_numbers<- sample(1:N0)

  
   # compute the size of each part (integer division)
   part_size<- N0/parts


   # Initialize a matrix to store the parts
   divided_parts <- matrix(0, nrow= part_size , ncol= parts)

   # Loop to split the numbers
   for (i in 1:parts) 
   {
    # Each part gets 'part_size' elements, with increasing order: 
    divided_parts[,i] <- sort(shuffled_numbers[((i - 1) * part_size + 1):(i * part_size)])
   }

    N<- part_size
    # for ath part:
    s<- rep(0,part_size)
    s<-  divided_parts[,a]

    # Call function to generate required information for each part to do analysis
    Data1 = Generate_information(s, Data)
    
    # Call function to conduct AI-REML algorithm for each part separately
    result_part[,a] = part.result.AI.ReML(initial.par, data=Data1, AI_ReML, f_V, AI_DL)
    
   }
  
   
    # Estimates with 6 digits  
    result_part = round(result_part, 6)
   
    # Truncation-adjusted meta-analysis to combine estimates from all parts and then derive final estimates: 
    # Conduct left-truncated meta-analysis for slope variance components except variance of residuals
    # Conduct doubly truncated meta-analysis for heritability metrics on dynamic aspects of longitudinal phenotypes
    
    
    est.result.g.velocity = left.truncation.meta.analysis(X=result_part[2,], sd= sqrt(result_part[10,]), tol=1e-6)
    est.result.b1 = left.truncation.meta.analysis(X=result_part[4,], sd= sqrt(result_part[14,]), tol=1e-6)
    est.result.h.velocity= doubly.truncation.meta.analysis(X=result_part[8,], sd= sqrt(result_part[12,]), tol=1e-4)
    
    # Combine results using truncated meta-analysis
    return(list(var.g.base = mean(result_part[1,]) ,
               var.g.velocity =  est.result.g.velocity$est.mu,
               var.b0 =  mean(result_part[3,]),
               var.b1 = est.result.b1$est.mu,
               var.e = mean(result_part[5,]),
               convergence = mean(result_part[6,]) , 
               h.base =  mean(result_part[7,]),
               h.velocity =  est.result.h.velocity$est.mu,
               sd.var.g.base =  mean(result_part[9,]),
               sd.var.g.velocity = est.result.g.velocity$est.sd,
               sd.h.base = mean(result_part[11,]), 
               sd.h.velocity = est.result.h.velocity$est.sd,
               cov.g.velocity.b1=  mean(result_part[13,]),
               sd.var.b0 = mean(result_part[14,]), 
               sd.var.b1 = est.result.b1$est.sd, 
               sd.var.e = mean(result_part[16,])
               ))
   
}





####################################################################################################
##### 1. Meta-analysis for left-truncated values 
left.truncation.meta.analysis <- function(X, sd, tol=1e-6)
{
  
  ###### objective function:  negative log likelihood to be minimized in further optim function ####
  objective  = function(par, X, sd)
  {
    # par : mean of the normal distribution to be estimate
    # sd: standard deviation of the normal distribution
    # X: random variables from a normal distribution
    # X_trun : truncated variables at zero
    
    
    # Total number of elements in the vector X with random variables: 
    n = length(X)
    
    # negative variables from X
    X_trun =  X[which(X<=tol)]
    
    # number of truncated values: nn
    nn = length(X_trun)
    
    # probability of X more than 0: 
    s0 = which(X<=tol)
    sd_trun = sd[s0]
    p0 = pnorm(par/sd_trun)
    
    
    # Remove the subset (with truncated values at zero) from the whole set: 
    X_rest = setdiff(X, X_trun)
    
    
    s = which(X>tol)
    sd_rest = sd[s]
    
    
    # log likelihood function:
    l =  sum(log(1- p0))  + length(s)* (- 0.5 *log(2*pi))  - 0.5*sum(log(sd_rest^2)) -  sum( ((X_rest - par)^2)/(2*sd_rest^2) )
    
    # return negative log likelihood to be minimized for MLE
    return(-l)
    
  }
  ###########################################################################################################################
  
  
  
  ###### The function to obtain the estimated mu and information related to convergence ##########
  f<- function(objective, X, sd)
  {
    
    # a: an order number for ath repetition for further simulation studies
    # objective: a function to characterize the negative log likelihood 
    # n : sample size 
    # mu: true value of mean 
    # sd:  values of standard deviation for each random variable 
    # nn:  number of truncated values 
    
    # sample size of estimates: 
    n= length(X)  
    #########################################################
    
    m = max(X)
    
    ########################### Optimization ######################################
    # Choose a random continuous value between   mu-1.96*sd  and mu + 1.96*sd
    initial.mu = runif(1, min =  0, max =m )
    # Define the lower and upper bounds for the "Brent" method
    lower_bound <-  0
    upper_bound <-  m
    
    # optim function to obtain the MLE for mu: 
    result <- optim(par= initial.mu, fn = objective, X=X, sd=sd,
                    method = "Brent", lower = lower_bound, upper = upper_bound)
    
    est.mu = result$par 
    conver = result$convergence
    
    df<- rep(0,2)
    df[1]<- est.mu
    df[2]<- conver 
    
    # Return estimated mu and the related value to indicate convergence
    return(df)  
    
  }
  
  
  
  # Second derivative for log likelihood: 
  l_2 <- function(mu, X, sd)
  {
    # probability of X more than 0: 
    s0 = which(X<=tol)
    sd_trun = sd[s0]
    # CDF:
    p0 = pnorm(mu/sd_trun)
    # density : 
    p1 = dnorm(mu/sd_trun)
    
    X_trun= X[which(X<=tol)]
    # Remove the subset (with truncated values at zero) from the whole set: 
    X_rest = setdiff(X, X_trun)
    
    s = which(X>tol)
    sd_rest = sd[s]
    
    l2 = sum(  (mu*p1*(1-p0) - p1^2)/(sd_trun*(1-p0))^2 ) - sum(1/sd_rest^2)
    
    var = -1/l2
    return(var)
  } 
  
  
  result = f(objective, X, sd)
  
  est.sd = sqrt(l_2(result[1], X, sd))
  
  # Estimated mu using the left-truncated meta-analysis approach:
  return(list(est.mu = result[1], convergence= result[2],  est.sd = est.sd))
  
}

####################################################################################################



####################################################################################################
##### 1. Meta-analysis for doubly truncated values 
doubly.truncation.meta.analysis <- function(X, sd, tol=1e-6)
{
   
  # Objective function: negative log likelihood to be minimized in further optim function
  objective1 <- function(par, X, sd) {
  # par : mean of the normal distribution to be estimated
  # sd: standard deviation of the normal distribution
  # X: random variables from a normal distribution

  # Total number of elements in the vector X with random variables
  n <- length(X)



  # Probability of X less than 0
  s00 <- which(X <= tol)
  if (length(s00) > 0) {
    sd_trun0 <- sd[s00]
    p0 <- pnorm(-par / sd_trun0)
    p0[p0 == 0] <- .Machine$double.eps  # Avoid log(0)
    log_p0 <- sum(log(p0))
  } else {
    log_p0 <- 0  # No values less than 1e-4, so contribution is zero
  }

  # Probability of X more than 1
  s1 <- which(X > 0.99)
  if (length(s1) > 0) {
    sd_trun1 <- sd[s1]
    p1 <- 1 - pnorm((1 - par) / sd_trun1)
    p1[p1 == 0] <- .Machine$double.eps  # Avoid log(0)
    log_p1 <- sum(log(p1))
  } else {
    log_p1 <- 0  # No values more than 0.98, so contribution is zero
  }

  # Remove the subset (with truncated values from 0 to 1) from the whole set
  X_rest <- X[X > tol & X <= 0.99]

  # Log likelihood calculation for non-truncated values
  s0 <- 0
  if(length(X_rest) > 0 ){
  for (i in 1:length(X_rest)) {
    density_value <- dnorm((X_rest[i] - par) / sd[i], log = TRUE)
     # if (is.infinite(density_value)) {
     #   density_value <- -Inf  # Handle potential -Inf from log(0)
     # }
    s0 <- s0 + density_value
    }
  }

  s <- which(X > tol & X <= 0.99)
  log_sd <- if (length(s) > 0) sum(-log(sd[s])) else 0  # Avoids errors if s is empty
  
  
  ## Log likelihood function
  # l <- sum(log(p0)) + sum(log(p1)) + sum(-log(sd[s])) + s0
  l <- log_p0 + log_p1 + log_sd + s0
  
  # Return negative log likelihood to be minimized for MLE
  return(-l)
}



  # Function to obtain the estimated mu and information related to convergence
  f2 <- function(objective, X, sd) {
  # objective: a function to characterize the negative log likelihood
  # n: sample size
  # sd: true value of standard deviation in the same normal distribution

  # Choose a random continuous value between mu - 1.96 * sd and mu + 1.96 * sd
  initial.mu <- runif(1, min = 0, max = 1)

  # Define the lower and upper bounds for the "L-BFGS-B" method
  lower_bound <- 0
  upper_bound <- 1

  # Optim function to obtain the MLE for mu
  result <- optim(par = initial.mu, fn = objective, X = X, sd = sd,
                  method = "L-BFGS-B", lower = lower_bound, upper = upper_bound)

  est.mu <- result$par
  conver <- result$convergence

  df <- rep(0, 2)
  df[1] <- est.mu
  df[2] <- conver

  # Return estimated mu and the related value to indicate convergence
  return(df)
}





# second derivative for log likelihood related to doubly-truncated meta analysis for ratio2:
l_2_lambda <- function(par, sd, X)
{


  # Total number of elements in the vector X with random variables
  n <- length(X)

  # Probability of X less than 0
  s00 <- which(X <= tol)
  if(length(s00) > 0)
  {
    sd_trun0 <- sd[s00]
    z0 = -par / sd_trun0
    p0 = pnorm(z0)
    p00 = dnorm(z0)
    p0[p0 == 0] <- .Machine$double.eps  # Avoid log(0)
    term1 =   - sum( (z0* p00*p0  +  p00^2 ) / (p0*sd_trun0)^2  ) 
  }else{term1=0}

  # Probability of X more than 1
  s1 <- which(X > 0.99)
  if(length(s1) >0)
  {
    sd_trun1 <- sd[s1]
    p1 <- 1 - pnorm((1 - par) / sd[s1])
    p11 = dnorm((1 - par) / sd[s1] )
    p1[p1 == 0] <- .Machine$double.eps  # Avoid log(0)
    term2 = sum((p11^2 - ((1 - par) / sd_trun1) * p11 * p1) / (p1 * sd_trun1)^2, na.rm = TRUE)
  }else{term2=0}
  
  s <- which(X > tol & X <= 0.99)
  if (length(s) > 0)
  {
    term3 <- - sum(1 / (sd[s]^2), na.rm = TRUE)
  }else{term3 =0}
  
  
  
  ## second derivative of log likelihood for mu:
  # l2 <- - sum( ( (-par/sd_trun0)* p00*p0  + p00^2   ) / (p0*sd_trun0)^2  )
  # + sum(  (p11^2 -  ((1 - par) / sd[s1])* p11*p1  ) / (p1*sd_trun1)^2  )
  # - sum(1/(sd[s]^2))
  l2 = term1 + term2 + term3

  # Handle NA or infinite case for the second derivative: 
    if (is.nan(l2) || is.infinite(l2)) {
      warning("Warning:  the second derivative of likelihood function is NA or infinite ")
        return(NA)  # Return NA for invalid or unstable values
    }
  # Handle non-negative for the second derivative:
    if(l2 >= 0) 
    { 
       warning("Warning: the second derivative of likelihood function is non-negative, implying that is convex at the estimated parameter value.")
       return(NA)  # Return NA for invalid or unstable values
    }
  
  
    # Calculate the variance as the negative inverse of l2
    var <- -1 / l2
  
    return(var)
}
  
    result = f2(objective1, X, sd)
    
    est.sd = sqrt(l_2_lambda(result[1], sd,X))
    
    # Estimated mu using the doubly truncated meta-analysis approach:
    return(list(est.mu = result[1], convergence= result[2],  est.sd = est.sd))

}
#########################################################################################################
#########################################################################################################


















   
   
   
   
   
   
   
   




part.result.AI.ReML <- function(initial.par, data, AI_ReML, f_V, AI_DL)
{
  ### initial.par : the input of variance components for estimating variance components via AI-REML algorithm 
  # where unknown variance components are sigma^2_g, sigma^2_g*, sigma^2_b0, sigma^2_b1 and sigma^2_e
  
  ### data: generated data information with a list of matrices:  
  # G : genetic relationship matrix 
  # t : the sum(n) x 1 column vector of time effect
  # n : the N x 1 column vector represents total number of measurements for each subjects
  # y : the  sum(n) x 1 column vector of response 
  # H : composite matrix for AI matrix and DL matrix in the iteration process of AI-ReML algorithm
  # A : covariate matrix for fixed effects
  
  # f_V: function for covariance matrix V 
  # AI_DL : function for computing AI and DL matrices
  N = length(data$n)
  
  result <- AI_ReML(par = initial.par, l_REML,  maxit = 1000, maxtol = 1e-4, data = data, f_V, AI_DL)
  
  # Estimated variance components = (sigma^2_g, sigma^2_g*, sigma^2_b0, sigma^2_b1, sigma^2_e)
  est.par= as.matrix(result$par, ncol=1)
  # Estimated two heritability metrics and their standard errors
  
  
  
  est.var.genetics<- numeric(2)
  
  est.AI_DL <- AI_DL(par=est.par,data, f_V)
  est.AI <- est.AI_DL$AI
  
  
  # inverse.est.AI<- solve(-est.AI,tol=1e-30)
  inverse.est.AI<- solve(-est.AI)
  
  est.var.genetics[1]<- inverse.est.AI[1,1]
  est.var.genetics[2]<- inverse.est.AI[2,2]
  
  
  
  
  # estimated ratios: lambda_base, lamba_slope: 
  est.ratio <- c(est.par[1]/(est.par[1] + est.par[3]),  est.par[2]/ (est.par[2] + est.par[4]) )
  
  est.xi <- c(est.ratio[1], est.ratio[2], est.par[3]+est.par[1], est.par[4]+est.par[2], est.par[5])
  
  
  
  ######################################
  ###  Jacobian matrix K ###
  n0=length(est.par)
  K =  matrix(0,nrow=n0,ncol=n0)
  K[1,1] = est.xi[3]
  K[1,3] = est.xi[1]
  K[2,2] = est.xi[4]
  K[2,4] = est.xi[2]
  K[3,1] = -est.xi[3]
  K[3,3] = 1-est.xi[1]
  K[4,2] = -est.xi[4]
  K[4,4] = 1- est.xi[2]
  K[5,5] = 1
  ######################################
  
  # estimated AI matrix for est.xi using delta method with Jacobian matrix K:
  est.AI.xi<- t(K)%*%est.AI%*%K
  
  
  ######################################
  est.inverse.AI.xi<- solve(- est.AI.xi)
  
  est.var.ratio <- numeric(2)
  est.var.ratio[1]<- est.inverse.AI.xi[1,1]
  est.var.ratio[2]<- est.inverse.AI.xi[2,2]
  ######################################
  
  df<- rep(0,n0+11)
  df[1:n0]<- est.par
  df[n0+1]<- result$convergence
  
  df[n0+2]<- est.ratio[1]
  df[n0+3]<- est.ratio[2]
  
  df[n0+4] <- est.var.genetics[1]
  df[n0+5] <- est.var.genetics[2]
  
  df[n0+6]<- est.var.ratio[1]
  df[n0+7]<- est.var.ratio[2]
  df[n0+8]<- inverse.est.AI[2,4]
  
  df[n0+9] =  inverse.est.AI[3,3]
  df[n0+10] = inverse.est.AI[4,4]
  df[n0+11] = inverse.est.AI[5,5]
  
  df = round(df, 6)
  
  #df <- data.frame(value = df)
  
  # Set custom column names
  #rownames(df) <- c("var.g", "var.g*", "var.b0", "var.b1", "var.e", 
                   # "convergence", 
                   #  "base.heritability","velocity.heritability", 
                   # "var.var.g", "var.var.g*", 
                   # "var.h1", "var.h2", "cov.g*.b1", 
                   # "var.var.b0", "var.var.b1", "var.var.e")
  
  return(df)
  
}












#################### adapted  AI-ReML algorithm #############
AI_ReML<- function(par, l_REML, maxit, maxtol, data, f_V, AI_DL)
{ 
  # par : unknown variance components (sigma^2_g, sigma^2_g*, sigma^2_b0, sigma^2_b1, sigma^2_e)
  # l_REML : function for computing restricted log likelihood
  # maxit : maximum of iterations 
  # maxtol : tolerance for convergence 
  
  ### data: generated data information with a list of matrices:  
  # G : genetic relationship matrix 
  # t : the sum(n) x 1 column vector of time effect
  # n : the N x 1 column vector represents total number of measurements for each subjects
  # y : the  sum(n) x 1 column vector of response 
  # H : composite matrix for AI matrix and DL matrix in the iteration process of AI-ReML algorithm.
  # A : covariate matrix for fixed effects
  
  # f_V: function for covariance matrix V 
  # AI_DL : function for computing AI and DL matrices
  
  G = data$G
  t = data$t
  n = data$n
  y = data$y
  H = data$H
  A = data$A
  
  
  
  
  
  # Assuming 'par' is a vector and n0 is its length
  n0 = length(par)
  
  # initial setup
  old_par<- par
  
  
  # Function to update parameters
  update_parameters <- function(par) {
    AD <- AI_DL(par, data, f_V)
    new_par <- par - solve(AD$AI) %*% AD$DL
    # Ensure parameters are non-negative
    new_par <- pmax(new_par, 1e-6)  # ensure positivity 
    return(new_par)
  }
  
  
  # Initial log-likelihood
  old_ll <- l_REML(old_par, H, A, y, n, f_V) 
  
  
  # Initial update
  new_par <- update_parameters(old_par)
  new_ll <- l_REML(new_par, H, A, y, n, f_V)
  
  
  tol <- abs(new_ll - old_ll)
  it <- 1
  
  
  # Iterative update until convergence
  while (it < maxit && tol > maxtol) {
    old_par <- new_par
    old_ll <- new_ll
    
    new_par <- update_parameters(old_par)
    new_ll <- l_REML(new_par, H, A, y, n, f_V)
    
    tol <- abs(new_ll - old_ll)
    it <- it + 1
  }
  
  # Check convergence status
  convergence <- ifelse(tol < maxtol, 0, 1)
  
  
  return(list(par=new_par, convergence=convergence, tol=tol))
  
}

############################################################



# function for computing covariance matrix  V  
f_V <- function(par, H, n)
{
  
  
  H1 =  H[1: sum(n), ] 
  H2 =  H[(sum(n) + 1): (2*sum(n)), ] 
  H3 =  H[ (2*sum(n) + 1): (3*sum(n)), ] 
  H4 =  H[ (3*sum(n) + 1): (4*sum(n)), ]
  H5 =  H[ (4*sum(n) + 1): (5*sum(n)), ] 
  
  V = par[1]*H1 + par[2]*H2 + par[3]*H3 + par[4]*H4 + par[5]*H5
  
  return(V)
}




##### Restricted log likelihood ######
l_REML<-function(par, H, A, y, n, f_V)
{
  # covariance matrix V
  V =  f_V(par, H, n)
  
  # upper triangular matrix by using the cheolsky decomposition for V
  L = chol(V)
  
  
  # inverse of V
  V1 = chol2inv(L)
  
  # a1 = log(det(V))
  a1 =  2*sum(log(diag(L)))
  
  
  R0 = V1%*%A
  
  # upper triangular matrix by using the cheolsky decomposition
  L3 = chol(t(A)%*%(R0))
  
  # inverse of t(A)%*%V1%*%A
  V2 = chol2inv(L3)
  
  # Projection matrix
  R1 =  V1 - (R0%*%V2)%*%t(R0) 
  
  # log(det(V2)) 
  a2 = 2*sum(log(diag(L3)))
  
  
  # The restricted log likelihood:   
  l= -(0.5)*( (t(y) %*% R1) %*% y +  a1 + a2)
  
  
  return(as.numeric(l))
  
}






#######################################################################################
## Function for first derivatives and average information related to second derivatives 
AI_DL<- function(par, data, f_V)
{
  
  G= data$G
  t= data$t
  n= data$n
  y= data$y
  H= data$H
  A= data$A
  
  V = f_V(par, H, n)
  
  L = chol(V)
  
  V1 =  chol2inv(L)
  
  R0 = V1%*%A
  
  V2 = solve(t(A)%*%(R0))
  
  # Projection matrix
  R1 =  V1 - (R0%*%V2)%*%t(R0) 
  
  
  
  n0 = length(par)
  AI_matrix = matrix(0, n0, n0)
  DL_matrix = numeric(n0) 
  
  
  
  R2 = R1 %*%y
  
  
  
  for (i in 1:n0)
  {
    ind1 = i*sum(n)
    ind2 = (i-1)*sum(n)
    ind_set_i = (ind2+1) : ind1
    
    DL_matrix[i] =  as.numeric(-0.5 *( sum(R1*H[ind_set_i, ]) - (t(R2) %*% H[ind_set_i, ]) %*%R2 ))
    
    
    for (j in 1:n0)
    {
      ind3 = j*sum(n)
      ind4 = (j-1)*sum(n)
      ind_set_j = (ind4+1) : ind3 
      
      AI_matrix[i,j] =  as.numeric(-0.5 * ((t(R2) %*% H[ind_set_i, ]) %*%R1) %*% (H[ind_set_j, ]%*%R2))
      
    }
  }
  
  
  
  return(list(AI= AI_matrix, DL=DL_matrix))
}


################################################################################################################
#########################################################################################################
#########################################################################################################    
# Function: Generate information for further implementing AI-REML algorithm for each part separately. 
Generate_information = function(s, Data)
{
  
  G1  = Data$G    # The genetic relationship matrix generated by genome-wide variants for all subjects from full dataset
  G01 = Data$G0   # The genetic relationship matrix generated by causal variants for all subjects from full dataset
  y   = Data$y     # The vector of observed values for phenotypic traits
  t   = Data$t     # The vector of observed/calculated time variable
  n   = Data$n     # Vector of elements and each element represents the total number of measurements per subject
  N   = length(s)  # Number of subjects included in each part
  
  
  ### The genotypic information matrix for each part:
  G<- matrix(0,nrow=N, ncol=N)
  
  for (i in 1:N)
  {
    for (j in 1:N)
    {
      G[i,j]<-  as.numeric(G1[s[i], s[j]])
    }
  }
  
  
  
  # count number of non-NA measurements for each male  
  n1 <- n[s]
  
  
  ## Helper function to calculate the cumulative sum of n up to the (i-1)th element
  cumulative_sum<- function(i,n)
  {
    ifelse(i==1, 0, sum(n[1:(i-1)]) )
  }
  
  y1 = numeric(sum(n1))
  t1 = numeric(sum(n1))
  
  ####### Access longitudinal settings ############
  for ( i in 1:N)
  {
    
    index_start_i <- cumulative_sum(i,n1) 
    index_end_i <- cumulative_sum(i+1,n1)
    index_set_i = (index_start_i +1) : index_end_i
    
    index_start <- cumulative_sum(s[i],n) 
    index_end <- cumulative_sum(s[i]+1,n)
    index_set = (index_start +1) : index_end
    
    y1[index_set_i] = y[index_set]
    t1[index_set_i] = t[index_set ]
  }
  
  
  # covariate matrix A1: 
  A0 <- matrix(1, nrow= sum(n1), ncol=2)
  A0[,2]<- t1
  
  
  n0 = 5
  H  <- matrix(0,nrow=n0*sum(n1), ncol=sum(n1))
  H1 <- matrix(0,nrow=sum(n1), ncol=sum(n1))
  H2 <- matrix(0,nrow=sum(n1), ncol=sum(n1))
  H3 <- matrix(0,nrow=sum(n1), ncol=sum(n1))
  H4 <- matrix(0,nrow=sum(n1), ncol=sum(n1))
  H5 <- diag(1,nrow=sum(n1), ncol=sum(n1))
  
  
  ##### H1, H2, H3, H4 ######
  # ith row and jth column block: 
  for (i in 1:N)
  {
    index1_i<- cumulative_sum(i,n1) 
    index2_i<- cumulative_sum(i+1,n1)
    index_set_i = (index1_i +1) : index2_i
    
    
    J1<- matrix(1,nrow=n1[i],ncol=n1[i])
    B1 <- t1[index_set_i]
    
    # For H3 and H4: 
    # diagonal block: 
    H3[ index_set_i,  index_set_i] <- J1
    H4[ index_set_i,  index_set_i] <- B1%*%t(B1)
    
    A1<- matrix(1,nrow=n1[i], ncol=1)
    
    for (j in 1:N)
    {
      
      index1_j<- cumulative_sum(j,n1) 
      index2_j<- cumulative_sum(j+1,n1)
      index_set_j = (index1_j +1) : index2_j
      
      # ith row and jth column block: 
      A2<- matrix(1, nrow=n1[j], ncol=1)
      B2<- t1[index_set_j]
      
      H1[index_set_i, index_set_j]<- (G[i,j]) *A1 %*% t(A2) 
      H2[index_set_i, index_set_j]<- (G[i,j]) *B1 %*% t(B2) 
    }
  }
  
  
  H[1: sum(n1), ] <- H1
  H[(sum(n1) + 1): (2*sum(n1)), ]<- H2
  H[ (2*sum(n1) + 1): (3*sum(n1)), ]<- H3
  H[ (3*sum(n1) + 1): (4*sum(n1)), ]<- H4
  H[ (4*sum(n1) + 1): (5*sum(n1)), ]<- H5
  
  return(list(n=n1,G=G, y=y1, t=t1, H=H, A= A0))
  
}  
################################################################################################################




















# Mimic PLCO data with gentic information 
# Function for generating n, G, G0, t, y: 
generate_data <- function(P, P0, N, J, theta, beta) {
  # Required library for mvrnorm function
  library(MASS)
  
  # P: number of genome-wide SNPs (millions)
  # P0: number of all causal SNPs  (~10,000)
  # N: number of subjects
  # J: maximum of measurements per subject among N subjects 
  # theta: column vector of variance components (sigma^2_g, sigma^2_g*, sigma^2_b0, sigma^2_b1, sigma^2_e)
  # beta: coefficients of fixed effects beta_0, beta_1
  # generate_information: a function to obtain matrices information based on known n, G,G0, t, y 
  
  
  ## Initialize fixed effects: entry age and number of measurements for ith individual
  set.seed(1)
  # Assign different probabilities for measurements from 1 to J for N individuals
  n <- sample(1:J, N, replace = TRUE, prob = c(0.01, 0.02, 0.02, 0.1, 0.15, 0.7))
  # Range of entry age from 55 to 74 years old to mimic the PLCO screening trial prostate cancer data set
  age <- sample(54:74, N, replace = TRUE)
  
  ##  Generate normalized genotypic values
  # Settings for minor allele frequency for all P SNPs 
  f0 <- runif(P, min = 0.01, max = 0.5)
  # Standardized genotypic values for pth SNP of ith individual
  Z <- matrix(0, nrow = N, ncol = P)
  for (p in 1:P) 
  {
    # Count number for pth SNP of ith individual
    x <- sample(0:2, N, replace = TRUE, prob = c((1 - f0[p])^2, 2 * f0[p] * (1 - f0[p]), f0[p]^2))
    # formula to calculate normalized genotypic value z_{ip}
    Z[, p] <- (x - 2 * f0[p]) / sqrt(2 * f0[p] * (1 - f0[p]))
  }
  
  
  ## Generate time effects
  t0 <- matrix(0, nrow = N, ncol = J)
  for (j in 1:J) 
  {
    # screening days for jth measurement of ith individual
    t0[, j] <- sample(((j - 1) * 365 + 1):(j * 365), N, replace = TRUE) / 365
  }
  
  
  
  ## Helper function to calculate the cumulative sum of n up to the (i-1)th element
  cumulative_sum<- function(i,n)
  {
    ifelse(i==1, 0, sum(n[1:(i-1)]) )
  }
  
  # Definition of time effect: estimated age minus 55
  t <- rep(0, sum(n))
  for (i in 1:N) {
    
    start_idx <- cumulative_sum(i, n) + 1
    end_idx <- cumulative_sum(i, n) + n[i]
    
    # time variable = (Estimated age -54)/30 : 
    t[start_idx:end_idx] <- (age[i] + sort(sample(t0[i, ], n[i])) -54)/30
  }
  
  
  
  
  
  # Simulate settings of genetic effects, random effects, and errors
  set.seed(1)
  # sigma^2_g = P0 * sigma^2_alpha,  sigma^2_g* = P0 * sigma^2_eta
  effect_sizes <- MASS::mvrnorm(P0, mu = rep(0,2), Sigma = matrix(c(theta[1]/P0, 0, 0, theta[2]/P0), nrow = 2))
  colnames(effect_sizes) <- c("alpha", "eta")
  alpha <- effect_sizes[, 1]
  eta <- effect_sizes[, 2]
  
  b0 <- rnorm(N, 0, sqrt(theta[3]))
  b1 <- rnorm(N, 0, sqrt(theta[4]))
  e <- rnorm(sum(n), 0, sqrt(theta[5]))
  
  
  # Generate response variable y
  y <- matrix(0, nrow = sum(n), ncol = 1)
  g1 <- matrix(0, nrow = N, ncol = 1)
  g2 <- matrix(0, nrow = N, ncol = 1)
  
  
  
  for (i in 1:N) 
  {
    # For notation simplicity, we make the first P0 SNPs as the causal SNPs
    
    # total genetic effects on baseline
    g1[i] <- Z[i, 1:P0] %*% alpha
    # total genetic effects on slope
    g2[i] <- Z[i, 1:P0] %*% eta
    
    for (j in 1:n[i]) 
    {
      idx <- cumulative_sum(i, n) + j
      y[idx] <- beta[1] + beta[2] * t[idx] + g1[i] + g2[i] * t[idx] + b0[i] + b1[i] * t[idx] + e[idx]
    }
  }
  
  # The genetic relationship matrix calculated using genome-wide variants
  G = (Z%*%t(Z))/P
  
  # The genetic relationship matrix calculated using all causal variants
  Z0 = Z[1:P0,1:P0]
  G0 = (Z0%*%t(Z0))/P0
  
  return(list(n=n, G=G, G0=G0, t=t, y=y))
  
}



###################################################################################################################
# Set a specific seed for reproducibility
set.seed(1)
# Define number of genome-wide common variants, number of causal variants, sample size etc.
#  P represents the number of genome-wide common variants
P <- 1000
# P0 represents the number of causal variants
P0 <- 100
#  N denotes the number of subjects
N <- 100000
#  J denotes maximum of expected measurements per subjects among N subjects
J <- 6
# theta: column vector of variance components = (sigma^2_g, sigma^2_g*, sigma^2_b0, sigma^2_b1, sigma^2_e)
theta = c(2, 2, 2, 2, 0.1)
# beta: coefficients of fixed effects (beta_0 and beta_1)
beta = c(-0.2118, 0.8415)
# Call function to generate data0 with a list of vectors: n, G, G0, t, y.
Data <- generate_data(P = P, P0 = P0, N = N, J = J, theta, beta)
##########################################################################################################


Final_results = Final.result.AI.ReML(initial.par=theta, parts=5, part.result.AI.ReML=part.result.AI.ReML, Data=Data)

print(Final_results)






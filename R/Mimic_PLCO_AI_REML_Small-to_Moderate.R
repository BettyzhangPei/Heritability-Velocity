final.result.AI.ReML <- function(initial.par, data, AI_ReML, f_V, AI_DL)
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
  
  df<- rep(0,n0+8)
  df[1:n0]<- est.par
  df[n0+1]<- result$convergence
  
  df[n0+2]<- est.ratio[1]
  df[n0+3]<- est.ratio[2]
  
  df[n0+4] <- est.var.genetics[1]
  df[n0+5] <- est.var.genetics[2]
  
  df[n0+6]<- est.var.ratio[1]
  df[n0+7]<- est.var.ratio[2]
  df[n0+8]<- inverse.est.AI[2,4]
  
  df = round(df, 3)
  
  df <- data.frame(value = df)
  
  # Set custom column names
  rownames(df) <- c("var.g", "var.g*", "var.b0", "var.b1", "var.e", 
                    "convergence", 
                    "base.heritability","velocity.heritability", 
                    "var.var.g", "var.var.g*", 
                    "var.h1", "var.h2", "cov.g.g*")
  
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
  
  
  # Estimated variance of phenotypic values
  var.ph = as.numeric(sd(y)^2)
  
  
  # Assuming 'par' is a vector and n0 is its length
  n0 = length(par)
  
  # initial setup
  old_par<- par
  
  
  # Function to update parameters
  update_parameters <- function(par) {
    AD <- AI_DL(par, data, f_V)
    new_par <- par - solve(AD$AI) %*% AD$DL
    # Ensure parameters are non-negative
    new_par <- pmax(new_par, var.ph*(1e-6))  # ensure positivity 
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





# Function for computing covariance matrix  V  
  f_V <- function(par, H, n)
  {
    
    K = length(par)
    
    V= matrix(0, nrow= sum(n), ncol=sum(n))
    
    for (i in 1:K)
    { 
      H_temp =   H[ ((i-1)*sum(n) + 1): (i*sum(n)), ] 
    
      V = V+ par[i]* H_temp 
    }
    
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



################### For real data analysis with known data0 as a list of column vectors n, G, G0, t, y ##################
# call function to obtain information n, G, G0, t, y, and matrices A and H based on known n, G, G0, t, y
generate_information<- function(n, G, t, y)
{
  # covariate matrix A: 
  A <- matrix(1, nrow= sum(n), ncol=2)
  A[,2]<- t
  
  ## Helper function to calculate the cumulative sum of n up to the (i-1)th element
  cumulative_sum<- function(i,n)
  {
    ifelse(i==1, 0, sum(n[1:(i-1)]) )
  }
  
  n0 = 5
  
  H  <- matrix(0,nrow=n0*sum(n), ncol=sum(n))
  H1 <- matrix(0,nrow=sum(n), ncol=sum(n))
  H2 <- matrix(0,nrow=sum(n), ncol=sum(n))
  H3 <- matrix(0,nrow=sum(n), ncol=sum(n))
  H4 <- matrix(0,nrow=sum(n), ncol=sum(n))
  H5 <- diag(1,nrow=sum(n), ncol=sum(n))
  
  
  ##### H1, H2, H3, H4 ######
  # ith row and jth column block: 
  for (i in 1:N)
  {
    index1_i<- cumulative_sum(i,n) 
    index2_i<- cumulative_sum(i+1,n)
    index_set_i = (index1_i +1) : index2_i
    
    
    J1<- matrix(1,nrow=n[i],ncol=n[i])
    B1 <- t[index_set_i]
    
    # For H3 and H4: 
    # diagonal block: 
    H3[ index_set_i,  index_set_i] <- J1
    H4[ index_set_i,  index_set_i] <- B1%*%t(B1)
    
    A1<- matrix(1,nrow=n[i], ncol=1)
    
    for (j in 1:N)
    {
      
      index1_j<- cumulative_sum(j,n) 
      index2_j<- cumulative_sum(j+1,n)
      index_set_j = (index1_j +1) : index2_j
      
      # ith row and jth column block: 
      A2<- matrix(1,nrow=n[j], ncol=1)
      B2<- t[index_set_j]
      
      H1[index_set_i, index_set_j]<- (G[i,j]) *A1 %*% t(A2) 
      H2[index_set_i, index_set_j]<- (G[i,j]) *B1 %*% t(B2) 
    }
  }
  
  
  H[1: sum(n), ] <- H1
  H[(sum(n) + 1): (2*sum(n)), ]<- H2
  H[ (2*sum(n) + 1): (3*sum(n)), ]<- H3
  H[ (3*sum(n) + 1): (4*sum(n)), ]<- H4
  H[ (4*sum(n) + 1): (5*sum(n)), ]<- H5
  
  
  return(list(G=G, t=t, n=n, y=y, H=H, A=A))
  
  
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
N <- 2000
#  J denotes maximum of expected measurements per subjects among N subjects
J <- 6
# theta: column vector of variance components = (sigma^2_g, sigma^2_g*, sigma^2_b0, sigma^2_b1, sigma^2_e)
theta = c(2, 2, 2, 2, 0.1)
# beta: coefficients of fixed effects (beta_0 and beta_1)
beta = c(-0.2118, 0.8415)
# Call function to generate data0 with a list of vectors: n, G, G0, t, y.
data0 <- generate_data(P = P, P0 = P0, N = N, J = J, theta, beta)
##########################################################################################################



################### For real data analysis with known data0 as a list of column vectors n, G, G0, t, y ##################
# call function to obtain information n, G, G0, t, y, and matrices A and H based on known n, G, G0, t, y
data = generate_information(n = data0$n, G = data0$G, t= data0$t, y=data0$y)
# Perform AI-ReML algorithm for estimation of unknown variance parameters, two heritability metrics and their standard errors 
# For instance, we choose an arbitrary input for unknown variance components for AI-ReML algorithm


final_results <- final.result.AI.ReML(initial.par=theta, data, AI_ReML, f_V, AI_DL)

print(final_results)



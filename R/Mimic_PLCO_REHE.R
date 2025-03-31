
#install.packages("quadprog")
library(quadprog)  # KKT optimization process for minimizing the convex quadratic function

library(Rcpp)
# Access the C++ Code:
# sourceCpp("~/Desktop/C++/sumS3.cpp")
# sourceCpp("~/Desktop/C++/SumD_c.cpp")
  


final.result.REHE <- function(par, REHE_results, data0)
{
  # par: unknown parameters = (beta0, beta1, var.g, var.g*, var.b0, var.b1, var.e)
  # REHE_results: function to implement the REHE method to obtain the estimates of par
  
  
  ### data: generated data information with a list of matrices:  
  # G : genetic relationship matrix 
  # t : the sum(n) x 1 column vector of time effect
  # n : the N x 1 column vector represents total number of measurements for each subjects
  # y : the  sum(n) x 1 column vector of response 
  # H : composite matrix for AI matrix and DL matrix in the iteration process of AI-ReML algorithm
  # A : covariate matrix for fixed effects
  
  
  par= as.numeric(REHE_results(data0))
  
  Rep = 1000
  
  result = matrix(0, nrow = 9, ncol=Rep)
  
  for (a in 1: Rep)
  { 
    
    ###################################################################################################################
    # Set a specific seed for reproducibility
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
    theta = as.numeric(par[3:7])
    # beta: coefficients of fixed effects (beta_0 and beta_1)
    beta  = as.numeric(par[1:2])
    # Call function to generate data0 with a list of vectors: n, G, G0, t, y.
    data <- generate_data(a=a,P = P, P0 = P0, N = N, J = J, theta, beta)
    ####################################################################################################################
    
    # Call function to estimate unknown parameters using REHE method: 
    result[,a] = as.numeric(REHE_results(data))
    
  }
  
  para.boot.REHE = matrix(0, nrow=9, ncol =3)
  para.boot.REHE[,1] = as.numeric(par)
  for(i in 1:9)
  {
    # Empirical SE for parameters
    para.boot.REHE[i, 2] = sd(result[i,])
    
    # Scaled MAD of parameters 
    para.boot.REHE[i, 3] = mad(result[i,])
  }
  
  
  para.boot.REHE = round( para.boot.REHE, 3)
  para.boot.REHE = data.frame(para.boot.REHE)
  rownames(para.boot.REHE) = c("beta0", "beta1", "var.g", "var.g*", "var.b0", "var.b1", "var.e", "base.heritability","velocity.heritability" )
  
  colnames(para.boot.REHE) = c("Parameters", "Emp SE", "Scaled MAD")
  
  return(para.boot.REHE)
  
}





REHE_results <- function(data0)
{
  # REHE_results: function to implement the REHE method to obtain the estimates of par
  
  
  ### data0: generated data information with a list of matrices:  
  # G : genetic relationship matrix 
  # t : the sum(n) x 1 column vector of time effect
  # n : the N x 1 column vector represents total number of measurements for each subjects
  # y : the  sum(n) x 1 column vector of response 
  
  library(Rcpp)
  # Access the C++ Code:
  sourceCpp("~/Desktop/C++/sumS3.cpp")
  sourceCpp("~/Desktop/C++/SumD_c.cpp")
  
  G = data0$G
  t = data0$t
  n = data0$n
  y = data0$y
  N = dim(G)[1]
  
  # covariate matrix A for fixed effects: 
  A <- matrix(1, nrow= sum(n), ncol=2)
  A[,2]<- t
  
  # Estimate coefficients of fixed effects and then center the response 
  # estimated fixed effects for beta0 and beta1:  
  est.beta<- solve(t(A)%*%A)%*%t(A)%*%y
  
  # resid_y with mean around zero: 
  resid_y <-  y - A%*% est.beta
  
  ##### Helper function to calculate the cumulative sum of n up to the (i-1)th element
  # Pre-compute the cumulative sums
  cum_n <- numeric(N + 1)
  cum_n[1] = 0
  
  for (i in 2:(N + 1)) {
    cum_n[i] <- cum_n[i - 1] + n[i - 1]
  }
  
  
  ################loss1##############
  Y1<- resid_y^2
  A1<- matrix(1, nrow=sum(n), ncol=5)
  A1[,4]<- t^2
  
  for (i in 1:N)
  {
    index_start <- cum_n[i] + 1
    index_end <- cum_n[i + 1]
    index_set = index_start : index_end
    
    A1[index_set, 1] <- G[i,i]
    A1[index_set, 2] <- G[i,i]* A1[index_set, 4]
  }
  
  
  ######################loss2##############
  # Calculate total pairs for pre-allocation
  total_pairs <- sum(sapply(n, function(x) choose(x, 2)))
  # total_pairs = 185368
  # length(which(n==1)) =210, then we know that there are 210 males with exactly 1 measurement of PSA level. 
  s<- which(n>1)
  
  # For Y2 and A2 in the Loss function: 
  Y2<- numeric(total_pairs)
  A2<- matrix(1, nrow=total_pairs, ncol=5)
  A2[,5] <- 0
  A2_1<- numeric(total_pairs)
  A2_2<- numeric(total_pairs)
  A2_4<- numeric(total_pairs)
  
  
  
  #### Method: Matrix form ###
  # Loop through each person with at least 2 measurements of PSA level: 
  idx <- 1
  for (i in s)
  {
    index_start <- cum_n[i] + 1
    index_end <- cum_n[i+1] 
    index_set = index_start : index_end
    
    
    new_resid_y = resid_y[index_set]
    new_t =  t[index_set]
    # new_yt =  yt[index_set]
    #new_t2 = t2[index_set]
    
    t_indices <- 1:n[i]
    
    comb <- combn(t_indices, 2)   # Get all pairwise combinations of indices
    
    
    # Compute the products and store them in the pre-allocated vector
    Y2[idx:(idx + ncol(comb) - 1)] <- new_resid_y[comb[1, ]] * new_resid_y[comb[2, ]]
    A2_1[idx:(idx + ncol(comb) - 1)] <- G[i,i]
    A2_2[idx:(idx + ncol(comb) - 1)] <- G[i,i]* (new_t[comb[1, ]] * new_t[comb[2, ]])
    A2_4[idx:(idx + ncol(comb) - 1)] <- new_t[comb[1, ]] * new_t[comb[2, ]]
    
    idx <- idx + ncol(comb)
  }
  
  
  A2[,1]<- A2_1
  A2[,2]<- A2_2
  A2[,4]<- A2_4
  
  
  
  # For D in the Gradient function: 
  
  ##### Gradient function of the loss function #####
  # Initialize vectors with appropriate sizes for efficiency: 
  
  # For D1 and c1:
  D1 = matrix(0, 5, 5)
  c1 = rep(0,5)
  
  
  
  Y11_1<- numeric(sum(n))
  Y12_1<- numeric(sum(n))
  Y13_1<- numeric(sum(n))
  Y14_1<- numeric(sum(n))
  Y15_1<- numeric(sum(n))
  Y22_1<- numeric(sum(n))
  Y23_1<- numeric(sum(n))
  Y24_1<- numeric(sum(n))
  Y25_1<- numeric(sum(n))
  Y33_1<- numeric(sum(n))
  Y34_1<- numeric(sum(n))
  Y35_1<- numeric(sum(n))
  Y44_1<- numeric(sum(n))
  Y45_1<- numeric(sum(n))
  Y55_1<- numeric(sum(n))
  c1_1 <-  numeric(sum(n))
  c2_1 <-  numeric(sum(n))
  c3_1 <-  numeric(sum(n))
  c4_1 <-  numeric(sum(n))
  c5_1 <-  numeric(sum(n))
  
  ##### Loop through each person
  for (i in 1:N) {
    index_start = cum_n[i] + 1
    index_end =  cum_n[i+1]
    index_set = index_start : index_end
    
    Y11_1[index_set] <-   G[i,i]^2
    Y12_1[index_set] <-  (G[i, i]^2) * (t[index_set]^2)
    Y13_1[index_set] <-   G[i,i]
    Y14_1[index_set] <-  G[i, i] * (t[index_set]^2)
    Y15_1[index_set] <-   G[i,i]
    
    Y22_1[index_set] <-  (G[i, i]^2) * (t[index_set]^4)
    Y23_1[index_set] <-  (G[i, i]) * (t[index_set]^2)
    Y24_1[index_set] <-  G[i,i] * (t[index_set]^4)
    Y25_1[index_set] <-  (G[i, i]) * (t[index_set]^2)
    
    Y33_1[index_set] <- 1
    Y34_1[index_set] <- (t[index_set]^2)
    Y35_1[index_set] <- 1
    
    Y44_1[index_set] <- (t[index_set]^4)
    Y45_1[index_set] <- (t[index_set]^2)
    
    Y55_1[index_set] <- 1
    
    c1_1[index_set]  <-  G[i, i] * (resid_y[index_set]^2)
    c2_1[index_set]  <-  G[i, i] * (resid_y[index_set]^2) * (t[index_set]^2)
    c3_1[index_set]  <-   (resid_y[index_set]^2)
    c4_1[index_set]  <-   (resid_y[index_set]^2) * (t[index_set]^2)
    c5_1[index_set]  <-   (resid_y[index_set]^2)
    
  }
  
  
  D1 <- matrix(c(
    sum(Y11_1), sum(Y12_1),  sum(Y13_1), sum(Y14_1),  sum(Y15_1),
    sum(Y12_1), sum(Y22_1),  sum(Y23_1), sum(Y24_1),  sum(Y25_1),
    sum(Y13_1), sum(Y23_1),  sum(Y33_1), sum(Y34_1),  sum(Y35_1),
    sum(Y14_1), sum(Y24_1),  sum(Y34_1), sum(Y44_1),  sum(Y45_1),
    sum(Y15_1), sum(Y25_1),  sum(Y35_1), sum(Y45_1),  sum(Y55_1)
  ), nrow = 5, ncol = 5, byrow = TRUE)
  
  c1 <- c( sum(c1_1), sum(c2_1),sum(c3_1), sum(c4_1), sum(c5_1))
  
  
  # For D2, c2: 
  D2 = matrix(0, 5, 5)
  c2 = rep(0,5)
  
  Y11_2<- numeric(total_pairs)
  Y12_2<- numeric(total_pairs)
  Y13_2<- numeric(total_pairs)
  Y14_2<- numeric(total_pairs)
  Y22_2<- numeric(total_pairs)
  Y23_2<- numeric(total_pairs)
  Y24_2<- numeric(total_pairs)
  Y33_2<- numeric(total_pairs)
  Y34_2<- numeric(total_pairs)
  Y44_2<- numeric(total_pairs)
  c1_2<- numeric(total_pairs)
  c2_2<- numeric(total_pairs)
  c3_2<- numeric(total_pairs)
  c4_2<- numeric(total_pairs)
  
  
  yt = resid_y*t
  t2 = t^2
  
  
  #### Method: Matrix form ###
  # Loop through each person with at least 2 measurements of PSA level: 
  idx <- 1
  for (i in s)
  {
    index_start <- cum_n[i] +1
    index_end <- cum_n[i+1]
    index_set = index_start : index_end
    
    
    new_resid_y = resid_y[index_set]
    new_t =  t[index_set]
    new_yt =  yt[index_set]
    new_t2 = t2[index_set]
    
    t_indices <- 1:n[i]
    
    comb <- combn(t_indices, 2)   # Get all pairwise combinations of indices
    
    Y11_2 [idx:(idx + ncol(comb) - 1)] <- (G[i,i]^2)
    Y12_2[idx:(idx + ncol(comb) - 1)] <- (G[i,i]^2)* (new_t[comb[1, ]] * new_t[comb[2, ]])
    Y13_2[idx:(idx + ncol(comb) - 1)] <- G[i,i]
    Y14_2[idx:(idx + ncol(comb) - 1)] <- (G[i,i])* (new_t[comb[1, ]] * new_t[comb[2, ]])
    
    Y22_2[idx:(idx + ncol(comb) - 1)] <- (G[i,i]^2)*  ((new_t2[comb[1, ]]* new_t2[comb[2, ]]))
    Y23_2[idx:(idx + ncol(comb) - 1)] <- (G[i,i])* (new_t[comb[1, ]] * new_t[comb[2, ]])
    Y24_2[idx:(idx + ncol(comb) - 1)] <- G[i,i] * (new_t2[comb[1, ]] * new_t2[comb[2, ]])
    
    Y33_2[idx:(idx + ncol(comb) - 1)] <- 1
    Y34_2[idx:(idx + ncol(comb) - 1)] <-  (new_t[comb[1, ]] * new_t[comb[2, ]])
    Y44_2[idx:(idx + ncol(comb) - 1)] <-  (new_t2[comb[1, ]] * new_t2[comb[2, ]])
    
    c1_2[idx:(idx + ncol(comb) - 1)] <-  G[i,i]* (new_resid_y[comb[1, ]]* new_resid_y[comb[2, ]])
    c2_2[idx:(idx + ncol(comb) - 1)] <-  G[i,i]* (new_yt[comb[1, ]]* new_yt[comb[2, ]])
    c3_2[idx:(idx + ncol(comb) - 1)] <-  new_resid_y[comb[1, ]]* new_resid_y[comb[2, ]]
    c4_2[idx:(idx + ncol(comb) - 1)] <-  (new_yt[comb[1, ]]* new_yt[comb[2, ]])
    
    idx <- idx + ncol(comb)
  }
  
  
  D2 <- matrix(c(
    sum(Y11_2), sum(Y12_2),  sum(Y13_2), sum(Y14_2),  0,
    sum(Y12_2), sum(Y22_2),  sum(Y23_2), sum(Y24_2),  0,
    sum(Y13_2), sum(Y23_2),  sum(Y33_2), sum(Y34_2),  0,
    sum(Y14_2), sum(Y24_2),  sum(Y34_2), sum(Y44_2), 0,
    0, 0 , 0, 0,  0
  ), nrow = 5, ncol = 5, byrow = TRUE)
  
  c2 <- c( sum(c1_2), sum(c2_2),sum(c3_2), sum(c4_2), 0)
  
  # For D3, c3: 
  D3 = matrix(0, 5, 5)
  c3 = rep(0,5)
  
  
  #sourceCpp("SumD_c.cpp")  # desktop 
  
  sum_D_c = SumD_c(G,resid_y, n, t)
  Y11_3 = sum_D_c[1]
  Y12_3 = sum_D_c[2]
  Y22_3 = sum_D_c[3]
  c1_3 = sum_D_c[4]
  c2_3 = sum_D_c[5]
  
  
  D3 <- matrix(c(
    Y11_3, Y12_3,  0, 0,  0,
    Y12_3, Y22_3,  0, 0,  0,
    0, 0, 0, 0,  0,
    0, 0, 0, 0,  0, 
    0, 0, 0, 0,  0), nrow = 5, ncol = 5, byrow = TRUE)
  
  
  c3 <- c(c1_3, c2_3,0,0, 0)
  
  D = D1 + 2*D2 + D3
  c = c1+ 2*c2 +c3  
  
  
  # KKT optimization process for minimizing the convex quadratic function
  #install.packages("quadprog")
  library(quadprog)
  
  ## Alternative way to express the loss function
  F1<- function(x)
  {
    0.5*t(x)%*%(2*D)%*%x + t(-2*c)%*%x
  } 
  
  
  # Solve the quadratic programming problem
  result_KKT <- solve.QP(Dmat=2*D, dvec= (2*c), Amat=diag(5), bvec= rep(0,5), meq=0)
  
  
  # Optimal solution
  x_opt <- result_KKT$solution
  
  df = rep(0, nrow=9, ncol=1)
  df[1:2] = est.beta
  df[3:7] = x_opt
  df[8] = x_opt[1]/(x_opt[1] + x_opt[3])
  df[9] = x_opt[2]/(x_opt[2] + x_opt[4])
  
  df = round(df, 3)
  # df <- data.frame(df)
  
  # Set custom column names
  # rownames(df) <- c("beta0", "beta1", "var.g", "var.g*", "var.b0", "var.b1", "var.e", "base.heritability","velocity.heritability" )
  return(df)
}


# Mimic PLCO data with gentic information 
# Function for generating n, G, G0, t, y:  
generate_data <- function(a=1,P, P0, N, J, theta, beta) {
  # Required library for mvrnorm function
  library(MASS)
  
  # P: number of genome-wide SNPs (millions)
  # P0: number of all causal SNPs  (~10,000)
  # N: number of subjects
  # J: maximum of measurements per subject among N subjects 
  # theta: column vector of variance components (sigma^2_g, sigma^2_g*, sigma^2_b0, sigma^2_b1, sigma^2_e)
  # beta: coefficients of fixed effects beta_0, beta_1
  # generate_information: a function to obtain matrices information based on known n, G, G0, t, y 
  
  
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
  set.seed(a)
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



# ### Application to a real data set 
# ###################################################################################################################
# # Set a specific seed for reproducibility
# set.seed(1)
# # Define number of genome-wide common variants, number of causal variants, sample size etc.
# #  P represents the number of genome-wide common variants
# P <- 1000
# # P0 represents the number of causal variants
# P0 <- 100
# #  N denotes the number of subjects
# N <- 2000
# #  J denotes maximum of expected measurements per subjects among N subjects
# J <- 6
# # theta: column vector of variance components = (sigma^2_g, sigma^2_g*, sigma^2_b0, sigma^2_b1, sigma^2_e)
# theta = c(2, 2, 2, 2, 0.1)
# # beta: coefficients of fixed effects (beta_0 and beta_1)
# beta = c(-0.2118, 0.8415)
# # Call function to generate data0 with a list of vectors: n, G, G0, t, y.
# data0 <- generate_data(a=1,P = P, P0 = P0, N = N, J = J, theta, beta)
# print(REHE_results(data0))
# 
# print(final.result.REHE(par, REHE_results, data0))


# Heritability-Velocity: Estimating the Heritability of Velocity in a Longitudinal Phenotypic Trait: Genetic Insights into PSA Velocity in Prostate Cancer-Free Individuals

Approaches for heritability analysis typically focus on analyzing cross-sectional phenotype, with less attention given to longitudinal phenotypes. This R code introduces a novel mixed model that incorporates genome-wide common variants to disentangle joint genetic effects on both baseline and velocity in a longitudinal phenpytpe. We employ two approaches: the average information restricted maximum likelihood (AI-REML) algorithm and the restricted Haseman-Elston (REHE) regression method. 

The AI-REML algorithm becomes computationally prohibitive for large-scale studies due to inverting the high dimensionality of covariance matrix. 
To overcome this challenge, we introduce a novel partitioning and combining strategy: first evenly partitioning subjects into multiple subgroups, estimating parameters separately via AI-REML, and then performing a meta-analysis using a left-truncated log-likelihood method to derive final estimates and their standard errors. However, this method sacrifices efficiency due to the loss of non-diagonal information in the covariance matrix. Alternatively, the REHE method with a fast algorithm that significantly reduces computational costs, leverages pairwise information to efficiently avoid the challenges of covariance matrix inversion in large-scale studies.  

Utilizing those two approaches, this code provides examples to estimate heritability metrics on both baseline and velocity of a longitudinal phenotype by mimicing prostate-specific antigen (PSA) trajectories among prostate cancer-free individuals in the Prostate, Lung, Colorectal, and Ovarian (PLCO) Cancer Screening Trial.  Simulation studies demonstrated that both methods yield unbiased estimates of joint heritability metrics but with reduced efficiency using genome-wide variants. Also, AI-REML algorithm outperforms REHE method when the dataset requires partitioning into only a few groups. 

Building on the proposed mixed model with jointly individual-level genetic effects on both baseline and velocity of a longitudinal phenotype, our R code employs the AI-ReML algorithm and REHE method to estimate two heritability metrics and their standard errors. This R code estimates:
- All variance components for genetic- and subject- specific random effects and residual error term.  
- Two heritability metrics on both baseline and velocity of a longitudinal phenotype.
- Standard errors for estimated variance components and two estimated heritability metrics. 


# Usage Examples
For real data analysis, we use `PLCO.R`. Among `PLCO.R`: 
- `generate_data` function
- `generate_information` function
- `AI-ReML` function
The `generate_data` function mimic the PLCO datasets to provide a list of observed data: n, Z, t and y. 
The `generate_information` function provides a list of obaserved data and required matrices: n, Z, t, y, A, S, G, W, and H.
The `AI-ReML` function provides estimates for variance components in a linear mixed model by integrating genetic effects in a longitudinal phenotype via AI-ReML algorithm based on inputed dataset. 
For instance: 
```r
###################################################################################################################
# Set a specific seed for reproducibility
set.seed(1)
# Define number of SNPs, sample size, etc.
#  P: number of SNPs
P <- 100
#  N: number of subjects
N <- 2000
#  J: maximum of expected measurements per person among N individuals
J <- 6
# theta: column vector of variance components = (sigma^2_alpha, sigma^2_eta, sigma^2_b0, sigma^2_b1, sigma^2_e)
theta = c(0.005, 0.002, 0.5, 0.2, 0.2)
# xi: coefficients of fixed effects beta_0, beta_1, mu_alpha, mu_eta
xi = c(0.5, 0.05, 0, 0)
# Call function to generate data0 with a list of vectors: n,Z,t,y
data0 <- generate_data(P = P, N = N, J = J, theta, xi)
####################################################################################################################

################### For real data analysis with known data0 as a list of column vectors n, Z, t,y ##################
# call function to obtain information n,Z,t,y, and matrices A, S, G, W,H based on known n,Z,t,y 
data = generate_information(n= data0$n, Z= data0$Z, t= data0$t, y=data0$y)
# Perform AI-ReML algorithm for estimation of unknown variance parameters
# For instance, we choose an arbitrary input for unknown variance components for AI-ReML algorithm
initial.par = theta
result <- AI_ReML(par = initial.par, l_ReML,  maxit = 1000, maxtol = 1e-4, data = data, f_V, AI_DL, woodbury_inverse_V)
# print results from AI-ReML algorithm
print(result)
# Estimated variance components = (sigma^2_alpha, sigma^2_eta, sigma^2_b0, sigma^2_b1, sigma^2_e)
est.par= as.matrix(result$par, ncol=1)
# Estimated two heritability metrics and their standard errors 
####################################################################################################################
```
where data0 has to be structured as a list of vectors n, Z, t, y. data has to be structured as a list of vectors and matrices n, Z, t, y, A, S, G, W and H. 
- n: A N x 1 column vector represents total number of measurements for each subject. 
- Z: A N x P matrix with standardized genotypic values.
- t: A sum(n) x 1 column vector of the time variable.
- y: A sum(n) x 1 column vector of the longitudinal response.
- A: A sum(n) x 2 covariate matrix of fixed effects beta_0 and beta_1.
- S: A sum(n) x 4 covariate matrix of fixed effects xi = (beta_0, beta1, mu_alpha, mu_eta).
- G: A N x N matrix as P times genetic relationship matrix.
- W: A sum(n) x (2P+2N) covariate matrix of random effects.
- H: composite matrix for AI matrix and DL matrix.

# Usage Notes
1. We recommend transforming response (e.g., using a log transformation) to approximate a normal distribution before applying our functions.
2. We recommend transforming time variable (e.g., the temporal effect t is defined using approximately min-max normalization based on age in years to
scale its range between 0 and 1. In specific, time variable is calculated as the age at each screening visit, adjusted by subtracting
54 and dividing by 30 for PLCO prostate cancer dataset) to control the scale of time variable before applying our functions. 
3. Our calculations assume there is unbalanced longitudinal data structure.
4. For `sim_II.R` and `sim_II_genetic_effects.R`, both of them are for simulation studies for a specific scenario with 1,000 repetitions, we conducted them on Biowulf by setting value of a from 1 to 1000.  
   
# References
Lynch, Michael and Walsh, Bruce (1998) Genetics and analysis of quantitative traits. Sinauer Sunderland. 

Haseman, Joseph K and Elston, Robert C (1972) The investigation of linkage between a quantitative trait and a marker locus. Behavior Genetics, 2(1): 3-19. 

Johnson, DL and Thompson, Robin (1995) Restricted maximum likelihood estimation of variance components for univariate animal models using sparse matrix techniques and average information.  Journal of dairy science, 78(2): 449-456.

Yang, Jian and Lee, S Hong and Goddard, Michael E and Visscher, Peter M (2011) GCTA: a tool for genome-wide complex trait analysis. The American Journal of Human Genetics, 88(1): 76-82.

Yue, Kun and Ma, Jing and Thornton, Timothy and Shojaie, Ali (2021) REHE: Fast variance components estimation for linear mixed models. Genetic Epidemiology, 45(8): 891-905.

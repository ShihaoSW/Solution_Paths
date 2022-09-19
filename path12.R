###### high dimensional model selection

library(MASS)
library(Matrix)
library(picasso)
library(glmnet)
source("functions.R")

###### set seed ########
## The following commented codes are for Monte Carlo repetitions
## The seeds are (1,2,...,100)
# seed <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# set.seed(seed)


####### unchanged parameters
p = 5000
s = 50
n = ceiling(2 * s * log(p))

set.seed(seed)

## 11a
# covariance for x
rho = 0
S_x = diag(p)
for (i1 in c(1:p)){
  for (i2 in c(1:p)){
    # choice 1: exp decay
    S_x[i1,i2] = rho ^ (abs(i1 - i2))
    # choice 2: constant cor
    #if (i1 != i2) S_x[i1,i2] = rho     
  }
}

# signal and noise
betamin = 0.1
sigma_e = 0.5

# decide lambda values in picasso (scad)
lambda_list_scad <- 10/(1.2 ^ seq(0, 99))

# decide lambda values for glmnet (lasso)
lambda_list_lasso <- 10/(1.2 ^ seq(0, 99))

# decide lambda values for picasso (mcp)
lambda_list_mcp <- 10/(1.2 ^ seq(0, 99))

# decide pi values for iht
lambda_list_iht = seq(2,100,2)

## generate data

data_i <- data_lm(n = n, mu_x = rep(0,p), S_x = S_x,
                  s = s, betamin = betamin, sigma_e = sigma_e)

X = data_i$X
Y = data_i$Y  
beta = data_i$beta
epsilon = data_i$epsilon
supp_true = data_i$supp_true

## optimization with nonconcave penalty
ptm <- proc.time()
scad_result = scad_select(X = X, Y = Y, supp_true = supp_true, lambda_list = lambda_list_scad)
scad_TPR = scad_result$TPR
scad_FDR = scad_result$FDR
scad_TPR_cv = scad_result$TPR_cv
scad_FDR_cv= scad_result$FDR_cv

## optimization with convex penalty
lasso_result = lasso_select(X = X, Y = Y, supp_true = supp_true, lambda_list = lambda_list_lasso)
lasso_TPR = lasso_result$TPR
lasso_FDR = lasso_result$FDR
lasso_TPR_cv = lasso_result$TPR_cv
lasso_FDR_cv = lasso_result$FDR_cv

## optimization with minimax concave penalty
mcp_result = mcp_select(X = X, Y = Y, supp_true = supp_true, lambda_list = lambda_list_mcp)
mcp_TPR = mcp_result$TPR
mcp_FDR = mcp_result$FDR
mcp_TPR_cv = mcp_result$TPR_cv
mcp_FDR_cv = mcp_result$FDR_cv

## optimization with iterative hard thresholding
iht_result = iht_m_select(X = X, Y = Y, supp_true = supp_true, pi_list = lambda_list_iht)
iht_TPR= iht_result$TPR
iht_FDR = iht_result$FDR
iht_TPR_cv = iht_result$TPR_cv
iht_FDR_cv = iht_result$FDR_cv
proc.time() - ptm
## save rdata
save(X, Y, beta, epsilon, supp_true, 
     scad_TPR, scad_FDR, scad_TPR_cv, scad_FDR_cv, 
     lasso_TPR, lasso_FDR, lasso_TPR_cv, lasso_FDR_cv,
     mcp_TPR, mcp_FDR, mcp_TPR_cv, mcp_FDR_cv,
     iht_TPR, iht_FDR, iht_TPR_cv, iht_FDR_cv, file = paste0("./11a/result",seed,".RData"))




set.seed(seed)


## 11b
# covariance for x
rho = 0.5
S_x = diag(p)
for (i1 in c(1:p)){
  for (i2 in c(1:p)){
    # choice 1: exp decay
    S_x[i1,i2] = rho ^ (abs(i1 - i2))
    # choice 2: constant cor
    #if (i1 != i2) S_x[i1,i2] = rho     
  }
}

# signal and noise
betamin = 0.1
sigma_e = 0.5

# decide lambda values in picasso (scad)
lambda_list_scad <- 1.5/(1.2 ^ seq(0, 99))

# decide lambda values for glmnet (lasso)
lambda_list_lasso <- 1.5/(1.2 ^ seq(0, 99))

# decide lambda values for picasso (mcp)
lambda_list_mcp <- 1.5/(1.2 ^ seq(0, 99))

# decide pi values for iht
lambda_list_iht = seq(2,100,2)

## generate data
data_i <- data_lm(n = n, mu_x = rep(0,p), S_x = S_x,
                  s = s, betamin = betamin, sigma_e = sigma_e)
X = data_i$X
Y = data_i$Y  
beta = data_i$beta
epsilon = data_i$epsilon
supp_true = data_i$supp_true

## optimization with nonconcave penalty
scad_result = scad_select(X = X, Y = Y, supp_true = supp_true, lambda_list = lambda_list_scad)
scad_TPR = scad_result$TPR
scad_FDR = scad_result$FDR
scad_TPR_cv = scad_result$TPR_cv
scad_FDR_cv= scad_result$FDR_cv

## optimization with convex penalty
lasso_result = lasso_select(X = X, Y = Y, supp_true = supp_true, lambda_list = lambda_list_lasso)
lasso_TPR = lasso_result$TPR
lasso_FDR = lasso_result$FDR
lasso_TPR_cv = lasso_result$TPR_cv
lasso_FDR_cv = lasso_result$FDR_cv

## optimization with minimax concave penalty
mcp_result = mcp_select(X = X, Y = Y, supp_true = supp_true, lambda_list = lambda_list_mcp)
mcp_TPR = mcp_result$TPR
mcp_FDR = mcp_result$FDR
mcp_TPR_cv = mcp_result$TPR_cv
mcp_FDR_cv = mcp_result$FDR_cv

## optimization with iterative hard thresholding
iht_result = iht_m_select(X = X, Y = Y, supp_true = supp_true, pi_list = lambda_list_iht)
iht_TPR= iht_result$TPR
iht_FDR = iht_result$FDR
iht_TPR_cv = iht_result$TPR_cv
iht_FDR_cv = iht_result$FDR_cv

## save rdata
save(X, Y, beta, epsilon, supp_true, 
     scad_TPR, scad_FDR, scad_TPR_cv, scad_FDR_cv, 
     lasso_TPR, lasso_FDR, lasso_TPR_cv, lasso_FDR_cv,
     mcp_TPR, mcp_FDR, mcp_TPR_cv, mcp_FDR_cv,
     iht_TPR, iht_FDR, iht_TPR_cv, iht_FDR_cv, file = paste0("./11b/result",seed,".RData"))





set.seed(seed)

## 11c
# covariance for x
rho = 0.8
S_x = diag(p)
for (i1 in c(1:p)){
  for (i2 in c(1:p)){
    # choice 1: exp decay
    S_x[i1,i2] = rho ^ (abs(i1 - i2))
    # choice 2: constant cor
    #if (i1 != i2) S_x[i1,i2] = rho     
  }
}

# signal and noise
betamin = 0.1
sigma_e = 0.5

# decide lambda values in picasso (scad)
lambda_list_scad <- 1.5/(1.2 ^ seq(0, 99))

# decide lambda values for glmnet (lasso)
lambda_list_lasso <- 1.5/(1.2 ^ seq(0, 99))

# decide lambda values for picasso (mcp)
lambda_list_mcp <- 1.5/(1.2 ^ seq(0, 99))

# decide pi values for iht
lambda_list_iht = seq(2,100,2)

## generate data
data_i <- data_lm(n = n, mu_x = rep(0,p), S_x = S_x,
                  s = s, betamin = betamin, sigma_e = sigma_e)
X = data_i$X
Y = data_i$Y  
beta = data_i$beta
epsilon = data_i$epsilon
supp_true = data_i$supp_true

## optimization with nonconcave penalty
scad_result = scad_select(X = X, Y = Y, supp_true = supp_true, lambda_list = lambda_list_scad)
scad_TPR = scad_result$TPR
scad_FDR = scad_result$FDR
scad_TPR_cv = scad_result$TPR_cv
scad_FDR_cv= scad_result$FDR_cv

## optimization with convex penalty
lasso_result = lasso_select(X = X, Y = Y, supp_true = supp_true, lambda_list = lambda_list_lasso)
lasso_TPR = lasso_result$TPR
lasso_FDR = lasso_result$FDR
lasso_TPR_cv = lasso_result$TPR_cv
lasso_FDR_cv = lasso_result$FDR_cv

## optimization with minimax concave penalty
mcp_result = mcp_select(X = X, Y = Y, supp_true = supp_true, lambda_list = lambda_list_mcp)
mcp_TPR = mcp_result$TPR
mcp_FDR = mcp_result$FDR
mcp_TPR_cv = mcp_result$TPR_cv
mcp_FDR_cv = mcp_result$FDR_cv

## optimization with iterative hard thresholding
iht_result = iht_m_select(X = X, Y = Y, supp_true = supp_true, pi_list = lambda_list_iht)
iht_TPR= iht_result$TPR
iht_FDR = iht_result$FDR
iht_TPR_cv = iht_result$TPR_cv
iht_FDR_cv = iht_result$FDR_cv

## save rdata
save(X, Y, beta, epsilon, supp_true, 
     scad_TPR, scad_FDR, scad_TPR_cv, scad_FDR_cv, 
     lasso_TPR, lasso_FDR, lasso_TPR_cv, lasso_FDR_cv,
     mcp_TPR, mcp_FDR, mcp_TPR_cv, mcp_FDR_cv,
     iht_TPR, iht_FDR, iht_TPR_cv, iht_FDR_cv, file = paste0("./11c/result",seed,".RData"))




set.seed(seed)


## 11d
# covariance for x
rho = 0
S_x = diag(p)
for (i1 in c(1:p)){
  for (i2 in c(1:p)){
    # choice 1: exp decay
    S_x[i1,i2] = rho ^ (abs(i1 - i2))
    # choice 2: constant cor
    #if (i1 != i2) S_x[i1,i2] = rho     
  }
}

# signal and noise
betamin = 0.1
sigma_e = 1

# decide lambda values in picasso (scad)
lambda_list_scad <- 1.5/(1.2 ^ seq(0, 99))

# decide lambda values for glmnet (lasso)
lambda_list_lasso <- 1.5/(1.2 ^ seq(0, 99))

# decide lambda values for picasso (mcp)
lambda_list_mcp <- 1.5/(1.2 ^ seq(0, 99))

# decide pi values for iht
lambda_list_iht = seq(2,100,2)

## generate data
data_i <- data_lm(n = n, mu_x = rep(0,p), S_x = S_x,
                  s = s, betamin = betamin, sigma_e = sigma_e)
X = data_i$X
Y = data_i$Y  
beta = data_i$beta
epsilon = data_i$epsilon
supp_true = data_i$supp_true

## optimization with nonconcave penalty
scad_result = scad_select(X = X, Y = Y, supp_true = supp_true, lambda_list = lambda_list_scad)
scad_TPR = scad_result$TPR
scad_FDR = scad_result$FDR
scad_TPR_cv = scad_result$TPR_cv
scad_FDR_cv= scad_result$FDR_cv

## optimization with convex penalty
lasso_result = lasso_select(X = X, Y = Y, supp_true = supp_true, lambda_list = lambda_list_lasso)
lasso_TPR = lasso_result$TPR
lasso_FDR = lasso_result$FDR
lasso_TPR_cv = lasso_result$TPR_cv
lasso_FDR_cv = lasso_result$FDR_cv

## optimization with minimax concave penalty
mcp_result = mcp_select(X = X, Y = Y, supp_true = supp_true, lambda_list = lambda_list_mcp)
mcp_TPR = mcp_result$TPR
mcp_FDR = mcp_result$FDR
mcp_TPR_cv = mcp_result$TPR_cv
mcp_FDR_cv = mcp_result$FDR_cv

## optimization with iterative hard thresholding
iht_result = iht_m_select(X = X, Y = Y, supp_true = supp_true, pi_list = lambda_list_iht)
iht_TPR= iht_result$TPR
iht_FDR = iht_result$FDR
iht_TPR_cv = iht_result$TPR_cv
iht_FDR_cv = iht_result$FDR_cv

## save rdata
save(X, Y, beta, epsilon, supp_true, 
     scad_TPR, scad_FDR, scad_TPR_cv, scad_FDR_cv, 
     lasso_TPR, lasso_FDR, lasso_TPR_cv, lasso_FDR_cv,
     mcp_TPR, mcp_FDR, mcp_TPR_cv, mcp_FDR_cv,
     iht_TPR, iht_FDR, iht_TPR_cv, iht_FDR_cv, file = paste0("./11d/result",seed,".RData"))




set.seed(seed)


## 11e
# covariance for x
rho = 0.5
S_x = diag(p)
for (i1 in c(1:p)){
  for (i2 in c(1:p)){
    # choice 1: exp decay
    S_x[i1,i2] = rho ^ (abs(i1 - i2))
    # choice 2: constant cor
    #if (i1 != i2) S_x[i1,i2] = rho     
  }
}

# signal and noise
betamin = 0.1
sigma_e = 1

# decide lambda values in picasso (scad)
lambda_list_scad <- 1.5/(1.2 ^ seq(0, 99))

# decide lambda values for glmnet (lasso)
lambda_list_lasso <- 1.5/(1.2 ^ seq(0, 99))

# decide lambda values for picasso (mcp)
lambda_list_mcp <- 1.5/(1.2 ^ seq(0, 99))

# decide pi values for iht
lambda_list_iht = seq(2,100,2)

## generate data
data_i <- data_lm(n = n, mu_x = rep(0,p), S_x = S_x,
                  s = s, betamin = betamin, sigma_e = sigma_e)
X = data_i$X
Y = data_i$Y  
beta = data_i$beta
epsilon = data_i$epsilon
supp_true = data_i$supp_true

## optimization with nonconcave penalty
scad_result = scad_select(X = X, Y = Y, supp_true = supp_true, lambda_list = lambda_list_scad)
scad_TPR = scad_result$TPR
scad_FDR = scad_result$FDR
scad_TPR_cv = scad_result$TPR_cv
scad_FDR_cv= scad_result$FDR_cv

## optimization with convex penalty
lasso_result = lasso_select(X = X, Y = Y, supp_true = supp_true, lambda_list = lambda_list_lasso)
lasso_TPR = lasso_result$TPR
lasso_FDR = lasso_result$FDR
lasso_TPR_cv = lasso_result$TPR_cv
lasso_FDR_cv = lasso_result$FDR_cv

## optimization with minimax concave penalty
mcp_result = mcp_select(X = X, Y = Y, supp_true = supp_true, lambda_list = lambda_list_mcp)
mcp_TPR = mcp_result$TPR
mcp_FDR = mcp_result$FDR
mcp_TPR_cv = mcp_result$TPR_cv
mcp_FDR_cv = mcp_result$FDR_cv

## optimization with iterative hard thresholding
iht_result = iht_m_select(X = X, Y = Y, supp_true = supp_true, pi_list = lambda_list_iht)
iht_TPR= iht_result$TPR
iht_FDR = iht_result$FDR
iht_TPR_cv = iht_result$TPR_cv
iht_FDR_cv = iht_result$FDR_cv

## save rdata
save(X, Y, beta, epsilon, supp_true, 
     scad_TPR, scad_FDR, scad_TPR_cv, scad_FDR_cv, 
     lasso_TPR, lasso_FDR, lasso_TPR_cv, lasso_FDR_cv,
     mcp_TPR, mcp_FDR, mcp_TPR_cv, mcp_FDR_cv,
     iht_TPR, iht_FDR, iht_TPR_cv, iht_FDR_cv, file = paste0("./11e/result",seed,".RData"))




set.seed(seed)


## 11f
# covariance for x
rho = 0.8
S_x = diag(p)
for (i1 in c(1:p)){
  for (i2 in c(1:p)){
    # choice 1: exp decay
    S_x[i1,i2] = rho ^ (abs(i1 - i2))
    # choice 2: constant cor
    #if (i1 != i2) S_x[i1,i2] = rho     
  }
}

# signal and noise
betamin = 0.1
sigma_e = 1

# decide lambda values in picasso (scad)
lambda_list_scad <- 1.5/(1.2 ^ seq(0, 99))

# decide lambda values for glmnet (lasso)
lambda_list_lasso <- 1.5/(1.2 ^ seq(0, 99))

# decide lambda values for picasso (mcp)
lambda_list_mcp <- 1.5/(1.2 ^ seq(0, 99))

# decide pi values for iht
lambda_list_iht = seq(2,100,2)

## generate data
data_i <- data_lm(n = n, mu_x = rep(0,p), S_x = S_x,
                  s = s, betamin = betamin, sigma_e = sigma_e)
X = data_i$X
Y = data_i$Y  
beta = data_i$beta
epsilon = data_i$epsilon
supp_true = data_i$supp_true

## optimization with nonconcave penalty
scad_result = scad_select(X = X, Y = Y, supp_true = supp_true, lambda_list = lambda_list_scad)
scad_TPR = scad_result$TPR
scad_FDR = scad_result$FDR
scad_TPR_cv = scad_result$TPR_cv
scad_FDR_cv= scad_result$FDR_cv

## optimization with convex penalty
lasso_result = lasso_select(X = X, Y = Y, supp_true = supp_true, lambda_list = lambda_list_lasso)
lasso_TPR = lasso_result$TPR
lasso_FDR = lasso_result$FDR
lasso_TPR_cv = lasso_result$TPR_cv
lasso_FDR_cv = lasso_result$FDR_cv

## optimization with minimax concave penalty
mcp_result = mcp_select(X = X, Y = Y, supp_true = supp_true, lambda_list = lambda_list_mcp)
mcp_TPR = mcp_result$TPR
mcp_FDR = mcp_result$FDR
mcp_TPR_cv = mcp_result$TPR_cv
mcp_FDR_cv = mcp_result$FDR_cv

## optimization with iterative hard thresholding
iht_result = iht_m_select(X = X, Y = Y, supp_true = supp_true, pi_list = lambda_list_iht)
iht_TPR= iht_result$TPR
iht_FDR = iht_result$FDR
iht_TPR_cv = iht_result$TPR_cv
iht_FDR_cv = iht_result$FDR_cv

## save rdata
save(X, Y, beta, epsilon, supp_true, 
     scad_TPR, scad_FDR, scad_TPR_cv, scad_FDR_cv, 
     lasso_TPR, lasso_FDR, lasso_TPR_cv, lasso_FDR_cv,
     mcp_TPR, mcp_FDR, mcp_TPR_cv, mcp_FDR_cv,
     iht_TPR, iht_FDR, iht_TPR_cv, iht_FDR_cv, file = paste0("./11f/result",seed,".RData"))



set.seed(seed)



## 12a
# covariance for x
rho = 0.5
S_x = matrix(rho, nrow = p, ncol = p) + diag(1-rho, nrow =p, ncol = p)

# signal and noise
betamin = 0.1
sigma_e = 0.5

# decide lambda values in picasso (scad)
lambda_list_scad <- 1.5/(1.2 ^ seq(0, 99))

# decide lambda values for glmnet (lasso)
lambda_list_lasso <- 1.5/(1.2 ^ seq(0, 99))

# decide lambda values for picasso (mcp)
lambda_list_mcp <- 1.5/(1.2 ^ seq(0, 99))

# decide pi values for iht
lambda_list_iht = seq(2,100,2)

## generate data
data_i <- data_lm(n = n, mu_x = rep(0,p), S_x = S_x,
                  s = s, betamin = betamin, sigma_e = sigma_e)
X = data_i$X
Y = data_i$Y  
beta = data_i$beta
epsilon = data_i$epsilon
supp_true = data_i$supp_true

## optimization with nonconcave penalty
scad_result = scad_select(X = X, Y = Y, supp_true = supp_true, lambda_list = lambda_list_scad)
scad_TPR = scad_result$TPR
scad_FDR = scad_result$FDR
scad_TPR_cv = scad_result$TPR_cv
scad_FDR_cv= scad_result$FDR_cv

## optimization with convex penalty
lasso_result = lasso_select(X = X, Y = Y, supp_true = supp_true, lambda_list = lambda_list_lasso)
lasso_TPR = lasso_result$TPR
lasso_FDR = lasso_result$FDR
lasso_TPR_cv = lasso_result$TPR_cv
lasso_FDR_cv = lasso_result$FDR_cv

## optimization with minimax concave penalty
mcp_result = mcp_select(X = X, Y = Y, supp_true = supp_true, lambda_list = lambda_list_mcp)
mcp_TPR = mcp_result$TPR
mcp_FDR = mcp_result$FDR
mcp_TPR_cv = mcp_result$TPR_cv
mcp_FDR_cv = mcp_result$FDR_cv

## optimization with iterative hard thresholding
iht_result = iht_m_select(X = X, Y = Y, supp_true = supp_true, pi_list = lambda_list_iht)
iht_TPR= iht_result$TPR
iht_FDR = iht_result$FDR
iht_TPR_cv = iht_result$TPR_cv
iht_FDR_cv = iht_result$FDR_cv

## save rdata
save(X, Y, beta, epsilon, supp_true, 
     scad_TPR, scad_FDR, scad_TPR_cv, scad_FDR_cv, 
     lasso_TPR, lasso_FDR, lasso_TPR_cv, lasso_FDR_cv,
     mcp_TPR, mcp_FDR, mcp_TPR_cv, mcp_FDR_cv,
     iht_TPR, iht_FDR, iht_TPR_cv, iht_FDR_cv, file = paste0("./12a/result",seed,".RData"))



set.seed(seed)

## 12b
# covariance for x
rho = 0.8
S_x = matrix(rho, nrow = p, ncol = p) + diag(1-rho, nrow =p, ncol = p)

# signal and noise
betamin = 0.1
sigma_e = 0.5

# decide lambda values in picasso (scad)
lambda_list_scad <- 1.5/(1.2 ^ seq(0, 99))

# decide lambda values for glmnet (lasso)
lambda_list_lasso <- 1.5/(1.2 ^ seq(0, 99))

# decide lambda values for picasso (mcp)
lambda_list_mcp <- 1.5/(1.2 ^ seq(0, 99))

# decide pi values for iht
lambda_list_iht = seq(2,100,2)

## generate data
data_i <- data_lm(n = n, mu_x = rep(0,p), S_x = S_x,
                  s = s, betamin = betamin, sigma_e = sigma_e)
X = data_i$X
Y = data_i$Y  
beta = data_i$beta
epsilon = data_i$epsilon
supp_true = data_i$supp_true

## optimization with nonconcave penalty
scad_result = scad_select(X = X, Y = Y, supp_true = supp_true, lambda_list = lambda_list_scad)
scad_TPR = scad_result$TPR
scad_FDR = scad_result$FDR
scad_TPR_cv = scad_result$TPR_cv
scad_FDR_cv= scad_result$FDR_cv

## optimization with convex penalty
lasso_result = lasso_select(X = X, Y = Y, supp_true = supp_true, lambda_list = lambda_list_lasso)
lasso_TPR = lasso_result$TPR
lasso_FDR = lasso_result$FDR
lasso_TPR_cv = lasso_result$TPR_cv
lasso_FDR_cv = lasso_result$FDR_cv

## optimization with minimax concave penalty
mcp_result = mcp_select(X = X, Y = Y, supp_true = supp_true, lambda_list = lambda_list_mcp)
mcp_TPR = mcp_result$TPR
mcp_FDR = mcp_result$FDR
mcp_TPR_cv = mcp_result$TPR_cv
mcp_FDR_cv = mcp_result$FDR_cv

## optimization with iterative hard thresholding
iht_result = iht_m_select(X = X, Y = Y, supp_true = supp_true, pi_list = lambda_list_iht)
iht_TPR= iht_result$TPR
iht_FDR = iht_result$FDR
iht_TPR_cv = iht_result$TPR_cv
iht_FDR_cv = iht_result$FDR_cv

## save rdata
save(X, Y, beta, epsilon, supp_true, 
     scad_TPR, scad_FDR, scad_TPR_cv, scad_FDR_cv, 
     lasso_TPR, lasso_FDR, lasso_TPR_cv, lasso_FDR_cv,
     mcp_TPR, mcp_FDR, mcp_TPR_cv, mcp_FDR_cv,
     iht_TPR, iht_FDR, iht_TPR_cv, iht_FDR_cv, file = paste0("./12b/result",seed,".RData"))



set.seed(seed)



## 12c
# covariance for x
rho = 0.5
S_x = matrix(rho, nrow = p, ncol = p) + diag(1-rho, nrow =p, ncol = p)

# signal and noise
betamin = 0.1
sigma_e = 1

# decide lambda values in picasso (scad)
lambda_list_scad <- 1.5/(1.2 ^ seq(0, 99))

# decide lambda values for glmnet (lasso)
lambda_list_lasso <- 1.5/(1.2 ^ seq(0, 99))

# decide lambda values for picasso (mcp)
lambda_list_mcp <- 1.5/(1.2 ^ seq(0, 99))


# decide pi values for iht
lambda_list_iht = seq(2,100,2)

## generate data
data_i <- data_lm(n = n, mu_x = rep(0,p), S_x = S_x,
                  s = s, betamin = betamin, sigma_e = sigma_e)
X = data_i$X
Y = data_i$Y  
beta = data_i$beta
epsilon = data_i$epsilon
supp_true = data_i$supp_true

## optimization with nonconcave penalty
scad_result = scad_select(X = X, Y = Y, supp_true = supp_true, lambda_list = lambda_list_scad)
scad_TPR = scad_result$TPR
scad_FDR = scad_result$FDR
scad_TPR_cv = scad_result$TPR_cv
scad_FDR_cv= scad_result$FDR_cv

## optimization with convex penalty
lasso_result = lasso_select(X = X, Y = Y, supp_true = supp_true, lambda_list = lambda_list_lasso)
lasso_TPR = lasso_result$TPR
lasso_FDR = lasso_result$FDR
lasso_TPR_cv = lasso_result$TPR_cv
lasso_FDR_cv = lasso_result$FDR_cv

## optimization with minimax concave penalty
mcp_result = mcp_select(X = X, Y = Y, supp_true = supp_true, lambda_list = lambda_list_mcp)
mcp_TPR = mcp_result$TPR
mcp_FDR = mcp_result$FDR
mcp_TPR_cv = mcp_result$TPR_cv
mcp_FDR_cv = mcp_result$FDR_cv

## optimization with iterative hard thresholding
iht_result = iht_m_select(X = X, Y = Y, supp_true = supp_true, pi_list = lambda_list_iht)
iht_TPR= iht_result$TPR
iht_FDR = iht_result$FDR
iht_TPR_cv = iht_result$TPR_cv
iht_FDR_cv = iht_result$FDR_cv

## save rdata
save(X, Y, beta, epsilon, supp_true, 
     scad_TPR, scad_FDR, scad_TPR_cv, scad_FDR_cv, 
     lasso_TPR, lasso_FDR, lasso_TPR_cv, lasso_FDR_cv,
     mcp_TPR, mcp_FDR, mcp_TPR_cv, mcp_FDR_cv,
     iht_TPR, iht_FDR, iht_TPR_cv, iht_FDR_cv, file = paste0("./12c/result",seed,".RData"))



set.seed(seed)

## 12d
# covariance for x
rho = 0.8
S_x = matrix(rho, nrow = p, ncol = p) + diag(1-rho, nrow =p, ncol = p)

# signal and noise
betamin = 0.1
sigma_e = 1

# decide lambda values in picasso (scad)
lambda_list_scad <- 1.5/(1.2 ^ seq(0, 99))

# decide lambda values for glmnet (lasso)
lambda_list_lasso <- 1.5/(1.2 ^ seq(0, 99))

# decide lambda values for picasso (mcp)
lambda_list_mcp <- 1.5/(1.2 ^ seq(0, 99))

# decide pi values for iht
lambda_list_iht = seq(2,100,2)

## generate data
data_i <- data_lm(n = n, mu_x = rep(0,p), S_x = S_x,
                  s = s, betamin = betamin, sigma_e = sigma_e)
X = data_i$X
Y = data_i$Y  
beta = data_i$beta
epsilon = data_i$epsilon
supp_true = data_i$supp_true

## optimization with nonconcave penalty
scad_result = scad_select(X = X, Y = Y, supp_true = supp_true, lambda_list = lambda_list_scad)
scad_TPR = scad_result$TPR
scad_FDR = scad_result$FDR
scad_TPR_cv = scad_result$TPR_cv
scad_FDR_cv= scad_result$FDR_cv

## optimization with convex penalty
lasso_result = lasso_select(X = X, Y = Y, supp_true = supp_true, lambda_list = lambda_list_lasso)
lasso_TPR = lasso_result$TPR
lasso_FDR = lasso_result$FDR
lasso_TPR_cv = lasso_result$TPR_cv
lasso_FDR_cv = lasso_result$FDR_cv

## optimization with minimax concave penalty
mcp_result = mcp_select(X = X, Y = Y, supp_true = supp_true, lambda_list = lambda_list_mcp)
mcp_TPR = mcp_result$TPR
mcp_FDR = mcp_result$FDR
mcp_TPR_cv = mcp_result$TPR_cv
mcp_FDR_cv = mcp_result$FDR_cv

## optimization with iterative hard thresholding
iht_result = iht_m_select(X = X, Y = Y, supp_true = supp_true, pi_list = lambda_list_iht)
iht_TPR= iht_result$TPR
iht_FDR = iht_result$FDR
iht_TPR_cv = iht_result$TPR_cv
iht_FDR_cv = iht_result$FDR_cv

## save rdata
save(X, Y, beta, epsilon, supp_true, 
     scad_TPR, scad_FDR, scad_TPR_cv, scad_FDR_cv, 
     lasso_TPR, lasso_FDR, lasso_TPR_cv, lasso_FDR_cv,
     mcp_TPR, mcp_FDR, mcp_TPR_cv, mcp_FDR_cv,
     iht_TPR, iht_FDR, iht_TPR_cv, iht_FDR_cv, file = paste0("./12d/result",seed,".RData"))








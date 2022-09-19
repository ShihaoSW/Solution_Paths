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



## 11a
# covariance for x
S_x = diag(p)



# signal and noise
betamin = 0.1
sigma_e = 0.5

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

## optimization with iterative hard thresholding
iht1_result = iht1_select(X = X, Y = Y, supp_true = supp_true, pi_list = lambda_list_iht)

iht1_TPR= iht1_result$TPR
iht1_FDR = iht1_result$FDR
iht1_TPR_cv = iht1_result$TPR_cv
iht1_FDR_cv = iht1_result$FDR_cv



## save rdata
save(X, Y, beta, epsilon, supp_true, 
     iht1_TPR, iht1_FDR, iht1_TPR_cv, iht1_FDR_cv, file = paste0("./11a/result_iht",seed,".RData"))

###### release memory
rm(list = ls())


############################################################
############################################################ again
############################################################

###### high dimensional model selection

library(MASS)
library(Matrix)
library(picasso)
library(glmnet)
source("IHT_real.R")

###### set seed
seed <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(seed)


####### unchanged parameters
p = 5000
s = 50
n = ceiling(2 * s * log(p))



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

## optimization with iterative hard thresholding
iht1_result = iht1_select(X = X, Y = Y, supp_true = supp_true, pi_list = lambda_list_iht)

iht1_TPR= iht1_result$TPR
iht1_FDR = iht1_result$FDR
iht1_TPR_cv = iht1_result$TPR_cv
iht1_FDR_cv = iht1_result$FDR_cv



## save rdata
save(X, Y, beta, epsilon, supp_true, 
     iht1_TPR, iht1_FDR, iht1_TPR_cv, iht1_FDR_cv, file = paste0("./11b/result_iht",seed,".RData"))

###### release memory
rm(list = ls())





############################################################
############################################################ again
############################################################

###### high dimensional model selection

library(MASS)
library(Matrix)
library(picasso)
library(glmnet)
source("IHT_real.R")

###### set seed
seed <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(seed)


####### unchanged parameters
p = 5000
s = 50
n = ceiling(2 * s * log(p))



## 11d
# covariance for x
S_x = diag(p)


# signal and noise
betamin = 0.1
sigma_e = 1

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

## optimization with iterative hard thresholding
iht1_result = iht1_select(X = X, Y = Y, supp_true = supp_true, pi_list = lambda_list_iht)

iht1_TPR= iht1_result$TPR
iht1_FDR = iht1_result$FDR
iht1_TPR_cv = iht1_result$TPR_cv
iht1_FDR_cv = iht1_result$FDR_cv



## save rdata
save(X, Y, beta, epsilon, supp_true, 
     iht1_TPR, iht1_FDR, iht1_TPR_cv, iht1_FDR_cv, file = paste0("./11d/result_iht",seed,".RData"))

###### release memory
rm(list = ls())




############################################################
############################################################ again
############################################################

###### high dimensional model selection

library(MASS)
library(Matrix)
library(picasso)
library(glmnet)
source("IHT_real.R")

###### set seed
seed <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(seed)


####### unchanged parameters
p = 5000
s = 50
n = ceiling(2 * s * log(p))



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

## optimization with iterative hard thresholding
iht1_result = iht1_select(X = X, Y = Y, supp_true = supp_true, pi_list = lambda_list_iht)

iht1_TPR= iht1_result$TPR
iht1_FDR = iht1_result$FDR
iht1_TPR_cv = iht1_result$TPR_cv
iht1_FDR_cv = iht1_result$FDR_cv



## save rdata
save(X, Y, beta, epsilon, supp_true, 
     iht1_TPR, iht1_FDR, iht1_TPR_cv, iht1_FDR_cv, file = paste0("./11e/result_iht",seed,".RData"))

###### release memory
rm(list = ls())






############################################################
############################################################ again
############################################################

###### high dimensional model selection

library(MASS)
library(Matrix)
library(picasso)
library(glmnet)
source("IHT_real.R")

###### set seed
seed <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(seed)


####### unchanged parameters
p = 5000
s = 50
n = ceiling(2 * s * log(p))



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

## optimization with iterative hard thresholding
iht1_result = iht1_select(X = X, Y = Y, supp_true = supp_true, pi_list = lambda_list_iht)

iht1_TPR= iht1_result$TPR
iht1_FDR = iht1_result$FDR
iht1_TPR_cv = iht1_result$TPR_cv
iht1_FDR_cv = iht1_result$FDR_cv



## save rdata
save(X, Y, beta, epsilon, supp_true, 
     iht1_TPR, iht1_FDR, iht1_TPR_cv, iht1_FDR_cv, file = paste0("./11f/result_iht",seed,".RData"))

###### release memory
rm(list = ls())




############################################################
############################################################ again
############################################################

###### high dimensional model selection

library(MASS)
library(Matrix)
library(picasso)
library(glmnet)
source("IHT_real.R")

###### set seed
seed <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(seed)


####### unchanged parameters
p = 5000
s = 50
n = ceiling(2 * s * log(p))



## 12a
# covariance for x
rho = 0.5
S_x = matrix(rho, nrow = p, ncol = p) + diag(1-rho, nrow =p, ncol = p)

# signal and noise
betamin = 0.1
sigma_e = 0.5

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

## optimization with iterative hard thresholding
iht1_result = iht1_select(X = X, Y = Y, supp_true = supp_true, pi_list = lambda_list_iht)

iht1_TPR= iht1_result$TPR
iht1_FDR = iht1_result$FDR
iht1_TPR_cv = iht1_result$TPR_cv
iht1_FDR_cv = iht1_result$FDR_cv



## save rdata
save(X, Y, beta, epsilon, supp_true, 
     iht1_TPR, iht1_FDR, iht1_TPR_cv, iht1_FDR_cv, file = paste0("./12a/result_iht",seed,".RData"))

###### release memory
rm(list = ls())


############################################################
############################################################ again
############################################################

###### high dimensional model selection

library(MASS)
library(Matrix)
library(picasso)
library(glmnet)
source("IHT_real.R")

###### set seed
seed <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(seed)


####### unchanged parameters
p = 5000
s = 50
n = ceiling(2 * s * log(p))



## 12b
# covariance for x
rho = 0.8
S_x = matrix(rho, nrow = p, ncol = p) + diag(1-rho, nrow =p, ncol = p)

# signal and noise
betamin = 0.1
sigma_e = 0.5

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

## optimization with iterative hard thresholding
iht1_result = iht1_select(X = X, Y = Y, supp_true = supp_true, pi_list = lambda_list_iht)

iht1_TPR= iht1_result$TPR
iht1_FDR = iht1_result$FDR
iht1_TPR_cv = iht1_result$TPR_cv
iht1_FDR_cv = iht1_result$FDR_cv



## save rdata
save(X, Y, beta, epsilon, supp_true, 
     iht1_TPR, iht1_FDR, iht1_TPR_cv, iht1_FDR_cv, file = paste0("./12b/result_iht",seed,".RData"))

###### release memory
rm(list = ls())




############################################################
############################################################ again
############################################################

###### high dimensional model selection

library(MASS)
library(Matrix)
library(picasso)
library(glmnet)
source("IHT_real.R")

###### set seed
seed <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(seed)


####### unchanged parameters
p = 5000
s = 50
n = ceiling(2 * s * log(p))



## 12c
# covariance for x
rho = 0.5
S_x = matrix(rho, nrow = p, ncol = p) + diag(1-rho, nrow =p, ncol = p)

# signal and noise
betamin = 0.1
sigma_e = 1

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

## optimization with iterative hard thresholding
iht1_result = iht1_select(X = X, Y = Y, supp_true = supp_true, pi_list = lambda_list_iht)

iht1_TPR= iht1_result$TPR
iht1_FDR = iht1_result$FDR
iht1_TPR_cv = iht1_result$TPR_cv
iht1_FDR_cv = iht1_result$FDR_cv



## save rdata
save(X, Y, beta, epsilon, supp_true, 
     iht1_TPR, iht1_FDR, iht1_TPR_cv, iht1_FDR_cv, file = paste0("./12c/result_iht",seed,".RData"))

###### release memory
rm(list = ls())



############################################################
############################################################ again
############################################################

###### high dimensional model selection

library(MASS)
library(Matrix)
library(picasso)
library(glmnet)
source("IHT_real.R")

###### set seed
seed <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(seed)


####### unchanged parameters
p = 5000
s = 50
n = ceiling(2 * s * log(p))



## 12d
# covariance for x
rho = 0.8
S_x = matrix(rho, nrow = p, ncol = p) + diag(1-rho, nrow =p, ncol = p)

# signal and noise
betamin = 0.1
sigma_e = 1

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

## optimization with iterative hard thresholding
iht1_result = iht1_select(X = X, Y = Y, supp_true = supp_true, pi_list = lambda_list_iht)

iht1_TPR= iht1_result$TPR
iht1_FDR = iht1_result$FDR
iht1_TPR_cv = iht1_result$TPR_cv
iht1_FDR_cv = iht1_result$FDR_cv



## save rdata
save(X, Y, beta, epsilon, supp_true, 
     iht1_TPR, iht1_FDR, iht1_TPR_cv, iht1_FDR_cv, file = paste0("./12d/result_iht",seed,".RData"))

###### release memory
rm(list = ls())

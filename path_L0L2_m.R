###### high dimensional model selection by L0Learn

library(MASS)
library(Matrix)
library(L0Learn)
source("functions.R")

###### set seed ########
seed <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(seed)
 

####### unchanged parameters. 
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

# decide lambda values and gamma in L0Learn
nlambda = 100
lambda_max = 0.1
lambda_ratio = 1.1
ratios = lambda_ratio ^ seq(0, nlambda-1)
lambda_list = list(lambda_max / ratios)

gamma_input = 1e-05


## generate data

data_i <- data_lm(n = n, mu_x = rep(0,p), S_x = S_x,
                  s = s, betamin = betamin, sigma_e = sigma_e)

X = data_i$X
Y = data_i$Y  
beta = data_i$beta
epsilon = data_i$epsilon
supp_true = data_i$supp_true

## L0Learn
result = L0Learn_path2(X, Y, supp_true, method = 'L0L2', max_size = 500,      
                       lambda_list = lambda_list, gamma_input = gamma_input)

## optimization with iterative hard thresholding
l0l2_FDR = result$FDR
l0l2_TPR = result$TPR

## save rdata
save(l0l2_FDR, l0l2_TPR, file = paste0("./11a/l0l2",seed,".RData"))




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

# decide lambda values and gamma in L0Learn
nlambda = 100
lambda_max = 0.1
lambda_ratio = 1.1
ratios = lambda_ratio ^ seq(0, nlambda-1)
lambda_list = list(lambda_max / ratios)

gamma_input = 1e-05


## generate data

data_i <- data_lm(n = n, mu_x = rep(0,p), S_x = S_x,
                  s = s, betamin = betamin, sigma_e = sigma_e)

X = data_i$X
Y = data_i$Y  
beta = data_i$beta
epsilon = data_i$epsilon
supp_true = data_i$supp_true

## L0Learn
result = L0Learn_path2(X, Y, supp_true, method = 'L0L2', max_size = 500,      
                       lambda_list = lambda_list, gamma_input = gamma_input)

## optimization with iterative hard thresholding
l0l2_FDR = result$FDR
l0l2_TPR = result$TPR

## save rdata
save(l0l2_FDR, l0l2_TPR, file = paste0("./11b/l0l2",seed,".RData"))


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


# decide lambda values and gamma in L0Learn
nlambda = 100
lambda_max = 0.1
lambda_ratio = 1.1
ratios = lambda_ratio ^ seq(0, nlambda-1)
lambda_list = list(lambda_max / ratios)

gamma_input = 1e-05


## generate data

data_i <- data_lm(n = n, mu_x = rep(0,p), S_x = S_x,
                  s = s, betamin = betamin, sigma_e = sigma_e)

X = data_i$X
Y = data_i$Y  
beta = data_i$beta
epsilon = data_i$epsilon
supp_true = data_i$supp_true

## L0Learn
result = L0Learn_path2(X, Y, supp_true, method = 'L0L2', max_size = 500,      
                       lambda_list = lambda_list, gamma_input = gamma_input)

## optimization with iterative hard thresholding
l0l2_FDR = result$FDR
l0l2_TPR = result$TPR

## save rdata
save(l0l2_FDR, l0l2_TPR, file = paste0("./11c/l0l2",seed,".RData"))


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


# decide lambda values and gamma in L0Learn
nlambda = 100
lambda_max = 0.1
lambda_ratio = 1.1
ratios = lambda_ratio ^ seq(0, nlambda-1)
lambda_list = list(lambda_max / ratios)

gamma_input = 1e-05


## generate data

data_i <- data_lm(n = n, mu_x = rep(0,p), S_x = S_x,
                  s = s, betamin = betamin, sigma_e = sigma_e)

X = data_i$X
Y = data_i$Y  
beta = data_i$beta
epsilon = data_i$epsilon
supp_true = data_i$supp_true

## L0Learn
result = L0Learn_path2(X, Y, supp_true, method = 'L0L2', max_size = 500,      
                       lambda_list = lambda_list, gamma_input = gamma_input)

## optimization with iterative hard thresholding
l0l2_FDR = result$FDR
l0l2_TPR = result$TPR

## save rdata
save(l0l2_FDR, l0l2_TPR, file = paste0("./11d/l0l2",seed,".RData"))




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


# decide lambda values and gamma in L0Learn
nlambda = 100
lambda_max = 0.1
lambda_ratio = 1.1
ratios = lambda_ratio ^ seq(0, nlambda-1)
lambda_list = list(lambda_max / ratios)

gamma_input = 1e-05


## generate data

data_i <- data_lm(n = n, mu_x = rep(0,p), S_x = S_x,
                  s = s, betamin = betamin, sigma_e = sigma_e)

X = data_i$X
Y = data_i$Y  
beta = data_i$beta
epsilon = data_i$epsilon
supp_true = data_i$supp_true

## L0Learn
result = L0Learn_path2(X, Y, supp_true, method = 'L0L2', max_size = 500,      
                       lambda_list = lambda_list, gamma_input = gamma_input)

## optimization with iterative hard thresholding
l0l2_FDR = result$FDR
l0l2_TPR = result$TPR

## save rdata
save(l0l2_FDR, l0l2_TPR, file = paste0("./11e/l0l2",seed,".RData"))





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


# decide lambda values and gamma in L0Learn
nlambda = 100
lambda_max = 0.1
lambda_ratio = 1.1
ratios = lambda_ratio ^ seq(0, nlambda-1)
lambda_list = list(lambda_max / ratios)

gamma_input = 1e-05


## generate data

data_i <- data_lm(n = n, mu_x = rep(0,p), S_x = S_x,
                  s = s, betamin = betamin, sigma_e = sigma_e)

X = data_i$X
Y = data_i$Y  
beta = data_i$beta
epsilon = data_i$epsilon
supp_true = data_i$supp_true

## L0Learn
result = L0Learn_path2(X, Y, supp_true, method = 'L0L2', max_size = 500,      
                       lambda_list = lambda_list, gamma_input = gamma_input)

## optimization with iterative hard thresholding
l0l2_FDR = result$FDR
l0l2_TPR = result$TPR

## save rdata
save(l0l2_FDR, l0l2_TPR, file = paste0("./11f/l0l2",seed,".RData"))

set.seed(seed)



## 12a
# covariance for x
rho = 0.5
S_x = matrix(rho, nrow = p, ncol = p) + diag(1-rho, nrow =p, ncol = p)

# signal and noise
betamin = 0.1
sigma_e = 0.5

# decide lambda values and gamma in L0Learn
nlambda = 100
lambda_max = 0.1
lambda_ratio = 1.1
ratios = lambda_ratio ^ seq(0, nlambda-1)
lambda_list = list(lambda_max / ratios)

gamma_input = 1e-05


## generate data

data_i <- data_lm(n = n, mu_x = rep(0,p), S_x = S_x,
                  s = s, betamin = betamin, sigma_e = sigma_e)

X = data_i$X
Y = data_i$Y  
beta = data_i$beta
epsilon = data_i$epsilon
supp_true = data_i$supp_true

## L0Learn
result = L0Learn_path2(X, Y, supp_true, method = 'L0L2', max_size = 500,      
                       lambda_list = lambda_list, gamma_input = gamma_input)

## optimization with iterative hard thresholding
l0l2_FDR = result$FDR
l0l2_TPR = result$TPR

## save rdata
save(l0l2_FDR, l0l2_TPR, file = paste0("./12a/l0l2",seed,".RData"))



set.seed(seed)

## 12b
# covariance for x
rho = 0.8
S_x = matrix(rho, nrow = p, ncol = p) + diag(1-rho, nrow =p, ncol = p)

# signal and noise
betamin = 0.1
sigma_e = 0.5

# decide lambda values and gamma in L0Learn
nlambda = 100
lambda_max = 0.1
lambda_ratio = 1.1
ratios = lambda_ratio ^ seq(0, nlambda-1)
lambda_list = list(lambda_max / ratios)

gamma_input = 1e-05


## generate data

data_i <- data_lm(n = n, mu_x = rep(0,p), S_x = S_x,
                  s = s, betamin = betamin, sigma_e = sigma_e)

X = data_i$X
Y = data_i$Y  
beta = data_i$beta
epsilon = data_i$epsilon
supp_true = data_i$supp_true

## L0Learn
result = L0Learn_path2(X, Y, supp_true, method = 'L0L2', max_size = 500,      
                       lambda_list = lambda_list, gamma_input = gamma_input)

## optimization with iterative hard thresholding
l0l2_FDR = result$FDR
l0l2_TPR = result$TPR

## save rdata
save(l0l2_FDR, l0l2_TPR, file = paste0("./12b/l0l2",seed,".RData"))




set.seed(seed)



## 12c
# covariance for x
rho = 0.5
S_x = matrix(rho, nrow = p, ncol = p) + diag(1-rho, nrow =p, ncol = p)

# signal and noise
betamin = 0.1
sigma_e = 1

# decide lambda values and gamma in L0Learn
nlambda = 100
lambda_max = 0.1
lambda_ratio = 1.1
ratios = lambda_ratio ^ seq(0, nlambda-1)
lambda_list = list(lambda_max / ratios)

gamma_input = 1e-05


## generate data

data_i <- data_lm(n = n, mu_x = rep(0,p), S_x = S_x,
                  s = s, betamin = betamin, sigma_e = sigma_e)

X = data_i$X
Y = data_i$Y  
beta = data_i$beta
epsilon = data_i$epsilon
supp_true = data_i$supp_true

## L0Learn
result = L0Learn_path2(X, Y, supp_true, method = 'L0L2', max_size = 500,      
                       lambda_list = lambda_list, gamma_input = gamma_input)

## optimization with iterative hard thresholding
l0l2_FDR = result$FDR
l0l2_TPR = result$TPR

## save rdata
save(l0l2_FDR, l0l2_TPR, file = paste0("./12c/l0l2",seed,".RData"))




set.seed(seed)

## 12d
# covariance for x
rho = 0.8
S_x = matrix(rho, nrow = p, ncol = p) + diag(1-rho, nrow =p, ncol = p)

# signal and noise
betamin = 0.1
sigma_e = 1

# decide lambda values and gamma in L0Learn
nlambda = 100
lambda_max = 0.1
lambda_ratio = 1.1
ratios = lambda_ratio ^ seq(0, nlambda-1)
lambda_list = list(lambda_max / ratios)

gamma_input = 1e-05


## generate data

data_i <- data_lm(n = n, mu_x = rep(0,p), S_x = S_x,
                  s = s, betamin = betamin, sigma_e = sigma_e)

X = data_i$X
Y = data_i$Y  
beta = data_i$beta
epsilon = data_i$epsilon
supp_true = data_i$supp_true

## L0Learn
result = L0Learn_path2(X, Y, supp_true, method = 'L0L2', max_size = 500,      
                       lambda_list = lambda_list, gamma_input = gamma_input)

## optimization with iterative hard thresholding
l0l2_FDR = result$FDR
l0l2_TPR = result$TPR

## save rdata
save(l0l2_FDR, l0l2_TPR, file = paste0("./12d/l0l2",seed,".RData"))




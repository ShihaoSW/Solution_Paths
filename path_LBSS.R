###### high dimensional model selection

library(MASS)
library(Matrix)
library(picasso)
library(glmnet)
library(gurobi)
library(bestsubset)
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

# decide s_list values for iht
s_list = c(c(1,2,3,4), seq(5,100,5))

## generate data

data_i <- data_lm(n = n, mu_x = rep(0,p), S_x = S_x,
                  s = s, betamin = betamin, sigma_e = sigma_e)

X = data_i$X
Y = data_i$Y  
beta = data_i$beta
epsilon = data_i$epsilon
supp_true = data_i$supp_true

## get the path
lassobss_result = lassobss(X, Y, supp_true = supp_true, s_list = s_list, scr_dim = 300, timelimit = 600)
lassobss_TPR= lassobss_result$TPR
lassobss_FDR = lassobss_result$FDR


save(beta, epsilon, supp_true, 
     lassobss_TPR, lassobss_FDR, file = paste0("./11a/result_lassobss",seed,".RData"))




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

# decide s_list values for iht
s_list = c(c(1,2,3,4), seq(5,100,5))

## generate data
data_i <- data_lm(n = n, mu_x = rep(0,p), S_x = S_x,
                  s = s, betamin = betamin, sigma_e = sigma_e)
X = data_i$X
Y = data_i$Y  
beta = data_i$beta
epsilon = data_i$epsilon
supp_true = data_i$supp_true

## get the path
lassobss_result = lassobss(X, Y, supp_true = supp_true, s_list = s_list, scr_dim = 300, timelimit = 600)
lassobss_TPR= lassobss_result$TPR
lassobss_FDR = lassobss_result$FDR


save(beta, epsilon, supp_true, 
     lassobss_TPR, lassobss_FDR, file = paste0("./11b/result_lassobss",seed,".RData"))












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


# decide s_list values for iht
s_list = c(c(1,2,3,4), seq(5,100,5))

## generate data

data_i <- data_lm(n = n, mu_x = rep(0,p), S_x = S_x,
                  s = s, betamin = betamin, sigma_e = sigma_e)

X = data_i$X
Y = data_i$Y  
beta = data_i$beta
epsilon = data_i$epsilon
supp_true = data_i$supp_true



## get the path
lassobss_result = lassobss(X, Y, supp_true = supp_true, s_list = s_list, scr_dim = 300, timelimit = 600)
lassobss_TPR= lassobss_result$TPR
lassobss_FDR = lassobss_result$FDR


save(beta, epsilon, supp_true, 
     lassobss_TPR, lassobss_FDR, file = paste0("./11c/result_lassobss",seed,".RData"))










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

# decide s_list values for iht
s_list = c(c(1,2,3,4), seq(5,100,5))

## generate data
data_i <- data_lm(n = n, mu_x = rep(0,p), S_x = S_x,
                  s = s, betamin = betamin, sigma_e = sigma_e)
X = data_i$X
Y = data_i$Y  
beta = data_i$beta
epsilon = data_i$epsilon
supp_true = data_i$supp_true

## get the path
lassobss_result = lassobss(X, Y, supp_true = supp_true, s_list = s_list, scr_dim = 300, timelimit = 600)
lassobss_TPR= lassobss_result$TPR
lassobss_FDR = lassobss_result$FDR


save(beta, epsilon, supp_true, 
     lassobss_TPR, lassobss_FDR, file = paste0("./11d/result_lassobss",seed,".RData"))




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

# decide s_list values for iht
s_list = c(c(1,2,3,4), seq(5,100,5))

## generate data
data_i <- data_lm(n = n, mu_x = rep(0,p), S_x = S_x,
                  s = s, betamin = betamin, sigma_e = sigma_e)
X = data_i$X
Y = data_i$Y  
beta = data_i$beta
epsilon = data_i$epsilon
supp_true = data_i$supp_true

## get the path
lassobss_result = lassobss(X, Y, supp_true = supp_true, s_list = s_list, scr_dim = 300, timelimit = 600)
lassobss_TPR= lassobss_result$TPR
lassobss_FDR = lassobss_result$FDR


save(beta, epsilon, supp_true, 
     lassobss_TPR, lassobss_FDR, file = paste0("./11e/result_lassobss",seed,".RData"))

















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


# decide s_list values for iht
s_list = seq(5,100,5)

## generate data

data_i <- data_lm(n = n, mu_x = rep(0,p), S_x = S_x,
                  s = s, betamin = betamin, sigma_e = sigma_e)

X = data_i$X
Y = data_i$Y  
beta = data_i$beta
epsilon = data_i$epsilon
supp_true = data_i$supp_true



## get the path
lassobss_result = lassobss(X, Y, supp_true = supp_true, s_list = s_list, scr_dim = 300, timelimit = 600)
lassobss_TPR= lassobss_result$TPR
lassobss_FDR = lassobss_result$FDR


save(beta, epsilon, supp_true, 
     lassobss_TPR, lassobss_FDR, file = paste0("./11f/result_lassobss",seed,".RData"))









set.seed(seed)



## 12a
# covariance for x
rho = 0.5
S_x = matrix(rho, nrow = p, ncol = p) + diag(1-rho, nrow =p, ncol = p)

# signal and noise
betamin = 0.1
sigma_e = 0.5

# decide s_list values for iht
s_list = c(c(1,2,3,4), seq(5,100,5))

## generate data
data_i <- data_lm(n = n, mu_x = rep(0,p), S_x = S_x,
                  s = s, betamin = betamin, sigma_e = sigma_e)
X = data_i$X
Y = data_i$Y  
beta = data_i$beta
epsilon = data_i$epsilon
supp_true = data_i$supp_true

## get the path
lassobss_result = lassobss(X, Y, supp_true = supp_true, s_list = s_list, scr_dim = 300, timelimit = 600)
lassobss_TPR= lassobss_result$TPR
lassobss_FDR = lassobss_result$FDR


save(beta, epsilon, supp_true, 
     lassobss_TPR, lassobss_FDR, file = paste0("./12a/result_lassobss",seed,".RData"))



set.seed(seed)

## 12b
# covariance for x
rho = 0.8
S_x = matrix(rho, nrow = p, ncol = p) + diag(1-rho, nrow =p, ncol = p)

# signal and noise
betamin = 0.1
sigma_e = 0.5

# decide s_list values for iht
s_list = c(c(1,2,3,4), seq(5,100,5))

## generate data
data_i <- data_lm(n = n, mu_x = rep(0,p), S_x = S_x,
                  s = s, betamin = betamin, sigma_e = sigma_e)
X = data_i$X
Y = data_i$Y  
beta = data_i$beta
epsilon = data_i$epsilon
supp_true = data_i$supp_true

## get the path
lassobss_result = lassobss(X, Y, supp_true = supp_true, s_list = s_list, scr_dim = 300, timelimit = 600)
lassobss_TPR= lassobss_result$TPR
lassobss_FDR = lassobss_result$FDR


save(beta, epsilon, supp_true, 
     lassobss_TPR, lassobss_FDR, file = paste0("./12b/result_lassobss",seed,".RData"))








set.seed(seed)

## 12c
# covariance for x
rho = 0.5
S_x = matrix(rho, nrow = p, ncol = p) + diag(1-rho, nrow =p, ncol = p)

# signal and noise
betamin = 0.1
sigma_e = 1


# decide s_list values for iht
s_list = c(c(1,2,3,4), seq(5,100,5))

## generate data

data_i <- data_lm(n = n, mu_x = rep(0,p), S_x = S_x,
                  s = s, betamin = betamin, sigma_e = sigma_e)

X = data_i$X
Y = data_i$Y  
beta = data_i$beta
epsilon = data_i$epsilon
supp_true = data_i$supp_true




## get the path
lassobss_result = lassobss(X, Y, supp_true = supp_true, s_list = s_list, scr_dim = 300, timelimit = 600)
lassobss_TPR= lassobss_result$TPR
lassobss_FDR = lassobss_result$FDR


save(beta, epsilon, supp_true, 
     lassobss_TPR, lassobss_FDR, file = paste0("./12c/result_lassobss",seed,".RData"))








set.seed(seed)

## 12d
# covariance for x
rho = 0.8
S_x = matrix(rho, nrow = p, ncol = p) + diag(1-rho, nrow =p, ncol = p)

# signal and noise
betamin = 0.1
sigma_e = 1


# decide s_list values for iht
s_list = c(c(1,2,3,4), seq(5,100,5))

## generate data

data_i <- data_lm(n = n, mu_x = rep(0,p), S_x = S_x,
                  s = s, betamin = betamin, sigma_e = sigma_e)

X = data_i$X
Y = data_i$Y  
beta = data_i$beta
epsilon = data_i$epsilon
supp_true = data_i$supp_true




## get the path
lassobss_result = lassobss(X, Y, supp_true = supp_true, s_list = s_list, scr_dim = 300, timelimit = 600)
lassobss_TPR= lassobss_result$TPR
lassobss_FDR = lassobss_result$FDR


save(beta, epsilon, supp_true, 
     lassobss_TPR, lassobss_FDR, file = paste0("./12d/result_lassobss",seed,".RData"))






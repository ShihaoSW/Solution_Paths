fm_cov <- function(p, eig_list){
    ## generate a covariance matrix from the factor model
    ## Input: p: dimension
    ##        eig_list: vector of eigenvalues of the low rank component Sigma_b
    ## Output: Sigma_f: covariance matrix (Sigma_f = Sigma_b + diag(p))
    
    # construct Sigma_b
    K = length(eig_list)
    V_0 = matrix(rnorm(n = p * K, mean = 0, sd = 1), nrow = p, ncol = K)
    V = qr.Q(qr(V_0))
    Sigma_b = V %*% diag(eig_list) %*% t(V)
    # construct Sigma_u
    Sigma_u = diag(p)
    # obtain covariance
    Sigma_f = Sigma_b + Sigma_u
    return(Sigma_f)
}

data_lm <- function(n, mu_x, S_x, s, betamin, sigma_e){
    ## generate linear model data
    ## Input: n: sample size
    ##        mu_x: mean of of each observation x_i
    ##        S_x: Covariance matrix of x_i
    ##        beta_min: minimum possible absolute value of entries of beta
    ##        sigma_e: standard deviation of the noise
    ## Output: data_lm: a list containing X, Y, beta, supp_true (true support), epsilon (noise)
    p = length(mu_x)
    X = mvrnorm(n, mu = mu_x, Sigma = S_x)
    beta0 = as.matrix(1 + rnorm(s)^2) * betamin 
    beta = matrix(0, nrow = p, ncol = 1)
    supp_true = sample(p, s)
    beta[supp_true] = beta0
    epsilon = as.matrix(rnorm(n, mean = 0, sd = sigma_e))
    Y = X %*% beta + epsilon
    data_lm = list(X = X, Y = Y, beta = beta, supp_true = supp_true, epsilon = epsilon)
}

iht <- function(X, Y, scale_X = TRUE, s, beta0, L = FALSE, maxiter = 1e3, prec = 1e-7){
  ## the one stage iterative hard threshoding algorithm
  ## Input: X, Y: data
  ##        scale_X (boolean): whether to scale X before training
  ##        s: the restricted sparsity
  ##        beta0: the initial value of beta
  ##        L: the approxinate maximal eigenvalue of XTX/n
  ##        maxiter: maximum number of iterations
  ##        prec: precision where the iterations stop
  n = dim(X)[1]
  p = dim(X)[2]
  if (scale_X)
    X = scale(X)
  Sxx = t(X) %*% X / n
  Sxy = t(X) %*% Y / n
  if (L == FALSE){
    vals = eigen(Sxx)$values
    L = vals[which(vals == max(vals))]
  }
  
  t = 0
  betat = beta0
  eta = 2 / (3 * L)  # step length
  while (t < maxiter) {
    beta_old = betat
    grad = -Sxy + Sxx %*% beta_old
    descent = beta_old - eta * grad
    descent_sort = sort(abs(descent), decreasing=TRUE)
    
    betat = rep(0,p)
    proj = which(abs(descent) >= descent_sort[s])
    
    betat[proj] = descent[proj]
    
    # identify convergence condition
    if (sum((betat - beta_old)^2) < prec * (sum((beta_old)^2) + 1)) 
      break
    t = t + 1
  }
  
  iht_result = list(beta = betat, step = t)
  return(iht_result)
}

iht1_select <- function(X, Y, scale_X = TRUE, supp_true, pi_list = c(), L = FALSE, 
                        maxiter = 1e4, prec = 1e-7,  nfold = 10){
  ## variable selection using IHT (prjected gradient descent, one stage)
  ## Input: X, Y: data
  ##        supp_true: the true support (not used during training)
  ##        pi_list: list of projection size for training
  ##        L: the approxinate maximal eigenvalue of XTX
  ##        maxiter: maximal number of iterations in IHT
  ##        prec: precision of stopping criterion in IHT
  ## Output: scad_result: list of TPR, FDR and their values at the best lambda value
  n = dim(X)[1]
  p = dim(X)[2]
  spars = length(supp_true)
  
  ## begin to do iterative hard thresholding
  if (scale_X)
    X = scale(X)
  Sxx = t(X) %*% X / n
  if (L == FALSE){
    vals = eigen(Sxx)$values
    L = vals[which(vals == max(vals))]
  }
  
  if(length(pi_list) == 0){
    pi_list = seq(2, ceiling(n/4) * 2, 2)
  }
  
  ## create beta table
  beta_iht = matrix(0, nrow = p, ncol = length(pi_list))
  step_iht = rep(0,length(pi_list))
  for (i in 1:length(pi_list)) {
    print(i)
    if (i == 1){
      res = iht(X = X, Y = Y, scale_X = scale_X, s = pi_list[1], L = L,
                beta0 = rep(0,p), maxiter = maxiter, prec = prec)
      beta_iht[,i] = res$beta
      step_iht[i] = res$step
    }
    else{
      res = iht(X = X, Y = Y, scale_X = scale_X, s = pi_list[i], L = L,
                beta0 = beta_iht[,i-1], maxiter = maxiter, prec = prec)
      beta_iht[,i] = res$beta
      step_iht[i] = res$step
    }
  }
  beta_supp = apply(beta_iht, 2, function(c) sum(c != 0))
  TP = apply(beta_iht, 2, function(c) length(intersect(which(c != 0), supp_true)))
  FP = beta_supp - TP
  TPR = TP / spars
  # deal with cases with zero supp
  beta_supp1 = beta_supp
  beta_supp1[which(beta_supp1 == 0)] = 1
  FDR = FP / beta_supp1
  
  
  ## cross validation 
  cv_score = rep(0, length(pi_list))
  n0 = floor(n / nfold) * nfold
  ind_sample = matrix(sample(n0), nrow = nfold, ncol = n0 / nfold)
  for (j in 1:nfold) {
    ind_test = ind_sample[j, ]
    ind_train = setdiff(c(1:n0), ind_test)
    X_test = X[ind_test, ]
    Y_test = Y[ind_test, ]
    X_train = X[ind_train, ]
    Y_train = Y[ind_train, ]
    beta_cv_iht = matrix(0, nrow = p, ncol = length(pi_list))
    for (k in 1:length(pi_list)) {
      cat(j,k,'\n')
      if (k == 1){
        beta_cv_iht[,1] = iht(X = X, Y = Y, scale_X = scale_X, s = pi_list[1], L = L,
                              beta0 = rep(0,p), maxiter = maxiter, prec = prec)$beta
      }
      else{
        beta_cv_iht[,k] = iht(X = X, Y = Y, scale_X = scale_X, s = pi_list[k], L = L,
                                    beta0 = beta_iht[,k-1], maxiter = maxiter, prec = prec)$beta
      }
    }
    for (k in 1:length(pi_list)) {
      cv_score[k] = cv_score[k] + sum((Y_test - X_test %*% beta_cv_iht[, k])^2)/nfold 
    }
  }
  ind_cv = which(cv_score == min(cv_score))[1]
  TPR_cv = TPR[ind_cv]
  FDR_cv = FDR[ind_cv]
  
  iht_result = list(TPR = TPR, FDR = FDR, TPR_cv = TPR_cv, FDR_cv = FDR_cv, steps = step_iht)
  return(iht_result)
}

iht_cor2 <- function(X, Y, scale_X = TRUE, s_list, g_list, beta0, maxiter = 1e3, prec = 1e-7){
    ## the iterative hard thresholding algorithm (two stage)
    ## Input: X, Y: data
    ##        scale_X (boolean): whether to scale X before training
    ##        s_list: list of possible values of parameter s
    ##        g_list (same length as s_list): list of possible values of parameter g
    ##        beta0: initial value of beta
    ##        maxiter: maximum number of iterations 
    ##        prec: precision where the iterations stop
    ## Output: iht_cor_result: list containing beta (coefficients), grad (gradient), steps_list (number of iterations performed)
    n = dim(X)[1]
    p = dim(X)[2]
    if (scale_X) 
        X = scale(X)
    Sxy = t(X) %*% Y / n
    s_list_len = length(s_list)
    steps_list = rep(0, s_list_len)
    beta_iht_cor = matrix(0, nrow=p, ncol=s_list_len)
    grad_iht_cor = matrix(0, nrow=p, ncol=s_list_len)

    for (k in 1:s_list_len){
        s_iht = s_list[k]
        g_iht = g_list[k]
        t = 0
        betat = beta0
        while (t < maxiter){
            beta_old = betat
            Sxbeta = X %*% beta_old
            if(is.na(Sxbeta[1])){
              grad = -Sxy
            }
            else{
              grad = - Sxy + ( t(X) %*% Sxbeta / n )  
            }
            indt1 = which(beta_old != 0)
            grad_1 = grad
            grad_sort = sort(abs(grad_1), decreasing=TRUE)
            indt2 = which(abs(grad_1) >= grad_sort[g_iht])
            indt = union(indt1, indt2)

            # refit 1
            Xt = X[, indt]
            betat = rep(0, p)
            betat[indt] = solve(t(Xt) %*% Xt) %*% (t(Xt) %*% Y)
            
            # truncation 
            betat_sort = sort(abs(betat), decreasing=TRUE)
            if ((betat_sort[s_iht] == 0) & (s_iht >= g_iht)){
              indt0 = which(betat_sort != 0)
            }
            else{
              indt0 = which(abs(betat) >= betat_sort[s_iht])
            }

            
            # refit 2
            Xt0 = X[, indt0]
            betat = rep(0, p)
            betat[indt0] = solve(t(Xt0) %*% Xt0) %*% (t(Xt0) %*% Y)
            
            # identify convergence condition
            if (sum((betat - beta_old)^2) < prec * (sum((beta_old)^2) + 1)) 
                break
            t = t + 1
        }
        beta_iht_cor[, k] = betat
        Sxbeta = X %*% betat
        grad_iht_cor[, k] = - Sxy + (t(X) %*% Sxbeta / n )
        steps_list[k] = t
    }

    iht_cor_result = list(beta=beta_iht_cor, grad= grad_iht_cor, steps_list=steps_list)
    return(iht_cor_result)
}

lambda_select_scad <- function(n, mu_x, S_x, s, betamin, sigma_e, nlambda = 100, lambda_ratio = 1.2, ntest){
    ## select a list of lambda for SCAD
    ## Input: mu_x: mean of x_i
    ##        S_x: variance of x_i
    ##        s: sparsity
    ##        betamin: minimum possible absolute value of beta
    ##        sigma_e: standard deviation of noise
    ##        nlambda: number of lambda's to be chosen
    ##        lambda_ratio: the ratio between adjacent lambda's
    ##        ntest: number of tests to perform to evaluate the range of lambda
    ## Output: lambda_list: list of lambda
    lambda_max = 0
    for (ntest in 1:ceiling(nrep/10)){
        data_i <- data_lm(n = n, mu_x = rep(0,p), S_x = S_x,
                          s = s, betamin = betamin, sigma_e = sigma_e)
        X = data_i$X
        Y = data_i$Y
        model_scad = picasso(X, Y, nlambda=nlambda, 
                             method="scad", intercept=FALSE)
        lambda_scad = model_scad$lambda
        lambda_max = max(max(lambda_scad),lambda_max)
    }
    lambda_max = lambda_max * 2
    ratios = lambda_ratio ^ seq(0, nlambda-1)
    lambda_list = lambda_max / ratios
    return(lambda_list)
}

lambda_select_scad_g <- function(X, Y, nlambda =100, lambda_ratio = 1.2, ntest){
  ## generate candidate lambda list for scad from only X and Y
  lambda_max = 0
  for (ntest in 1:ceiling(nrep/10)){
    model_scad = picasso(X, Y, nlambda=nlambda, 
                         method="scad", intercept=FALSE)
    lambda_scad = model_scad$lambda
    lambda_max = max(max(lambda_scad),lambda_max)
  }
  lambda_max = lambda_max * 2
  ratios = lambda_ratio ^ seq(0, nlambda-1)
  lambda_list = lambda_max / ratios
  return(lambda_list)
}

lambda_select_lasso <- function(n, mu_x, S_x, s, betamin, sigma_e, nlambda = 100, lambda_ratio = 1.2, ntest){
    ## select a list of lambda for LASSO
    ## Input: mu_x: mean of x_i
    ##        S_x: variance of x_i
    ##        s: sparsity
    ##        betamin: minimum possible absolute value of beta
    ##        sigma_e: standard deviation of noise
    ##        nlambda: number of lambda's to be chosen
    ##        lambda_ratio: the ratio between adjacent lambda's
    ##        ntest: number of tests to perform to evaluate the range of lambda
    ## Output: lambda_list: list of lambda
    lambda_max = 0
    for (ntest in 1:ceiling(nrep/10)){
        data_i <- data_lm(n = n, mu_x = rep(0,p), S_x = S_x,
                          s = s, betamin = betamin, sigma_e = sigma_e)
        X = data_i$X
        Y = data_i$Y
        model_scad = picasso(X, Y, nlambda=nlambda, 
                             method="l1", intercept=FALSE)
        lambda_scad = model_scad$lambda
        lambda_max = max(max(lambda_scad),lambda_max)
    }
    lambda_max = lambda_max * 2
    ratios = lambda_ratio ^ seq(0, nlambda-1)
    lambda_list = lambda_max / ratios
    return(lambda_list)
}

lambda_select_lasso_g <- function(X, Y, nlambda = 100, lambda_ratio = 1.2, ntest){
  ## generate candidate lambda list for lasso from only X and Y
  lambda_max = 0
  for (ntest in 1:ceiling(nrep/10)){
    model_scad = picasso(X, Y, nlambda=nlambda, 
                         method="l1", intercept=FALSE)
    lambda_scad = model_scad$lambda
    lambda_max = max(max(lambda_scad),lambda_max)
  }
  lambda_max = lambda_max * 2
  ratios = lambda_ratio ^ seq(0, nlambda-1)
  lambda_list = lambda_max / ratios
  return(lambda_list)
}

lambda_select_mcp <- function(n, mu_x, S_x, s, betamin, sigma_e, nlambda = 100, lambda_ratio = 1.2, ntest){
  ## select a list of lambda for MCP
  ## Input: mu_x: mean of x_i
  ##        S_x: variance of x_i
  ##        s: sparsity
  ##        betamin: minimum possible absolute value of beta
  ##        sigma_e: standard deviation of noise
  ##        nlambda: number of lambda's to be chosen
  ##        lambda_ratio: the ratio between adjacent lambda's
  ##        ntest: number of tests to perform to evaluate the range of lambda
  ## Output: lambda_list: list of lambda
  lambda_max = 0
  for (ntest in 1:ceiling(nrep/10)){
    data_i <- data_lm(n = n, mu_x = rep(0,p), S_x = S_x,
                      s = s, betamin = betamin, sigma_e = sigma_e)
    X = data_i$X
    Y = data_i$Y
    model_mcp = picasso(X, Y, nlambda=nlambda, 
                        method="mcp", intercept=FALSE)
    lambda_mcp = model_mcp$lambda
    lambda_max = max(max(lambda_mcp),lambda_max)
  }
  lambda_max = lambda_max * 2
  ratios = lambda_ratio ^ seq(0, nlambda-1)
  lambda_list = lambda_max / ratios
  return(lambda_list)
}

lambda_select_mcp_g <- function(X, Y, nlambda = 100, lambda_ratio = 1.2, ntest){
  ## generate candidate lambda list for mcp from only X and Y
  lambda_max = 0
  for (ntest in 1:ceiling(nrep/10)){
    model_mcp = picasso(X, Y, nlambda=nlambda, 
                        method="mcp", intercept=FALSE)
    lambda_mcp = model_mcp$lambda
    lambda_max = max(max(lambda_mcp),lambda_max)
  }
  lambda_max = lambda_max * 2
  ratios = lambda_ratio ^ seq(0, nlambda-1)
  lambda_list = lambda_max / ratios
  return(lambda_list)
}

scad_select <- function(X, Y, supp_true, lambda_list, nfold = 10){
    ## variable selection using SCAD
    ## Input: X, Y: data
    ##        supp_true: the true support (not used during training)
    ##        lambda_list: list of lambda values for training
    ##        nfold: number of folds for cross validation
    ## Output: scad_result: list of TPR, FDR and their values at the best lambda value
    s = length(supp_true)
    n = dim(X)[1]
    p = dim(X)[2]
    
    ## TPR and FDR
    model_scad = picasso(X, Y, lambda=lambda_list, 
                         method="scad", intercept=FALSE)
    beta_scad = model_scad$beta
    beta_supp = apply(beta_scad, 2, function(c) sum(c != 0))
    TP = apply(beta_scad, 2, function(c) length(intersect(which(c != 0), supp_true)))
    FP = beta_supp - TP
    TPR = TP / s
    # deal with cases with zero supp
    beta_supp1 = beta_supp
    beta_supp1[which(beta_supp1 == 0)] = 1
    FDR = FP / beta_supp1
    
    ## cross validation 
    cv_score = rep(0, length(lambda_list))
    n0 = floor(n / nfold) * nfold
    ind_sample = matrix(sample(n0), nrow = nfold, ncol = n0 / nfold)
    for (k in 1:length(lambda_list)){
        lambda = lambda_list[k]
        for (j in 1:nfold){
            ind_test = ind_sample[j, ]
            ind_train = setdiff(c(1:n0), ind_test)
            X_test = X[ind_test, ]
            Y_test = Y[ind_test, ]
            X_train = X[ind_train, ]
            Y_train = Y[ind_train, ]
            model_train = picasso(X_train, Y_train, lambda=lambda, 
                                  method="scad", intercept=FALSE)
            beta_train = model_train$beta
            cv_score[k] = cv_score[k] + sum((Y_test - X_test %*% beta_train)^2)/nfold 
        }
    }
    ind_cv = which(cv_score == min(cv_score))
    TPR_cv = TPR[ind_cv]
    FDR_cv = FDR[ind_cv]

    scad_result = list(TPR = TPR, FDR = FDR, TPR_cv = TPR_cv, FDR_cv = FDR_cv)
    return(scad_result)
}

lasso_select <- function(X, Y, supp_true, lambda_list, nfold = 10){
    ## variable selection using LASSO
    ## Input: X, Y: data
    ##        supp_true: the true support (not used during training)
    ##        lambda_list: list of lambda values for training
    ##        nfold: number of folds for cross validation
    ## Output: scad_result: list of TPR, FDR and their values at the best lambda value
  
    s = length(supp_true)
    n = dim(X)[1]
    p = dim(X)[2]
  
    ## TPR and FDR
    model_lasso = picasso(X, Y, lambda=lambda_list, 
                          method="l1", intercept=FALSE)
    beta_lasso = model_lasso$beta
    beta_supp = apply(beta_lasso, 2, function(c) sum(c != 0))
    TP = apply(beta_lasso, 2, function(c) length(intersect(which(c != 0), supp_true)))
    FP = beta_supp - TP
    TPR = TP / s
    # deal with cases with zero supp
    beta_supp1 = beta_supp
    beta_supp1[which(beta_supp1 == 0)] = 1
    FDR = FP / beta_supp1
  
    ## cross validation 
    cv_score = rep(0, length(lambda_list))
    n0 = floor(n / nfold) * nfold
    ind_sample = matrix(sample(n0), nrow = nfold, ncol = n0 / nfold)
    for (k in 1:length(lambda_list)){
        lambda = lambda_list[k]
        for (j in 1:nfold){
            ind_test = ind_sample[j, ]
            ind_train = setdiff(c(1:n0), ind_test)
            X_test = X[ind_test, ]
            Y_test = Y[ind_test, ]
            X_train = X[ind_train, ]
            Y_train = Y[ind_train, ]
            model_train = picasso(X_train, Y_train, lambda=lambda, 
                                  method="l1", intercept=FALSE)
            beta_train = model_train$beta
            cv_score[k] = cv_score[k] + sum((Y_test - X_test %*% beta_train)^2)/nfold 
        }
    }
    ind_cv = which(cv_score == min(cv_score))
    TPR_cv = TPR[ind_cv]
    FDR_cv = FDR[ind_cv]
  
    lasso_result = list(TPR = TPR, FDR = FDR, TPR_cv = TPR_cv, FDR_cv = FDR_cv)
    return(lasso_result)
}

mcp_select <- function(X, Y, supp_true, lambda_list, nfold = 10){
  ## variable selection using SCAD
  ## Input: X, Y: data
  ##        supp_true: the true support (not used during training)
  ##        lambda_list: list of lambda values for training
  ##        nfold: number of folds for cross validation
  ## Output: scad_result: list of TPR, FDR and their values at the best lambda value
  s = length(supp_true)
  n = dim(X)[1]
  p = dim(X)[2]
  
  ## TPR and FDR
  model_mcp = picasso(X, Y, lambda=lambda_list, 
                      method="mcp", intercept=FALSE)
  beta_mcp = model_mcp$beta
  beta_supp = apply(beta_mcp, 2, function(c) sum(c != 0))
  TP = apply(beta_mcp, 2, function(c) length(intersect(which(c != 0), supp_true)))
  FP = beta_supp - TP
  TPR = TP / s
  # deal with cases with zero supp
  beta_supp1 = beta_supp
  beta_supp1[which(beta_supp1 == 0)] = 1
  FDR = FP / beta_supp1
  
  ## cross validation 
  cv_score = rep(0, length(lambda_list))
  n0 = floor(n / nfold) * nfold
  ind_sample = matrix(sample(n0), nrow = nfold, ncol = n0 / nfold)
  for (k in 1:length(lambda_list)){
    lambda = lambda_list[k]
    for (j in 1:nfold){
      ind_test = ind_sample[j, ]
      ind_train = setdiff(c(1:n0), ind_test)
      X_test = X[ind_test, ]
      Y_test = Y[ind_test, ]
      X_train = X[ind_train, ]
      Y_train = Y[ind_train, ]
      model_train = picasso(X_train, Y_train, lambda=lambda, 
                            method="mcp", intercept=FALSE)
      beta_train = model_train$beta
      cv_score[k] = cv_score[k] + sum((Y_test - X_test %*% beta_train)^2)/nfold 
    }
  }
  ind_cv = which(cv_score == min(cv_score))
  TPR_cv = TPR[ind_cv]
  FDR_cv = FDR[ind_cv]
  
  mcp_result = list(TPR = TPR, FDR = FDR, TPR_cv = TPR_cv, FDR_cv = FDR_cv)
  return(mcp_result)
}


cosbss <- function(X, Y, scale_X = TRUE, supp_true, s_list = c(), proj_size,
                   expan_size, timelimit = 3600, maxiter = 1e3, prec = 1e-7){
  ## variable selection using CoSaMP plus best subset selection
  ## Input: X, Y: data
  ##        supp_true: the true support (not used during training)
  ##        s_list: list of estimated sparsity for training
  ##        proj_size: projection size in cosamp
  ##        expan_size: expansion size in cosamp
  ##        timelimit: time limit for each round of bss
  ##        maxiter: maximal number of iterations in cosamp
  ##        prec: precision of stopping criterion in cosamp
  ## Output: cosbss_result: list of TPR, FDR
  n = dim(X)[1]
  p = dim(X)[2]
  s = length(supp_true)
  
  # screening step
  scr_beta = iht_cor2(X, Y, scale_X = scale_X, s_list = proj_size,
           g_list = expan_size, beta0 = rep(0,p), maxiter = maxiter, prec = prec)$beta
  
  scr_supp = which(scr_beta != 0)
  
  # screened design
  X_scr = X[,scr_supp]
  
  # create beta_table
  beta_table = matrix(0, nrow = p, ncol = length(s_list))
  
  # do bss
  beta_on_scr = coef( bs(X_scr, Y, intercept=FALSE, k=s_list, verbose=TRUE, time.limit = timelimit) )
  beta_table[scr_supp, ] = beta_on_scr
  
  # FDR and TPR
  beta_supp = apply(beta_table, 2, function(c) sum(c != 0))
  TP = apply(beta_table, 2, function(c) length(intersect(which(c != 0), supp_true)))
  FP = beta_supp - TP
  TPR = TP / s
  # deal with cases with zero supp
  beta_supp1 = beta_supp
  beta_supp1[which(beta_supp1 == 0)] = 1
  FDR = FP / beta_supp1
  
  cosbss_result = list(TPR = TPR, FDR = FDR)
  return(cosbss_result)
}

L0Learn_grid <- function(X, Y, supp_true, method = 'L0L2', max_size){
  ## variable selection (with multiple solutions in grid) via package 'L0Learn'
  # prepare
 
  n = dim(X)[1]
  p = dim(X)[2]
  s = length(supp_true)
  
  fit <- L0Learn.fit(X, Y, penalty = method, maxSuppSize = max_size)    

  lambda_ns = lengths(fit$lambda)
  gamma_size = length(lengths(fit$lambda))
  sol_num = sum(lengths(fit$lambda))
  
  beta_table = matrix(0, nrow = p, ncol = sol_num)
  
  col_index = 0
  
  for (i in 1:gamma_size) {
    gamma = fit$gamma[i]
    for (j in 1:lambda_ns[i]) {
      col_index = col_index + 1
      lambda = unlist(fit$lambda[i])[j]
      beta_table[,col_index] = coef(fit, lambda = lambda, gamma = gamma)[2:(p+1)]
    }
  }
  
  # FDR and TPR
  beta_supp = apply(beta_table, 2, function(c) sum(c != 0))
  TP = apply(beta_table, 2, function(c) length(intersect(which(c != 0), supp_true)))
  FP = beta_supp - TP
  TPR = TP / s
  # deal with cases with zero supp
  beta_supp1 = beta_supp
  beta_supp1[which(beta_supp1 == 0)] = 1
  FDR = FP / beta_supp1
  
  cosbss_result = list(TPR = TPR, FDR = FDR)
  return(cosbss_result)
}



L0Learn_path <- function(X, Y, supp_true, method = 'L0L2', max_size, lambda_ratio = 1.1, nlambda = 100){
  ## variable selection (with multiple solutions on a path) via package 'L0Learn'
  # prepare
  
  n = dim(X)[1]
  p = dim(X)[2]
  s = length(supp_true)
  
  first_fit <- L0Learn.fit(X, Y, penalty = method)
  lambda_max = max(unlist(first_fit$lambda)) * 3
  ratios = lambda_ratio ^ seq(0, nlambda-1)
  lambda_list = list(lambda_max / ratios)
  gamma_input = min(first_fit$gamma)
  
  fit <- L0Learn.fit(X, Y, penalty = method, maxSuppSize = max_size, nGamma = 1, gammaMax = gamma_input, 
                     gammaMin = gamma_input, lambdaGrid = lambda_list)   
  
  
  lambda_ns = lengths(fit$lambda)
  gamma_size = length(lengths(fit$lambda))
  sol_num = sum(lengths(fit$lambda))
  
  beta_table = matrix(0, nrow = p, ncol = sol_num)
  
  col_index = 0
  
  for (i in 1:gamma_size) {
    gamma = fit$gamma[i]
    for (j in 1:lambda_ns[i]) {
      col_index = col_index + 1
      lambda = unlist(fit$lambda[i])[j]
      beta_table[,col_index] = coef(fit, lambda = lambda, gamma = gamma)[2:(p+1)]
    }
  }
  
  # FDR and TPR
  beta_supp = apply(beta_table, 2, function(c) sum(c != 0))
  TP = apply(beta_table, 2, function(c) length(intersect(which(c != 0), supp_true)))
  FP = beta_supp - TP
  TPR = TP / s
  # deal with cases with zero supp
  beta_supp1 = beta_supp
  beta_supp1[which(beta_supp1 == 0)] = 1
  FDR = FP / beta_supp1
  
  cosbss_result = list(TPR = TPR, FDR = FDR)
  return(cosbss_result)
}


L0Learn_path2 <- function(X, Y, supp_true, method = 'L0L2', max_size, lambda_list, gamma_input){
  ## variable selection (with multiple solutions on a path) via package 'L0Learn'
  ## lambda and gamma are predetermined
  # prepare
  
  n = dim(X)[1]
  p = dim(X)[2]
  s = length(supp_true)
  
  
  fit <- L0Learn.fit(X, Y, penalty = method, maxSuppSize = max_size, nGamma = 1, gammaMax = gamma_input, 
                     gammaMin = gamma_input, lambdaGrid = lambda_list, maxIters = 2000, rtol = 1e-09)   
  
  
  lambda_ns = lengths(fit$lambda)
  gamma_size = length(lengths(fit$lambda))
  sol_num = sum(lengths(fit$lambda))
  
  beta_table = matrix(0, nrow = p, ncol = sol_num)
  
  col_index = 0
  
  for (i in 1:gamma_size) {
    gamma = fit$gamma[i]
    for (j in 1:lambda_ns[i]) {
      col_index = col_index + 1
      lambda = unlist(fit$lambda[i])[j]
      beta_table[,col_index] = coef(fit, lambda = lambda, gamma = gamma)[2:(p+1)]
    }
  }
  
  # FDR and TPR
  beta_supp = apply(beta_table, 2, function(c) sum(c != 0))
  TP = apply(beta_table, 2, function(c) length(intersect(which(c != 0), supp_true)))
  FP = beta_supp - TP
  TPR = TP / s
  # deal with cases with zero supp
  beta_supp1 = beta_supp
  beta_supp1[which(beta_supp1 == 0)] = 1
  FDR = FP / beta_supp1
  
  cosbss_result = list(TPR = TPR, FDR = FDR)
  return(cosbss_result)
}



lasso_screen <- function(X, Y, scr_dim, nlambda = 100, lambda_ratio = 1.1){
  model_lasso = picasso(X, Y, nlambda=nlambda, method="l1", intercept=FALSE)
  lambda_lasso = model_lasso$lambda
  lambda_max = lambda_lasso * 3
  ratios = lambda_ratio ^ seq(0, nlambda-1)
  lambda_list = lambda_max / ratios
  
  model_lasso = picasso(X, Y, lambda=lambda_list, method="l1", intercept=FALSE)
  
  count_list = rep(0,nlambda)
  for (i in 1:nlambda){
    count_list[i] = length(which(model_lasso$beta[,i] != 0))
  }
  dif_list = abs(count_list - scr_dim)
  selected = which(dif_list == min(dif_list))[1] 
  beta = model_lasso$beta[, selected]
  supp = which(beta != 0)
  spars = count_list[selected]
  result = list(beta = beta, supp = supp, spars = spars)
  return(result)
}

lassobss <- function(X, Y, supp_true, s_list, scr_dim, timelimit, maxiter, prec){
  ## variable selection using CoSaMP plus best subset selection
  ## Input: X, Y: data
  ##        supp_true: the true support (not used during training)
  ##        s_list: list of estimated sparsity for training
  ##        scr_dim: screened dimension
  ##        timelimit: time limit for each round of bss
  ##        maxiter: maximal number of iterations in cosamp
  ##        prec: precision of stopping criterion in cosamp
  ## Output: cosbss_result: list of TPR, FDR
  n = dim(X)[1]
  p = dim(X)[2]
  s = length(supp_true)
  
  # screening step
  scr_beta = lasso_screen(X, Y, scr_dim = scr_dim)$beta
  
  scr_supp = which(scr_beta != 0)
  
  # screened design
  X_scr = X[,scr_supp]
  
  # create beta_table
  beta_table = matrix(0, nrow = p, ncol = length(s_list))
  
  # do bss
  beta_on_scr = coef( bs(X_scr, Y, intercept=FALSE, k=s_list, verbose=TRUE, time.limit = timelimit) )
  beta_table[scr_supp, ] = beta_on_scr
  
  # FDR and TPR
  beta_supp = apply(beta_table, 2, function(c) sum(c != 0))
  TP = apply(beta_table, 2, function(c) length(intersect(which(c != 0), supp_true)))
  FP = beta_supp - TP
  TPR = TP / s
  # deal with cases with zero supp
  beta_supp1 = beta_supp
  beta_supp1[which(beta_supp1 == 0)] = 1
  FDR = FP / beta_supp1
  
  cosbss_result = list(TPR = TPR, FDR = FDR)
  return(cosbss_result)
}


iht_m_select <- function(X, Y, scale_X = TRUE, supp_true, pi_list = c(), nlambda = 100, 
                       lambda_ratio = 1.2, maxiter = 1e3, prec = 1e-7, nfold = 10){
  ## variable selection using IHT
  ## Input: X, Y: data
  ##        supp_true: the true support (not used during training)
  ##        pi_list: list of projection size for training
  ##        nlambda: number of lambda values in mcp
  ##        lambda_ratio: determine the gap between lambda values
  ##        maxiter: maximal number of iterations in IHT
  ##        prec: precision of stopping criterion in IHT
  ##        nfold: number of folds for cross validation
  ## Output: scad_result: list of TPR, FDR and their values at the best lambda value
  n = dim(X)[1]
  p = dim(X)[2]
  s = length(supp_true)
  
  ## tune the expansion size from cv point of mcp
  # get a lambda list
  model_mcp = picasso(X, Y, nlambda=nlambda, 
                      method="mcp", intercept=FALSE)
  lambda_max = (model_mcp$lambda) * 2
  ratios = lambda_ratio ^ seq(0, nlambda-1)
  lambda_list = lambda_max / ratios

  # cross validation to search for expansion size
  cv_score = rep(0, length(lambda_list))
  n0 = floor(n / nfold) * nfold
  ind_sample = matrix(sample(n0), nrow = nfold, ncol = n0 / nfold)
  for (k in 1:length(lambda_list)){
    lambda = lambda_list[k]
    for (j in 1:nfold){
      ind_test = ind_sample[j, ]
      ind_train = setdiff(c(1:n0), ind_test)
      X_test = X[ind_test, ]
      Y_test = Y[ind_test, ]
      X_train = X[ind_train, ]
      Y_train = Y[ind_train, ]
      model_train = picasso(X_train, Y_train, lambda=lambda, 
                            method="mcp", intercept=FALSE)
      beta_train = model_train$beta
      cv_score[k] = cv_score[k] + sum((Y_test - X_test %*% beta_train)^2)/nfold 
    }
  }
  ind_cv = which(cv_score == min(cv_score))
  expan_size = length(which(picasso(X, Y, lambda=lambda_list[ind_cv], 
                                    method="mcp", intercept=FALSE)$beta != 0 ))
  
  
  
  if(length(pi_list) == 0){
    pi_list = seq(2, ceiling(n/4) * 2, 2)
  }
  
  if(expan_size >= n/2){
    expan_size = pi_list[floor(length(pi_list)/2)]
  }
  
  
  ## create beta table
  beta_iht = matrix(0, nrow = p, ncol = length(pi_list))
  for (i in 1:length(pi_list)) {
    cat(i)
    if (i == 1){
      beta_iht[,1] = iht_cor2(X, Y, scale_X = scale_X, s_list = pi_list[i],
                              g_list = expan_size, beta0 = rep(0,p), maxiter = maxiter, prec = prec)$beta
    }
    else{
      beta_iht[,i] = iht_cor2(X, Y, scale_X = scale_X, s_list = pi_list[i],
                              g_list = expan_size, beta0 = beta_iht[, (i-1)], maxiter = maxiter, prec = prec)$beta
    }
  }
  beta_supp = apply(beta_iht, 2, function(c) sum(c != 0))
  TP = apply(beta_iht, 2, function(c) length(intersect(which(c != 0), supp_true)))
  FP = beta_supp - TP
  TPR = TP / s
  # deal with cases with zero supp
  beta_supp1 = beta_supp
  beta_supp1[which(beta_supp1 == 0)] = 1
  FDR = FP / beta_supp1
  
  ## cross validation 
  cv_score = rep(0, length(pi_list))
  n0 = floor(n / nfold) * nfold
  ind_sample = matrix(sample(n0), nrow = nfold, ncol = n0 / nfold)
  for (j in 1:nfold) {
    ind_test = ind_sample[j, ]
    ind_train = setdiff(c(1:n0), ind_test)
    X_test = X[ind_test, ]
    Y_test = Y[ind_test, ]
    X_train = X[ind_train, ]
    Y_train = Y[ind_train, ]
    beta_cv_iht = matrix(0, nrow = p, ncol = length(pi_list))
    for (k in 1:length(pi_list)) {
      cat(j,k,'\n')
      if (k == 1){
        beta_cv_iht[,1] = iht_cor2(X_train, Y_train, scale_X = scale_X, s_list = pi_list[1],
                                   g_list = expan_size, beta0 = rep(0,p), maxiter = maxiter, prec = prec)$beta
      }
      else{
        beta_cv_iht[,k] = iht_cor2(X_train, Y_train, scale_X = scale_X, s_list = pi_list[k],
                                   g_list = expan_size, beta0 = beta_cv_iht[, (k-1)], maxiter = maxiter, prec = prec)$beta
      }
    }
    for (k in 1:length(pi_list)) {
      cv_score[k] = cv_score[k] + sum((Y_test - X_test %*% beta_cv_iht[, k])^2)/nfold 
    }
  }
  ind_cv = which(cv_score == min(cv_score))[1]
  TPR_cv = TPR[ind_cv]
  FDR_cv = FDR[ind_cv]
  
  iht_result = list(TPR = TPR, FDR = FDR, TPR_cv = TPR_cv, FDR_cv = FDR_cv)
  return(iht_result)
}




iht_select <- function(X, Y, scale_X = TRUE, supp_true, lambda_list = c(), 
                       maxiter = 1e3, prec = 1e-7, nfold = 10){
  ## variable selection using IHT
  ## Input: X, Y: data
  ##        supp_true: the true support (not used during training)
  ##        lambda_list: list of lambda values for training
  ##        maxiter: maximal number of iterations in IHT
  ##        prec: precision of stopping criterion in IHT
  ##        nfold: number of folds for cross validation
  ## Output: scad_result: list of TPR, FDR and their values at the best lambda value
  n = dim(X)[1]
  p = dim(X)[2]
  s = length(supp_true)
  
  if(length(lambda_list) == 0){
    lambda_list = seq(2, ceiling(n/4) * 2, 2)
  }
  
  expan_size = lambda_list[floor(length(lambda_list)/2)]
  
  ## create beta table
  beta_iht = matrix(0, nrow = p, ncol = length(lambda_list))
  for (i in 1:length(lambda_list)) {
    cat(i)
    if (i == 1){
      beta_iht[,1] = iht_cor2(X, Y, scale_X = scale_X, s_list = lambda_list[i],
                              g_list = expan_size, beta0 = rep(0,p), maxiter = maxiter, prec = prec)$beta
    }
    else{
      beta_iht[,i] = iht_cor2(X, Y, scale_X = scale_X, s_list = lambda_list[i],
                              g_list = expan_size, beta0 = beta_iht[, (i-1)], maxiter = maxiter, prec = prec)$beta
    }
  }
  beta_supp = apply(beta_iht, 2, function(c) sum(c != 0))
  TP = apply(beta_iht, 2, function(c) length(intersect(which(c != 0), supp_true)))
  FP = beta_supp - TP
  TPR = TP / s
  # deal with cases with zero supp
  beta_supp1 = beta_supp
  beta_supp1[which(beta_supp1 == 0)] = 1
  FDR = FP / beta_supp1
  
  ## cross validation 
  cv_score = rep(0, length(lambda_list))
  n0 = floor(n / nfold) * nfold
  ind_sample = matrix(sample(n0), nrow = nfold, ncol = n0 / nfold)
  for (j in 1:nfold) {
    ind_test = ind_sample[j, ]
    ind_train = setdiff(c(1:n0), ind_test)
    X_test = X[ind_test, ]
    Y_test = Y[ind_test, ]
    X_train = X[ind_train, ]
    Y_train = Y[ind_train, ]
    beta_cv_iht = matrix(0, nrow = p, ncol = length(lambda_list))
    for (k in 1:length(lambda_list)) {
      if (k == 1){
        beta_cv_iht[,1] = iht_cor2(X = X_train, Y = Y_train, scale_X = scale_X, s_list = lambda_list[1],
                                   g_list = expan_size, beta0 = rep(0,p), maxiter = maxiter, prec = prec)$beta
      }
      else{
        beta_cv_iht[,k] = iht_cor2(X = X_train, Y = Y_train, scale_X = scale_X, s_list = lambda_list[k],
                                   g_list = expan_size, beta0 = beta_cv_iht[, (k-1)], maxiter = maxiter, prec = prec)$beta
      }
    }
    for (k in 1:length(lambda_list)) {
      cv_score[k] = cv_score[k] + sum((Y_test - X_test %*% beta_cv_iht[, k])^2)/nfold 
    }
  }
  ind_cv = which(cv_score == min(cv_score))[1]
  TPR_cv = TPR[ind_cv]
  FDR_cv = FDR[ind_cv]
  
  iht_result = list(TPR = TPR, FDR = FDR, TPR_cv = TPR_cv, FDR_cv = FDR_cv)
  return(iht_result)
}

iht2_select <- function(X, Y, supp_true, s_list, g_list, beta0, maxiter = 1e3, prec = 1e-7){
    ## variable selection using iterative hard thresholding +
    ## Input: X, Y: data
    ##        supp_true: the true support (not used during training)
    ##        s_list: list of s values for training
    ##        g_list (same length as s_list): list of g values for training
    ##        beta0: initialization
    ##        maxiter: maximum number of iterations
    ##        prec: precision to determine when to stop the iterations
    ## Output: scad_result: list of TPR, FDR and number of iterations
  
    s = length(supp_true)
    p = dim(X)[2]
    model_iht2 = iht_cor2(X, Y, s_list=s_list, g_list=g_list, beta0=beta0, maxiter = maxiter, prec = prec)
    beta_iht2 = model_iht2$beta
    grad_iht2 = model_iht2$grad
    steps_iht2 = model_iht2$steps_list
    
    beta_order = matrix(0, nrow=p, ncol=length(s_list))
    grad_order = matrix(0, nrow=p, ncol=length(s_list))
    for (j in 1:length(s_list)){
        beta_order[, j] = rev(order(abs(beta_iht2[, j])))
        grad_iht2[beta_order[1:s_list[j], j], j] = 0
        grad_order[, j] = rev(order(abs(grad_iht2[, j])))
        beta_order[(s_list[j] + 1):p, j] = grad_order[1:(p - s_list[j]), j]
    }

    TP = matrix(0, nrow=p, ncol=length(s_list))
    FP = matrix(0, nrow=p, ncol=length(s_list))
    TPR = matrix(0, nrow=p, ncol=length(s_list))
    FDR = matrix(0, nrow=p, ncol=length(s_list))
    for (j in 1:(length(s_list))){
        TP_ind = as.numeric(is.element(beta_order[, j], supp_true))
        TP[, j] = cumsum(TP_ind)
        TPR[, j] = TP[, j] / s
        FP[, j] = (1:p) - TP[, j]
        FDR[, j] = FP[, j] / (1:p)
    }
    iht2_result = list(TPR=TPR, FDR=FDR, steps_iht2=steps_iht2)
    return(iht2_result)
}

ms_select <- function(X, Y, supp_true, scale_X = TRUE){
    ## variable selection using SIS (marginal screening)
    ## Input: X, Y: data
    ##        supp_true: the true support (not used during training)
    ##        scale_X: whether to scale X before training
    ## Output: ms_result: list of TPR and FDR 
  
    s = length(supp_true)
    n = dim(X)[1]
    p = dim(X)[2]
    if (scale_X) 
        X = scale(X)
    Sxx = t(X) %*% X / n
    Sxy = t(X) %*% Y / n
    grad = - Sxy
    grad_sort_idx = order(abs(grad), decreasing = TRUE)
    TP_ms = rep(0, p)
    TP_ms[1] = as.numeric(grad_sort_idx[1] %in% supp_true)
    for (i in c(2:p)){
        TP_ms[i] = TP_ms[i-1] + as.numeric(grad_sort_idx[i] %in% supp_true)
    }
    TPR_ms = TP_ms / s
    FP_ms = c(1:p) - TP_ms
    FDR_ms = rep(0, p)
    for (i in c(1:p)){
        FDR_ms[i] = FP_ms[i] / i
    }
    ms_result = list(TPR=TPR_ms, FDR=FDR_ms)
    return(ms_result)
}

scad_cv <- function(X, Y, lambda_list, nfold = 10){
    ## cross validation for SCAD (picasso) using a fixed list of lambdas
    ## Input: X, Y: data
    ##        lambda_list: list of lambda values
    ##        nfold: number of folds for cross validation
    ## Output: list SCADcv, containing: lambda.min, coef_min
    n = dim(X)[1]
    p = dim(X)[2]
  
    cv_score = rep(0, length(lambda_list))
    n0 = floor(n / nfold) * nfold
    ind_sample = matrix(sample(n0), nrow = nfold, ncol = n0 / nfold)
    for (k in 1:length(lambda_list)){
        lambda = lambda_list[k]
        for (j in 1:nfold){
            ind_test = ind_sample[j, ]
            ind_train = setdiff(c(1:n0), ind_test)
            X_test = X[ind_test, ]
            Y_test = Y[ind_test, ]
            X_train = X[ind_train, ]
            Y_train = Y[ind_train, ]
            model_train = picasso(X_train, Y_train, lambda=lambda, 
                                  method="scad", intercept=FALSE)
            beta_train = model_train$beta
            cv_score[k] = cv_score[k] + sum((Y_test - X_test %*% beta_train)^2)/nfold 
        }
    }
    ind_cv = which(cv_score == min(cv_score))
    lambda.min = lambda_list[ind_cv]
    model_min = picasso(X, Y, lambda=lambda.min, standardize = FALSE,
                        method="scad", intercept=FALSE)
    coef_min = model_min$beta
    SCADcv = list(lambda.min = lambda.min, coef_min = coef_min)
    return(SCADcv)
}

iht_cv <- function(X, Y, s_list, g_list, beta0, nfold = 10, n_cv = 100, maxiter = 1e3, prec = 1e-7){
    ## cross validation for IHT
    ## Input: X, Y: data
    ##        s_list, g_list: list of s and g values
    ##        beta0: initialization
    ##        nfold: number of folds for cross validation
    ##        n_cv: number of cross validations performed
    ## Output: list IHTcv, containing: s.min, g.min, coef_min, cv_score
    n = dim(X)[1]
    p = dim(X)[2]
    cv_score = rep(0, length(s_list))
    n0 = floor(n / nfold) * nfold
    for (ind_cv in 1:n_cv){
    ind_sample = matrix(sample(n0), nrow = nfold, ncol = n0 / nfold)
    for (j in 1:nfold){
        ind_test = ind_sample[j, ]
        ind_train = setdiff(c(1:n0), ind_test)
        X_test = X[ind_test, ]
        Y_test = Y[ind_test, ]
        X_train = X[ind_train, ]
        Y_train = Y[ind_train, ]
        model_iht = iht_cor2(X_train, Y_train, s_list=s_list, g_list=g_list, beta0=beta0, maxiter = maxiter, prec = prec)
        beta_iht = model_iht$beta
        for (k in 1:length(s_list)){
            cv_score[k] = cv_score[k] + sum((Y_test - X_test %*% beta_iht[, k])^2)/nfold/n_cv
        }
    }
    }
    ind_cv = which(cv_score == min(cv_score))
    s_iht.min = s_list[ind_cv]
    g_iht.min = g_list[ind_cv]
    model_iht = iht_cor2(X, Y, s_list=c(s_iht.min), g_list=c(g_iht.min), beta0=beta0, maxiter = maxiter, prec = prec)
    coef_min = model_iht$beta
    IHTcv = list(s.min = s_iht.min, g.min=g_iht.min, coef_min = coef_min, cv_score = cv_score)
    return(IHTcv)
}





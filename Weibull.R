
################## Weibull ################## 

weibull <- function(W_u, W_p, time, delta, epsilon = 0.001, tol = 1e-05,
                    nIter=2e4)
{  
  # W_u: N by nonp, non-penalized covariate matrix associated with latency
  # W_p:  N by M, penalized covariate matrix associated with incidence
  # time: vector of length N, observed survival time
  # delta: vector of length N, censoring status (not censored = 0)
  # epsilon: incremental size
  # tol: difference between log-likelihood
  
  N = length(time)
  M = ncol(W_p) # number of penalized latency covariates
  W_p = cbind(W_p, -W_p) # N by 2M
  
  ######## initialization ########
  step = 1
  beta_p = rep(0, 2*M) # penalized latency  
  alpha = uniroot(function(a) # MOM estimates
    log(gamma(1+2/a))-2*log(gamma(1+1/a))-log(var(time)+(mean(time))^2)+2*log(mean(time)),
    c(0.01,10))$root
  log_alpha = log(max(alpha, 1e-15))
  lambda = gamma(1+1/alpha)/mean(time)
  log_lambda = log(max(lambda, 1e-15))
  beta_u = rep(0, ncol(W_u))
  init = optim(par = c(log_alpha, log_lambda, beta_u), fn = negloglik_weibull, gr = gradient_weibull, 
               beta_p = beta_p, W_u = W_u, W_p = W_p, 
               time = time, delta = delta, method="BFGS")
  log_alpha = init$par[1]
  alpha = exp(log_alpha)
  log_lambda = init$par[2]
  lambda = exp(log_lambda)
  beta_u = init$par[3:(ncol(W_u)+2)] # unpenalized latency 
  LL0 = -init$value
  
  beta_p_path <- alpha_path <- lambda_path <- beta_u_path <- NULL
  logLikelihood <- numeric()
  ######## loop ########  
  
  repeat{
    
    #### update penalized parameters
    upd = update_weibull(alpha, lambda, beta_p, beta_u, W_u, W_p, time, delta, epsilon)
    beta_p = upd$beta_p
    beta_p_path = rbind(beta_p_path, beta_p)
    
    #### update other parameters
    out = optim(par = c(log_alpha, log_lambda, beta_u), fn = negloglik_weibull, gr = gradient_weibull, 
                beta_p = beta_p, W_u = W_u, W_p = W_p, 
                time = time, delta = delta, method="BFGS")
    log_alpha = out$par[1]
    alpha = exp(log_alpha)
    log_lambda = out$par[2]
    lambda = exp(log_lambda)
    beta_u = out$par[3:(ncol(W_u)+2)] # unpenalized latency  
    
    alpha_path = c(alpha_path, alpha)
    lambda_path = c(lambda_path, lambda)
    beta_u_path = rbind(beta_u_path, beta_u)
    
    LL1 = -out$value
    logLikelihood = c(logLikelihood, LL1)
    
    # cat("step=", step, "\n")
    if (step > 1 && (abs(LL1 - LL0) < tol | step >= nIter)) { # 
      break
    }
    LL0 <- LL1
    step <- 1 + step
  }
  
  ######## output ########
  beta_p_path = beta_p_path[,1:M] - beta_p_path[,(M+1):(2*M)]
  output <- list(beta_p_path = beta_p_path, alpha_path = alpha_path,
                 lambda_path = lambda_path, 
                 beta_u_path = beta_u_path, logLikelihood = logLikelihood)
  output
}

negloglik_weibull <- function(theta, beta_p, W_u, W_p, time, delta)
{
  N = length(time)
  P = ncol(W_u) + ncol(W_p)/2 # total number of latency covariates
  M = ncol(W_p)/2 # number of penalized latency covariates
  alpha = exp(theta[1])
  lambda = exp(theta[2])
  beta_u = theta[3:(P-M+2)] # unpenalized latency
  beta_nonzero = which(beta_p!=0)
  logC_beta = W_u %*% beta_u + W_p[,beta_nonzero,drop=FALSE] %*% beta_p[beta_nonzero]
  ll1 = delta * (log(max(lambda,1e-100)) + log(max(alpha,1e-100)) + 
                   (alpha-1)*log(pmax(lambda * time,rep(1e-100,N))) + 
                   logC_beta)
  ll2 = -(lambda * time)^alpha * exp(logC_beta)
  return(-sum(ll1+ll2))
}

gradient_weibull <- function(theta, beta_p, W_u, W_p, time, delta)
{
  N = length(time)
  P = ncol(W_u) + ncol(W_p)/2 # total number of latency covariates
  M = ncol(W_p)/2 # number of penalized latency covariates
  alpha = exp(theta[1])
  lambda = exp(theta[2])
  beta_u = theta[3:(P-M+2)] # unpenalized latency
  beta_nonzero = which(beta_p!=0)
  C_beta = exp(W_u %*% beta_u + W_p[,beta_nonzero,drop=FALSE] %*% beta_p[beta_nonzero])
  temp1 = log (pmax(lambda*time,rep(1e-100,N)))
  temp2 = (lambda * time)^alpha * C_beta
  grad1 = sum(delta * (1 + temp1*alpha) - temp2 * temp1 * alpha)
  grad2 = alpha * sum(delta - temp2)
  grad3 = matrix(delta - temp2, 1) %*% W_u
  return(-c(grad1, grad2, grad3))
}

update_weibull <- function(alpha, lambda, beta_p, beta_u, W_u, W_p, time, delta, epsilon)
{
  beta_nonzero = which(beta_p!=0)
  C_beta = exp(W_u %*% beta_u + W_p[,beta_nonzero,drop=FALSE] %*% beta_p[beta_nonzero])
  temp2 = (lambda * time)^alpha * C_beta
  ### update beta_p
  grad_beta = matrix(delta - temp2,1) %*% W_p # length M
  j_beta = which.max(grad_beta)
  beta_p[j_beta] = beta_p[j_beta] + epsilon
  
  return(list(beta_p = beta_p))
}

Weibull_fit <- function(W_u, W_p, time, delta, nIter=1e4, n_folds=5){
  # Cross-validation
  folds_i = sample(rep(1:n_folds, length.out = length(time)))
  Cstat = matrix(NA, nIter, n_folds)
  for (k in 1:n_folds) {
    test_i = which(folds_i == k)
    train_out = weibull(scale(W_u[-test_i,]), scale(W_p[-test_i,]), time[-test_i], delta[-test_i], nIter=nIter)
    test_delta = delta[test_i]
    test_time = time[test_i]
    test_W_u = scale(W_u[test_i,])
    test_W_p = scale(W_p[test_i,])
    Cstat[1:length(train_out$logLikelihood),k] = 
      sapply(1:length(train_out$logLikelihood), 
             function(x) 
               c_stat_beta(train_out$beta_u_path[x,], train_out$beta_p_path[x,],
                           test_delta, test_time, test_W_u, test_W_p))
  }
  c_stat = rowMeans(Cstat, na.rm = T)
  model_select = which.max(c_stat)
  output = weibull(scale(W_u), scale(W_p), time, delta, nIter=c_stat)
  return(output)
}
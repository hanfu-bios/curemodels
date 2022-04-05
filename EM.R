
################## EM ################## 
## Version 0319: change initialization

library(survival)
library(glmnet)

EM <- function(X_u, X_p, W_u, W_p, time, delta, mu, nIter, tol = 1e-4){
  # mu: penalty parameter 
  # tol: difference between log-likelihood
  N = length(time)
  J = ncol(X_p) # number of penalized incidence covariates
  M = ncol(W_p) # number of penalized latency covariates
  
  ######## initialization ########
  step = 1
  b_p <- beta_p <- rep(0,J)
  b_p_ext <- beta_p_ext <- rep(0,2*J)
  fit2 = coxph(as.formula(paste("Surv(time, delta) ~", paste(colnames(as.data.frame(W_u)),collapse=" + "))), 
               data =as.data.frame(W_u),ties = "breslow")
  b_u <- beta_u <- coef(fit2)
  # beta_u <- coef(fit2)
  # fit1 = glm(as.formula(paste("delta ~", paste(colnames(as.data.frame(X_u)),collapse=" + "))), 
  #            data =as.data.frame(X_u))
  # b_u = fit1$coefficients[-1]
  # itct = fit1$coefficients[1]
  
  alpha = uniroot(function(a) # MOM estimates
    log(gamma(1+2/a))-2*log(gamma(1+1/a))-log(var(time)+(mean(time))^2)+2*log(mean(time)),
    c(0.01,10))$root
  log_alpha = log(max(alpha, 1e-15))
  lambda = gamma(1+1/alpha)/mean(time)
  log_lambda = log(max(lambda, 1e-15))
  itct = 0
  pir = rep(1, N)
  pir[delta==0] = 1/(1+exp(-(itct+X_u[delta==0,,drop=FALSE] %*% b_u + X_p[delta==0,b_p!=0,drop=FALSE] %*% b_p[b_p!=0])+
                             (lambda * time[delta==0])^alpha*
                             exp(W_u[delta==0,,drop=FALSE] %*% beta_u + W_p[delta==0,beta_p!=0,drop=FALSE] %*% beta_p[beta_p!=0])))
  init_inc = optim(par = c(itct,b_u), fn = negloglik_inc, gr = gradient_inc,  
                   b_p = b_p, X_u =X_u, X_p = X_p, pir = pir, mu = mu, method="BFGS") # incidence
  itct = init_inc$par[1]
  b_u = init_inc$par[2:(ncol(X_u)+1)]
  llp1_0 <- llp1_init <- -init_inc$value
  init_lat = optim(par = c(log_alpha, log_lambda, beta_u), fn = negloglik_lat, gr = gradient_lat, 
                   beta_p=beta_p, W_u=W_u, W_p=W_p, time=time, delta=delta, pir=pir, mu=mu,
                   method="BFGS") # latency
  log_alpha = init_lat$par[1]
  alpha = exp(log_alpha)
  log_lambda = init_lat$par[2]
  lambda = exp(log_lambda)
  beta_u = init_lat$par[3:(ncol(W_u)+2)] # unpenalized incidence 
  llp2_0 <- llp2_init <- -init_lat$value
  
  llp1 <- llp2 <- numeric()
  conv1 <- conv2 <- FALSE
  
  ######## loop ########  
  
  repeat{
    
    #### E-step
    pir[delta==0] = 1/(1+exp(-(itct+X_u[delta==0,,drop=FALSE] %*% b_u + X_p[delta==0,b_p!=0,drop=FALSE] %*% b_p[b_p!=0])+
                               (lambda * time[delta==0])^alpha*
                               exp(W_u[delta==0,,drop=FALSE] %*% beta_u + W_p[delta==0,beta_p!=0,drop=FALSE] %*% beta_p[beta_p!=0])))
    
    #### update penalized parameters
    if (!conv1){
      b_p_ext = optim(par = b_p_ext, fn = negloglik_inc_b, gr = grad_b, itct = itct, b_u = b_u, X_u = X_u, X_p = X_p, 
                      pir = pir, mu = mu, method = 'L-BFGS-B', lower = rep(0, 2*J))$par
      b_p = b_p_ext[1:J]-b_p_ext[(J+1):(2*J)]
    }
    
    if (!conv2){
      beta_p_ext = optim(par = beta_p_ext, fn = negloglik_lat_beta, gr = grad_beta, alpha = alpha, lambda = lambda, 
                         beta_u = beta_u, W_u = W_u, W_p = W_p, time = time, delta = delta, pir = pir, mu = mu,
                         method = 'L-BFGS-B', lower = rep(0, 2*M))$par
      beta_p = beta_p_ext[1:M]-beta_p_ext[(M+1):(2*M)]
    }
    
    #### update nonpenalized parameters
    if (!conv1){
      out_inc = optim(par = c(itct,b_u), fn = negloglik_inc, gr = gradient_inc,  
                      b_p = b_p, X_u =X_u, X_p = X_p, pir = pir, mu = mu, method="BFGS") # incidence
      itct = out_inc$par[1]
      b_u = out_inc$par[2:(ncol(X_u)+1)]
      llp1_1 = -out_inc$value
    }
    llp1 = c(llp1, llp1_1)
    
    if (!conv2){
      out_lat = optim(par = c(log_alpha, log_lambda, beta_u), fn = negloglik_lat, gr = gradient_lat, 
                      beta_p=beta_p, W_u=W_u, W_p=W_p, time=time, delta=delta, pir=pir, mu=mu,
                      method="BFGS") # latency
      log_alpha = out_lat$par[1]
      alpha = exp(log_alpha)
      log_lambda = out_lat$par[2]
      lambda = exp(log_lambda)
      beta_u = out_lat$par[3:(ncol(W_u)+2)] # unpenalized incidence
      llp2_1 = -out_lat$value
    }
    llp2 = c(llp2, llp2_1)
    
    # cat("step=", step, "\n")
    if (!conv1 & abs(llp1_1- llp1_0)< tol) conv1 <- TRUE
    if (!conv2 & abs(llp2_1- llp2_0)< tol) conv2 <- TRUE
    if (step > 1 & (conv1 & conv2) | step >= nIter) { # 
      break
    }
    llp1_0 <- llp1_1
    llp2_0 <- llp2_1
    step <- 1 + step
  }
  
  ######## output ########
  output <- list(b_p = b_p, beta_p = beta_p, alpha = alpha, lambda = lambda, b_u = b_u, itct = itct,
                 beta_u = beta_u, llp1 = llp1, llp2 = llp2)
  output
}


############ unpenalized parameters #################


negloglik_inc <- function(theta, b_p, X_u , X_p, pir, mu) # incidence
{
  itct = theta[1]
  b_u = theta[2:(ncol(X_u)+1)]
  N = nrow(X_p)
  b_nonzero = which(b_p!=0)
  bx = itct + X_u %*% b_u + X_p[,b_nonzero, drop=FALSE] %*% b_p[b_nonzero]
  ll = sum(pir*bx - log(1+exp(bx)))-N*mu*sum(abs(b_p))
  return(-ll)
}

gradient_inc <- function(theta, b_p, X_u , X_p, pir, mu) # incidence
{
  itct = theta[1]
  b_u = theta[2:(ncol(X_u)+1)]
  b_nonzero = which(b_p!=0)
  bx = itct + X_u %*% b_u + X_p[,b_nonzero, drop=FALSE] %*% b_p[b_nonzero]
  pix = 1/(1+exp(-bx))
  grad1 = sum(pir-pix)
  grad2 = matrix(pir-pix,1) %*% X_u 
  return(-c(grad1, grad2))
}

negloglik_lat <- function(theta, beta_p, W_u, W_p, time, delta, pir, mu){ # latency
  N = nrow(W_p)
  log_alpha = theta[1]
  log_lambda = theta[2]
  alpha = exp(log_alpha)
  lambda = exp(log_lambda)
  beta_u = theta[3:(ncol(W_u)+2)]
  beta_nonzero = which(beta_p!=0)
  betaw = W_u %*% beta_u + W_p[,beta_nonzero, drop=FALSE] %*% beta_p[beta_nonzero]
  ll1 = delta*( log_alpha+log_lambda +(alpha-1)*log(pmax(lambda * time,1e-100))+betaw )
  ll2 = -pir*(lambda*time)^alpha*exp(betaw)
  return(-sum(ll1+ll2)+N*mu*sum(abs(beta_p)))
}

gradient_lat <- function(theta, beta_p, W_u, W_p, time, delta, pir, mu){ # latency
  log_alpha = theta[1]
  log_lambda = theta[2]
  alpha = exp(log_alpha)
  lambda = exp(log_lambda)
  beta_u = theta[3:(ncol(W_u)+2)]
  beta_nonzero = which(beta_p!=0)
  betaw = W_u %*% beta_u + W_p[,beta_nonzero, drop=FALSE] %*% beta_p[beta_nonzero]
  temp1 = pir*(lambda*time)^alpha*exp(betaw)
  temp2 = delta-temp1
  grad1 = sum(delta+log(pmax(lambda * time,1e-100))*temp2*alpha)
  grad2 = sum(alpha*temp2)
  grad3 = matrix(temp2,1) %*% W_u
  return(-c(grad1, grad2, grad3))
}

############ penalized parameters #################

negloglik_inc_b <- function(theta, itct, b_u, X_u , X_p, pir, mu) # incidence
{
  N = nrow(X_p)
  J = ncol(X_p)
  b_p_ext = theta
  b_p = b_p_ext[1:J]-b_p_ext[(J+1):(2*J)]
  b_nonzero = which(b_p!=0)
  bx = itct + X_u %*% b_u + X_p[,b_nonzero, drop=FALSE] %*% b_p[b_nonzero]
  ll = sum(pir*bx - log(1+exp(bx)))-N*mu*sum(b_p_ext)
  return(-ll)
}

grad_b <- function(theta, itct, b_u, X_u , X_p, pir, mu){
  N = nrow(X_p)
  J = ncol(X_p)
  b_p_ext = theta
  b_p = b_p_ext[1:J]-b_p_ext[(J+1):(2*J)]
  b_nonzero = which(b_p!=0)
  bx = itct + X_u %*% b_u + X_p[,b_nonzero, drop=FALSE] %*% b_p[b_nonzero]
  pix = 1/(1+exp(-bx))
  grad = matrix(pir-pix,1) %*% X_p
  return(c(-grad+N*mu, grad+N*mu))
}

negloglik_lat_beta <- function(theta, alpha, lambda, beta_u, W_u, W_p, time, delta, pir, mu){ # latency
  N = nrow(W_p)
  M = ncol(W_p)
  beta_p_ext = theta
  beta_p = beta_p_ext[1:M]-beta_p_ext[(M+1):(2*M)]
  beta_nonzero = which(beta_p!=0)
  betaw = W_u %*% beta_u + W_p[,beta_nonzero, drop=FALSE] %*% beta_p[beta_nonzero]
  ll1 = delta*( log(alpha)+log(lambda) +(alpha-1)*log(pmax(lambda * time,1e-100))+betaw )
  ll2 = -pir*(lambda*time)^alpha*exp(betaw)
  return(-sum(ll1+ll2)+N*mu*sum(beta_p_ext))
}

grad_beta <- function(theta, alpha, lambda, beta_u, W_u, W_p, time, delta, pir, mu){
  N = nrow(W_p)
  M = ncol(W_p)
  beta_p_ext = theta
  beta_p = beta_p_ext[1:M]-beta_p_ext[(M+1):(2*M)]
  beta_nonzero = which(beta_p!=0)
  betaw = W_u %*% beta_u + W_p[,beta_nonzero, drop=FALSE] %*% beta_p[beta_nonzero]
  temp1 = pir*(lambda*time)^alpha*exp(betaw)
  grad = matrix(delta-temp1,1)%*% W_p
  return(c(-grad+ N*mu, grad+ N*mu))
}

em_fit <- function(X_u, X_p, W_u, W_p, time, delta, nIter=100, n_folds=5, n_mu=50){
  # Cross-validation
  mu_max = max(colSums(abs(scale(X_p))))/(2*length(time))
  mu_min = mu_max * 0.01
  k1 = (0:(n_mu-1)) / n_mu
  mu_list = mu_max * (mu_min/mu_max)^k1
  folds_i = sample(rep(1:n_folds, length.out = length(time)))
  Cstat = matrix(NA, n_mu, n_folds)
  for (k in 1:n_folds) { 
    test_i = which(folds_i == k)
    test_delta = delta[test_i]
    test_time = time[test_i]
    test_X_u = scale(X_u[test_i,])
    test_X_p = scale(X_p[test_i,])
    test_W_u = scale(W_u[test_i,])
    test_W_p = scale(W_p[test_i,])
    Cstat[,k] = sapply(mu_list,
                       function(mu) {
                         train_out = EM(scale(X_u[-test_i,]), scale(X_p[-test_i,]),
                                        scale(W_u[-test_i,]), scale(W_p[-test_i,]),
                                        time[-test_i], delta[-test_i], mu=mu, nIter, tol = 1e-3)
                         cs = C.stat(cure_cutoff=5,
                                     train_out$b_u, train_out$itct, train_out$b_p,
                                     train_out$beta_u, train_out$beta_p,
                                     test_delta, test_time, test_X_u, test_X_p, test_W_u, test_W_p)
                         return(cs)
                       })
    cat("Fold", k, "training finished\n")
  }
  c_stat = rowMeans(Cstat, na.rm = T)
  model_select = which.max(c_stat)
  cat("Selected mu:",round(mu_list[model_select],4),"\n")
  output = EM(scale(X_u), scale(X_p), scale(W_u), scale(W_p),
              time, delta, mu=mu_list[model_select], nIter=nIter)
  return(output)
}

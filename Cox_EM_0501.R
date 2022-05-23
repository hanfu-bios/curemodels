################## EM for Cox ################## 
###
## Version 05/01
## Changes: 1. grid tuning
##          2. allow X_u and W_u to be null, X_p and W_p to be different

library(glmnet)
library(survival)
library(smcure)
library(knockoff)

cox.cure <- function(X_u=NULL, X_p, W_u=NULL, W_p, time, delta, mu_inc, mu_lat, inits, nIter=100, tol = 1e-4)
{  
  # mu: penalty parameter 
  # tol: difference between log-likelihood
  N = length(time)
  J = ncol(X_p) # number of penalized incidence covariates
  M = ncol(W_p) # number of penalized latency covariates
  
  CAP = 10
  
  event_time = time[delta==1]
  uniq_event_time = unique(event_time)
  n0 = length(uniq_event_time)
  I0 = matrix(time, ncol = 1) %*% matrix(1, 1, n0) >= 
    matrix(1, N, 1) %*% matrix(uniq_event_time, nrow = 1)    # N*n0
  d_j = as.numeric(table(event_time)[rank(unique(event_time))])   #  number of events at time T_j
  T_n0 = max(event_time)   # last event
  tail_ind = which(time>=T_n0 & delta==0)
  
  ######## initialization ########
  step = 1
  b_p <- b_p0 <- rep(0,J)
  beta_p <- beta_p0 <- rep(0,M)
  itct = inits$itct; b_u = inits$b_u; beta_u = inits$beta_u; survprob = inits$survprob
  pir = rep(1, N)
  if(is.null(X_u)) b_x = rep(itct,N) else b_x = itct+X_u %*% b_u
  lik_inc <- lik_lat <- lik_inc_p_list <- lik_lat_p_list <- c()
  b_p_path <- beta_p_path <- b_u_path <- beta_u_path <- itct_path <- NULL
  if(is.null(X_u)) pen_fac_inc = rep(1, ncol(X_p)) 
  else pen_fac_inc = c(rep(0, ncol(X_u)), rep(1, ncol(X_p)))
  if(is.null(W_u)) pen_fac_lat = rep(1, ncol(W_p)) 
  else pen_fac_lat = c(rep(0, ncol(W_u)), rep(1, ncol(W_p)))
  lik_p0 = 0
  diff_b <- diff_beta <- c()
  ######## loop ########  
  
  repeat{
    
    #### E-step
    pi_x = 1/(1+exp(-b_x))
    numerator = pi_x[delta==0] * survprob[delta==0]
    pir[delta==0] = numerator/(1-pi_x[delta==0] + numerator)
    
    #### M-step
    ### incidence
    fit_inc = glmnet(x = cbind(X_u, X_p), y = cbind(1-pir, pir), family = "binomial", 
                     penalty.factor = pen_fac_inc, lambda=mu_inc)
    coef_inc = coef(fit_inc)
    itct = max(-CAP, min(CAP, coef_inc[1]))
    if(is.null(X_u)) {
      b_p = pmax(-CAP, pmin(CAP, coef_inc[-1]))
      b_x = itct + X_p[,b_p!=0,drop=FALSE] %*% b_p[b_p!=0]
    } else{
      b_u = pmax(-CAP, pmin(CAP, coef_inc[2:(ncol(X_u)+1)]))
      b_p = pmax(-CAP, pmin(CAP, coef_inc[-(1:(ncol(X_u)+1))]))
      b_x = itct + X_u %*% b_u + X_p[,b_p!=0,drop=FALSE] %*% b_p[b_p!=0]
    }
    lik_inc = c(lik_inc,
                sum(pir*b_x - log(1+exp(b_x))) )#- N * mu *sum(abs(b_p)))
    lik_inc_p = sum(pir*b_x - log(1+exp(b_x))) - N * mu_inc *sum(abs(b_p))
    lik_inc_p_list = c(lik_inc_p_list, lik_inc_p)
    ### latency
    nz_pir = which(pir>0)
    fit_lat = glmnet(x = cbind(W_u, W_p)[nz_pir,], 
                     y = Surv(time = time[nz_pir], event = delta[nz_pir]),
                     family = "cox",
                     offset = log(pir[nz_pir]),
                     penalty.factor = pen_fac_lat,
                     lambda = mu_lat)
    coef_lat = coef(fit_lat)
    if(is.null(W_u)) {
      beta_p = pmax(-CAP, pmin(CAP, as.numeric(coef_lat)))
      beta_w = W_p[,beta_p!=0,drop=FALSE] %*% beta_p[beta_p!=0]
    }else{
      beta_u = pmax(-CAP, pmin(CAP, coef_lat[1:ncol(W_u)]))
      beta_p = pmax(-CAP, pmin(CAP, coef_lat[-(1:ncol(W_u))]))
      beta_w = W_u %*% beta_u + W_p[,beta_p!=0,drop=FALSE] %*% beta_p[beta_p!=0]
    }
    exp_beta_w_p = exp(beta_w) * pir
    denom = matrix(exp_beta_w_p,1) %*% I0   # 1*n0
    hazard = rep(0, N)
    hazard[delta==1] = sapply(event_time, function(x) (d_j/denom)[uniq_event_time==x])
    accum_hazard = sapply(time, function(x) sum((d_j/denom)[uniq_event_time<=x]))
    surv_baseline = exp(-accum_hazard)
    if (length(tail_ind)>0){
      wtail_alpha = uniroot(function(a)
        sum(-exp_beta_w_p*max(accum_hazard)*log(time/T_n0)*(time/T_n0)^a + delta/a + delta*log(time/T_n0)),
        c(0.01,10), extendInt = "yes")$root
      wtail_lambda = (max(accum_hazard))^(1/wtail_alpha)/T_n0
      surv_baseline[tail_ind] = exp(-(wtail_lambda*time[tail_ind])^wtail_alpha)
    }
    survprob = surv_baseline^exp(beta_w)
    lik_lat = c(lik_lat, sum(beta_w[delta==1]) - sum(d_j*log(denom)))
    lik_lat_p = sum(log(hazard[delta==1]) + beta_w[delta==1]) - sum(exp_beta_w_p*accum_hazard) -
      N * mu_lat *sum(abs(beta_p))
    lik_lat_p_list = c(lik_lat_p_list, lik_lat_p)
    
    diff_b = c(diff_b, max(abs(b_p - b_p0)))
    diff_beta = c(diff_beta, max(abs(beta_p - beta_p0)))
    itct_path = c(itct_path, itct)
    if(!is.null(X_u)) b_u_path = rbind(b_u_path, b_u)
    b_p_path = rbind(b_p_path, b_p)
    if(!is.null(W_u)) beta_u_path = rbind(beta_u_path, beta_u)
    beta_p_path = rbind(beta_p_path, beta_p)
    # cat("step=", step, "\n")
    if (is.na(lik_inc_p+lik_lat_p)){
      cat("WARNING: Repetition failed at mu_b =",mu_inc, "and mu_beta =", mu_lat, "\n")
      cat("lik_inc_p=", lik_inc_p, "\n")
      cat("lik_lat_p=", lik_lat_p, "\n")
      break
    }
    if (step >= nIter | 
        (step > 1 & abs(lik_inc_p+lik_lat_p - lik_p0) < tol)) {
      # (step > 1 & max(abs(b_p - b_p0)) < tol & max(abs(beta_p - beta_p0)) < tol)) { #
      break
    }
    step <- 1 + step
    b_p0 = b_p
    beta_p0 = beta_p
    lik_p0 = lik_inc_p+lik_lat_p
  }
  
  ######## output ########
  output <- list(b_p_path = b_p_path, beta_p_path = beta_p_path, b_u_path = b_u_path, itct_path = itct_path,
                 beta_u_path = beta_u_path, 
                 lik_inc = lik_inc, lik_lat = lik_lat, 
                 lik_inc_p_list = lik_inc_p_list, lik_lat_p_list = lik_lat_p_list, 
                 diff_b=diff_b, diff_beta=diff_beta,
                 surv_baseline = surv_baseline)
  output
}


initialization <- function(X_u=NULL, W_u=NULL, time, delta){
  N = length(time)
  pir = rep(1, N)
  data_init = as.data.frame(cbind(time, delta, X_u, W_u)) 
  if(is.null(X_u) & is.null(W_u)){
    fit = glm(delta~1, family="binomial")
    pi_x = mean(delta)
    surv = survfit(Surv(time, delta)~1)
    survprob = sapply(time, function(x) summary(surv, times = x, extend=T)$surv)
    numerator = pi_x * survprob[delta==0]
    pir[delta==0] = numerator/(1-pi_x + numerator)
    inits = list(itct=coef(fit), b_u=NULL, beta_u=NULL, survprob=survprob, pir=pir)
  } else if(is.null(X_u)){
    colnames(data_init) <- c("time","delta",paste0("lat",1:ncol(W_u)))
    formula_lat = as.formula(paste("Surv(time, delta) ~", paste(paste0("lat",1:ncol(W_u)),collapse=" + ")))
    fit1 = coxph(formula_lat, data = data_init)
    surv = survfit(fit1)
    survprob = sapply(time, function(x) summary(surv, times = x, extend=T)$surv)
    fit2 = glm(delta~1, family="binomial")
    pi_x = mean(delta)
    numerator = pi_x * survprob[delta==0]
    pir[delta==0] = numerator/(1-pi_x + numerator)
    inits = list(itct=coef(fit2), b_u=NULL, beta_u=coef(fit1), survprob=survprob, pir=pir)
  } else if(is.null(W_u)){
    colnames(data_init) <- c("time","delta",paste0("inc",1:ncol(X_u)))
    formula_inc = as.formula(paste("delta ~", paste(paste0("inc",1:ncol(X_u)),collapse=" + ")))
    fit = glm(formula_inc, family="binomial", data = data_init)
    surv = survfit(Surv(time, delta)~1)
    survprob = sapply(time, function(x) summary(surv, times = x, extend=T)$surv)
    itct = coef(fit)[1]
    b_u = coef(fit)[-1]
    b_x = itct+X_u %*% b_u
    pi_x = 1/(1+exp(-b_x))
    numerator = pi_x[delta==0] * survprob[delta==0]
    pir[delta==0] = numerator/(1-pi_x[delta==0] + numerator)
    inits = list(itct=itct, b_u=b_u, beta_u=NULL, survprob=survprob, pir=pir)
  } else{
    colnames(data_init) <- c("time","delta",paste0("inc",1:ncol(X_u)), paste0("lat",1:ncol(W_u)))
    formula_lat = as.formula(paste("Surv(time, delta) ~", paste(paste0("lat",1:ncol(W_u)),collapse=" + ")))
    formula_inc = as.formula(paste("delta ~", paste(paste0("inc",1:ncol(X_u)),collapse=" + ")))
    # invisible(capture.output(fit <- smcure(formula_lat, formula_inc, data=data_init, model = "ph")))
    # itct = fit$b[1]
    # b_u = fit$b[2:length(fit$b)]
    # b_x = itct+X_u %*% b_u
    # beta_u = fit$beta
    # surv_baseline = fit$s
    # beta_w = W_u %*% beta_u
    # survprob = surv_baseline^exp(beta_w)
    fit1 = glm(formula_inc, family="binomial", data = data_init)
    itct = coef(fit1)[1]
    b_u = coef(fit1)[-1]
    b_x = itct+X_u %*% b_u
    fit2 = coxph(formula_lat, data = data_init)
    beta_u=coef(fit2)
    surv = survfit(fit2)
    survprob = sapply(time, function(x) summary(surv, times = x, extend=T)$surv)
    
    pi_x = 1/(1+exp(-b_x))
    numerator = pi_x[delta==0] * survprob[delta==0]
    pir[delta==0] = numerator/(1-pi_x[delta==0] + numerator)
    inits = list(itct=itct, b_u=b_u, beta_u=beta_u, survprob=survprob, pir=pir)
  }
  return(inits)
}

self_scale <- function(X){
  if(is.null(X)) return(NULL)
  n = nrow(X)
  Xs = apply(X, 2, function(x) if(sd(x)==0) return(rep(0,n)) else return(scale(x)))
  return(Xs)
} 

cox_em_fit <- function(X_u=NULL, X_p, W_u=NULL, W_p, time, delta, cure_cutoff=5,
                       nIter=100, n_folds=5, n_mu=50){
  X_u = self_scale(X_u)
  X_p = self_scale(X_p)
  W_u = self_scale(W_u)
  W_p = self_scale(W_p)
  
  inits = initialization(X_u, W_u, time, delta)
  N = length(time)
  pir = inits$pir
  nz_pir = which(pir>0)
  mu_b_max = max(abs(matrix(pir-0.5,1,N)%*%X_p))/N
  mu_b_min = mu_b_max * 0.1
  k1 = (0:(n_mu-1)) / n_mu
  mu_b_list = mu_b_max * (mu_b_min/mu_b_max)^k1
  mu_beta_max = glmnet:::get_cox_lambda_max(x = cbind(W_u, W_p)[nz_pir,], 
                                            y = Surv(time = time[nz_pir], event = delta[nz_pir]),
                                            offset = log(pir[nz_pir]), alpha = 1)
  # mu_beta_max = mu_beta_max * 1.2
  mu_beta_min = mu_beta_max * 0.1
  mu_beta_list = mu_beta_max * (mu_beta_min/mu_beta_max)^k1
  mu_grid = expand.grid(mu_b_list, mu_beta_list)
  
  # Cross-validation
  
  folds_i = sample(rep(1:n_folds, length.out = length(time)))
  Cstat <- AUC <- matrix(NA, nrow(mu_grid), n_folds)
  
  for (k in 1:n_folds) { 
    test_i = which(folds_i == k)
    test_delta = delta[test_i]
    test_time = time[test_i]
    test_X_p = X_p[test_i,]
    test_W_p = W_p[test_i,]
    if(is.null(X_u)) test_X_u = NULL else test_X_u = X_u[test_i,]
    if(is.null(W_u)) test_W_u = NULL else test_W_u = W_u[test_i,]
    if(is.null(X_u)) X_u_train = NULL else X_u_train = X_u[-test_i,]
    if(is.null(W_u)) W_u_train = NULL else W_u_train = W_u[-test_i,]
    initial_val = inits
    initial_val$survprob = initial_val$survprob[-test_i]
    for (i in 1:nrow(mu_grid)){
      mu_b = mu_grid[i,1]
      mu_beta = mu_grid[i,2]
      train_out = cox.cure(X_u_train, X_p[-test_i,],
                           W_u_train, W_p[-test_i,],
                           time[-test_i], delta[-test_i], mu_inc=mu_b, mu_lat=mu_beta,
                           inits=initial_val, nIter, tol = 1e-3)
      model_select = which.max(train_out$lik_inc_p_list+train_out$lik_lat_p_list)
      itct = train_out$itct_path[model_select]
      if(is.null(X_u)) b_u = NULL else b_u = train_out$b_u_path[model_select,]
      b_p = train_out$b_p_path[model_select,]
      if(is.null(W_u)) beta_u = NULL else beta_u = train_out$beta_u_path[model_select,]
      beta_p = train_out$beta_p_path[model_select,]
      Cstat[i,k] = C.stat(cure_cutoff=cure_cutoff,
                  b_u, itct, b_p,
                  beta_u, beta_p,
                  test_delta, test_time, test_X_u, test_X_p, test_W_u, test_W_p)
      AUC[i,k] = AUC_msi(cure_cutoff = cure_cutoff, b_u, itct, b_p, 
                         test_delta, test_time, test_X_u, test_X_p)
      
    }
    cat("Fold", k, "training finished\n")
  }
  c_stat = rowMeans(Cstat, na.rm = T)
  auc_val = rowMeans(AUC, na.rm = T)
  model_select_b = which.max(auc_val)
  model_select_beta = which.max(c_stat)
  optimal_mu_b = mu_grid[model_select_b,1]
  optimal_mu_beta = mu_grid[model_select_beta,2]
  cat("selected mu_b:",round(optimal_mu_b,3),", mu_beta:", round(optimal_mu_beta,3),"\n")
  cat("selected index mu_b:",ifelse(model_select_b%%n_mu==0, n_mu, model_select_b%%n_mu),
      ", mu_beta:", (model_select_beta-1)%/%n_mu+1,"\n")
  cat("optimal c_stat:",max(c_stat),", AUC:", max(auc_val),"\n")
  
  output = cox.cure(X_u, X_p, W_u, W_p,
                    time, delta, mu_inc=optimal_mu_b, mu_lat=optimal_mu_beta, 
                    inits=inits,
                    nIter=nIter)
  return(output)
}

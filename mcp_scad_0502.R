################## Cox MCM (MCP/SCAD) ################## 
###
## Version 05/02
## Change: grid tuning, allow X_u and W_u to be null
###

library(smcure)
library(knockoff)

cox_mcp_scad <- function(X_u=NULL, X_p, W_u=NULL, W_p, time, delta, penalty = c("MCP", "SCAD"),
                         lambda_inc, lambda_lat, gamma_inc = 3, gamma_lat = 3, inits = NULL, 
                         nIter=100, tol = 1e-4)
{  
  # lambda: penalty parameter 
  # tol: difference between log-likelihood
  N = length(time)
  J = ncol(X_p) # number of penalized incidence covariates
  M = ncol(W_p) # number of penalized latency covariates
  
  CAP = 20
  
  event_time = time[delta==1]
  uniq_event_time = unique(event_time)
  n0 = length(uniq_event_time)
  I0 = matrix(time, ncol = 1) %*% matrix(1, 1, n0) >= 
    matrix(1, N, 1) %*% matrix(uniq_event_time, nrow = 1)    # N*n0
  d_j = as.numeric(table(event_time)[rank(uniq_event_time)])   #  number of events at time T_j
  Z = sapply(uniq_event_time, function(t) colSums(W_p[time==t & delta==1,,drop=FALSE])) # M*n0
  T_n0 = max(event_time)   # last event
  tail_ind = which(time>=T_n0 & delta==0)
  
  ######## initialization ########
  step = 1
  b_p <- rep(0,J)
  beta_p <- rep(0,M)
  itct = inits$itct; b_u = inits$b_u; beta_u = inits$beta_u; survprob = inits$survprob
  if(is.null(X_u)) b_x = rep(itct,N) else b_x = itct+X_u %*% b_u
  if(is.null(W_u)) beta_w = rep(0,N) else beta_w = W_u %*% beta_u
  pir = rep(1, N)
  exp_beta_w_p = exp(beta_w) * pir
  lik_inc <- lik_lat <- c()
  b_p_path <- beta_p_path <- b_u_path <- beta_u_path <- itct_path <- NULL
  llp1_0 <- llp2_0 <- 0
  conv1 <- conv2 <- FALSE
  ######## loop ########  
  
  repeat{
    
    #### E-step
    pi_x = 1/(1+exp( pmin(CAP, -b_x) ))
    numerator = pi_x[delta==0] * survprob[delta==0]
    pir[delta==0] = numerator/pmax(1e-50,1-pi_x[delta==0] + numerator)
    
    #### M-step
    ### incidence
    if (!conv1){
      for (l in 1:J){
        bl = b_p[l]
        b_p[l] = max(-CAP, min(CAP, cd_upd_inc(penalty=penalty, pir, pi_x, X_p[,l], bl, gamma_inc,
                                               lambda_inc)))
        b_x = b_x + X_p[,l,drop=FALSE] %*% (b_p[l] - bl)
        pi_x = 1/(1+exp( pmin(CAP, -b_x) ))
      }
      out_inc = optim(par = c(itct,b_u), fn = negloglik_inc, gr = gradient_inc, 
                      b_p = b_p, X_u =X_u, X_p = X_p, pir = pir, CAP=CAP, method="BFGS") # incidence
      itct = max(-CAP, min(CAP, out_inc$par[1]))
      if(!is.null(X_u)) b_u = pmax(-CAP, pmin(CAP, out_inc$par[2:(ncol(X_u)+1)]))
      pen = ifelse(penalty=="MCP", mcp_penalty(b_p, gamma_inc, lambda_inc), 
                   scad_penalty(b_p, gamma_inc, lambda_inc))
      llp1_1 = -out_inc$value - N*pen
      if(is.null(X_u)) b_x = itct + X_p[,b_p!=0,drop=FALSE] %*% b_p[b_p!=0]
      else b_x = itct + X_u %*% b_u + X_p[,b_p!=0,drop=FALSE] %*% b_p[b_p!=0]
    }
    lik_inc = c(lik_inc, llp1_1)
    
    ### latency
    if (!conv2){
      for (l in 1:M){
        betal = beta_p[l]
        beta_p[l] = max(-CAP, min(CAP, cd_upd_lat(penalty, exp_beta_w_p, W_p[,l], betal, 
                                                  Z[l,], I0, d_j, gamma_lat, lambda_lat)))
        change = W_p[,l,drop=FALSE] %*% (beta_p[l] - betal)  # N*1
        beta_w = beta_w + change  # N*1
        exp_beta_w_p = exp( pmin(CAP, beta_w) ) * pir  # N*1
        l = l + 1
      }
      if(!is.null(W_u)){
        out_lat = optim(par = beta_u, fn = negloglik_lat, gr = gradient_lat, 
                        beta_p=beta_p, W_u=W_u, W_p=W_p, delta=delta, pir=pir, I0=I0, d_j=d_j,
                        CAP=CAP, method="BFGS") # latency
        beta_u = pmax(-CAP, pmin(CAP, out_lat$par)) # unpenalized incidence
        beta_w = W_u %*% beta_u + W_p[,beta_p!=0,drop=FALSE] %*% beta_p[beta_p!=0]
        exp_beta_w_p = exp( pmin(CAP, beta_w) ) * pir
      }
      
      ### survival
      denom = matrix(exp_beta_w_p,1) %*% I0   # 1*n0
      hazard = rep(0, N)
      hazard[delta==1] = sapply(event_time, function(x) (d_j/denom)[uniq_event_time==x])
      accum_hazard = sapply(time, function(x) sum((d_j/denom)[uniq_event_time<=x]))
      surv_baseline = exp(-accum_hazard)
      if (length(tail_ind)>0){
        tmp <- try(wtail_alpha <- uniroot(function(a)
          sum(-exp_beta_w_p*max(accum_hazard)*log(time/T_n0)*(time/T_n0)^a + delta/a + delta*log(time/T_n0)),
          c(0.01,10), extendInt = "yes")$root, silent = TRUE)
        if (!inherits(tmp, "try-error")){
          wtail_lambda = (max(accum_hazard))^(1/wtail_alpha)/T_n0
          surv_baseline[tail_ind] = exp(-(wtail_lambda*time[tail_ind])^wtail_alpha)
        }
      }
      survprob = surv_baseline^exp(pmin(CAP, beta_w))
      pen = ifelse(penalty=="MCP", mcp_penalty(beta_p, gamma_lat, lambda_lat), 
                   scad_penalty(beta_p, gamma_lat, lambda_lat))
      llp2_1 = sum(log(hazard[delta==1]) + beta_w[delta==1]) 
      - sum(exp_beta_w_p*accum_hazard) - N*pen
    }
    lik_lat = c(lik_lat, llp2_1)
    
    ## record updated parameters
    itct_path = c(itct_path, itct)
    if(!is.null(X_u)) b_u_path = rbind(b_u_path, b_u)
    b_p_path = rbind(b_p_path, b_p)
    if(!is.null(W_u)) beta_u_path = rbind(beta_u_path, beta_u)
    beta_p_path = rbind(beta_p_path, beta_p)
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
  output <- list(b_p_path = b_p_path, b_u_path = b_u_path, itct_path = itct_path,
                 beta_p_path = beta_p_path, beta_u_path = beta_u_path, 
                 lik_inc = lik_inc, lik_lat = lik_lat)
  output
}


######## helper functions ########

soft <- function(z, gamma){
  if(gamma >= abs(z)) return(0)
  else return(ifelse(z>0, z-gamma, z+gamma))
}

mcp_penalty <- function(b_p, gamma, lambda){
  idx = which(b_p !=0 & abs(b_p) <= gamma*lambda)
  s1 = lambda*sum(abs(b_p[idx]))-sum((b_p[idx])^2)/(2*gamma)
  s2 = gamma*lambda^2/2 * sum(abs(b_p) > gamma*lambda)
  return(s1 + s2)
}

scad_penalty <- function(b_p, gamma, lambda){
  idx1 = which(b_p !=0 & abs(b_p) <= lambda)
  idx2 = which(abs(b_p) > lambda & abs(b_p)<=gamma*lambda)
  s1 = lambda*sum(abs(b_p[idx1]))
  s2 = (gamma*lambda*sum(abs(b_p[idx2])) - 0.5*(sum((b_p[idx2])^2 + lambda^2)))/(gamma-1)
  s3 = lambda^2*(gamma^2-1)/(2*(gamma-1)) * sum(abs(b_p) > gamma*lambda)
  return(s1 + s2 + s3)
}


######## incidence functions ########

cd_upd_inc <- function(penalty, pir, pi_x, x_l, bl, gamma, lambda){
  d1 = mean((pir - pi_x)*x_l)
  vl = mean(pi_x*(1-pi_x)*x_l^2)
  zl = d1 + vl * bl
  if(penalty == "MCP"){
    if(abs(bl) <= gamma*lambda) return(soft(zl, lambda)/(vl-1/gamma))
    else return(ifelse(vl==0,0,zl/vl))
  }  else if(penalty == "SCAD"){
    if(abs(bl)<=lambda) return(ifelse(vl==0,0,soft(zl, lambda)/vl))
    else if(abs(bl)<=gamma*lambda) return(soft(zl, gamma*lambda/(gamma-1))/(vl-1/(gamma-1)))
    else return(ifelse(vl==0,0,zl/vl))
  }
}

negloglik_inc <- function(theta, b_p, X_u , X_p, pir, CAP) # incidence
{
  itct = theta[1]
  if(!is.null(X_u)) b_u = theta[2:(ncol(X_u)+1)]
  N = nrow(X_p)
  b_nonzero = which(b_p!=0)
  if(is.null(X_u)) bx = itct + X_p[,b_nonzero, drop=FALSE] %*% b_p[b_nonzero]
  else bx = itct + X_u %*% b_u + X_p[,b_nonzero, drop=FALSE] %*% b_p[b_nonzero]
  ll = sum(pir*bx - log(1+exp(pmin(CAP, bx))))
  return(-ll)
}

gradient_inc <- function(theta, b_p, X_u , X_p, pir, CAP) # incidence
{
  itct = theta[1]
  if(!is.null(X_u)) b_u = theta[2:(ncol(X_u)+1)]
  b_nonzero = which(b_p!=0)
  if(is.null(X_u)) bx = itct + X_p[,b_nonzero, drop=FALSE] %*% b_p[b_nonzero]
  else bx = itct + X_u %*% b_u + X_p[,b_nonzero, drop=FALSE] %*% b_p[b_nonzero]
  pix = 1/(1+exp( pmin(CAP, -bx) ))
  grad1 = sum(pir-pix)
  if(!is.null(X_u)) grad2 = matrix(pir-pix,1) %*% X_u else grad2 = NULL
  return(-c(grad1, grad2))
}


######## latency functions ########

cd_upd_lat <- function(penalty, exp_beta_w_p, w_l, betal, Zl, I0, d_j, gamma, lambda){
  N = length(exp_beta_w_p)
  sum0 = pmax(1e-50, matrix(exp_beta_w_p,1) %*% I0)   # 1*n0
  sum1 = matrix(exp_beta_w_p * w_l,1) %*% I0   # 1*n0
  sum2 = matrix(exp_beta_w_p * w_l^2,1) %*% I0   # 1*n0
  d1 = sum(Zl - d_j * sum1/sum0)/N
  vl = -sum(d_j * (sum1^2/sum0^2 - sum2/sum0))/N
  zl = d1 + vl * betal
  if(penalty == "MCP"){
    if(abs(betal) <= gamma*lambda) return(soft(zl, lambda)/(vl-1/gamma))
    else return(ifelse(vl==0,0,zl/vl))
  }  else if(penalty == "SCAD"){
    if(abs(betal)<=lambda) return(ifelse(vl==0,0,soft(zl, lambda)/vl))
    else if(abs(betal)<=gamma*lambda) return(soft(zl, gamma*lambda/(gamma-1))/(vl-1/(gamma-1)))
    else return(ifelse(vl==0,0,zl/vl))
  }
} 

negloglik_lat <- function(theta, beta_p, W_u, W_p, delta, pir, I0, d_j, CAP){ # latency
  beta_u = theta
  N = nrow(W_p)
  beta_nonzero = which(beta_p!=0)
  beta_w = W_u %*% beta_u + W_p[,beta_nonzero, drop=FALSE] %*% beta_p[beta_nonzero]
  exp_beta_w_p = exp( pmin(CAP, beta_w) ) * pir
  sum0 = pmax(1e-50, matrix(exp_beta_w_p,1) %*% I0)  # 1*n0
  ll = sum(beta_w[delta==1]) - sum(d_j * log(sum0))
  return(-ll)
}

gradient_lat <- function(theta, beta_p, W_u, W_p, delta, pir, I0, d_j, CAP){ # latency
  beta_u = theta
  N = nrow(W_p)
  beta_nonzero = which(beta_p!=0)
  beta_w = W_u %*% beta_u + W_p[,beta_nonzero, drop=FALSE] %*% beta_p[beta_nonzero]
  exp_beta_w_p = exp( pmin(CAP, beta_w) ) * pir
  sum0 = pmax(1e-50, matrix(exp_beta_w_p,1) %*% I0)   # 1*n0
  sum1 = t(W_u) %*% diag(as.numeric(exp_beta_w_p)) %*% I0   # ncol(W_u) * n0
  du = colSums(W_u[delta==1,,drop=FALSE])- rowSums(sum1 %*% diag(as.numeric(d_j/sum0))) # ncol(W_u) * 1
  return(-du)
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

mcp_scad_fit <- function(X_u=NULL, X_p, W_u=NULL, W_p, time, delta, penalty, gamma_inc, gamma_lat, 
                         cure_cutoff=5, nIter=100, n_folds=5, n_mu=50){
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
      train_out = cox_mcp_scad(X_u_train, X_p[-test_i,],
                               W_u_train, W_p[-test_i,],
                               time[-test_i], delta[-test_i], penalty = penalty,
                               lambda_inc=mu_b, lambda_lat=mu_beta, 
                               gamma_inc = gamma_inc, gamma_lat = gamma_lat, 
                               inits=initial_val, nIter, tol = 1e-3)
      model_select = which.max(train_out$lik_inc+train_out$lik_lat)
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
  
  output = cox_mcp_scad(X_u, X_p, W_u, W_p,
                        time, delta, penalty = penalty, 
                        lambda_inc = optimal_mu_b, lambda_lat=optimal_mu_beta, 
                        gamma_inc = gamma_inc, gamma_lat = gamma_lat, inits=inits, nIter=nIter)
  return(output)
  
}

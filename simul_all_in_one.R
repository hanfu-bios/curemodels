rm(list = ls())

######## input
## A = 0.4, 0.6, 0.8, 1.0, ..., 1.8
## J = 100, 500
## seed = 1, 2, 3, 4, 5
## method_ind = 1, 2, 3, 4, 5, 6, 7 ("gmifs", "EM", "SCinCRM", "cmix", "Weibull", "CoxEM", "Cox")
## knockoff = TRUE, FALSE
## same_signs = TRUE, FALSE ("same", "random")
## model_ind = 1, 2, 3, 4, 5 
##     (1 for "Weibull", 2 for "GG", 3 for "Gompertz", 4 for "nonparametric", 5 for "GG baseline")
## itct_mean = 0.5, 1.5 (cure rate 0.38, 0.18)
## cens_ub = 10, 20
## cure_cutoff = 5, 15

args=(commandArgs(TRUE))
for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
}

method = c("gmifs", "EM", "SCinCRM", "cmix", "Weibull", "Cox_EM", "glmnet_Cox")[method_ind]

library(doMC)
library(knockoff)
library(Rcpp)
library(glmnet)
library(survival)
library(smcure)

setwd("~/Weibull_cure_tidy_up_code/")
source("gmifs.R")
source("EM.R")
source("cmix.R")
source("Weibull.R")
source("glmnet_cox.R")
# source("Cox_EM.R")
source("evaluation.R")

source("~/result220201_Gompertz/Cox_EM_0304.R")

source("~/result220201_Gompertz/data_generation.R")

if(method=="cmix"){
  setwd("~/c-mix/C-mix/")
}

set.seed(seed)

ncores = 10
nIteration = 10

registerDoMC(ncores)
start.time = Sys.time()
if(method %in% c("cmix","glmnet_Cox")){ # do not use parallel computing
  res = list()
  for (i in 1:nIteration){
    data = data.gener(N=400, J=J, nonp=2, nTrue=10, A=A, rho=0.2, itct_mean, cens_ub,
                      alpha=1, lambda=2, same_signs, 
                      c("Weibull","GG","Gompertz","nonparametric","GG_baseline")[model_ind])
    Tr_data = data$Tr_data
    Te_data = data$Te_data
    if(knockoff){ # with knockoff
      if(method == "cmix"){
        X_k = create.second_order(Tr_data$X_p, method = "asdp", shrink = T)
        Xaug = cbind(Tr_data$X_p, X_k)
        X = cbind(Tr_data$X_u, Xaug)
        Y = pmax(Tr_data$time*365.25,1)
        delta = Tr_data$delta
        py_run_string(py_string)
        Z_b = abs(py$coeffs[-(1:(ncol(Tr_data$X_u)+1))])
        orig = 1:J
        W_b = Z_b[orig] - Z_b[orig+J]
        T_b = knockoff.threshold(W_b, fdr = 0.2, offset = 1)  # data-dependent threshold
        selected_b = (1:J)[W_b > T_b]
        res[[i]] = matrix(c(fdp(selected_b, Tr_data$nonzero_b),
                            power(selected_b, Tr_data$nonzero_b), 0, 0),
                          ncol = 1)
      } else if(method == "glmnet_Cox"){
        W_k = create.second_order(Tr_data$W_p, method = "asdp", shrink = T)
        Waug = cbind(Tr_data$W_p, W_k)
        fit = cox_glmnet(Tr_data$W_u, Waug, Tr_data$time, Tr_data$delta)
        Z_beta = abs(fit$beta_p)
        orig = 1:J
        W_beta = Z_beta[orig] - Z_beta[orig+J]
        T_beta = knockoff.threshold(W_beta, fdr = 0.2, offset = 1)  # data-dependent threshold
        selected_beta = (1:J)[W_beta > T_beta]
        res[[i]] = matrix(c(0, 0, fdp(selected_beta, Tr_data$nonzero_beta),
                            power(selected_beta, Tr_data$nonzero_beta)),
                          ncol = 1)
      }
      
    }else{ # without knockoff
      if (method == "cmix"){
        X = cbind(Tr_data$X_u, Tr_data$X_p)
        Y = pmax(Tr_data$time*365.25,1)
        delta = Tr_data$delta
        py_run_string(py_string)
        output = py$coeffs
        res[[i]] = list(data = data, output = output)
      } else if(method == "glmnet_Cox"){
        output = cox_glmnet(Tr_data$W_u, Tr_data$W_p, Tr_data$time, Tr_data$delta)
        res[[i]] = list(data = data, output = output)
      }
    }
  }
}else{ # use parallel computing
  res = foreach(i = 1:nIteration, .packages = c("knockoff","Rcpp")) %dopar% {
    data = data.gener(N=400, J=J, nonp=2, nTrue=10, A=A, rho=0.2, itct_mean, cens_ub,
                      alpha=1, lambda=2, same_signs, 
                      c("Weibull","GG","Gompertz","nonparametric","GG_baseline")[model_ind])
    Tr_data = data$Tr_data
    Te_data = data$Te_data
    
    if(knockoff){ # with knockoff
      
      X_k = create.second_order(Tr_data$X_p, method = "asdp", shrink = T)
      Xaug = cbind(Tr_data$X_p, X_k)
      Waug = cbind(Tr_data$W_p, X_k)
      if(method == "SCinCRM") X = cbind(Tr_data$X_u, Xaug)
      
      if(method=="gmifs"){
        fit = gmifs_fit(Tr_data$X_u, Xaug, Tr_data$W_u, Waug, Tr_data$time, Tr_data$delta, 
                        cure_cutoff=cure_cutoff, #### temporary
                        nIter=1e4, n_folds=5)
        Z_b = abs(fit$b_p_path[nrow(fit$b_p_path),])
        Z_beta = abs(fit$beta_p_path[nrow(fit$beta_p_path),])
      } else if(method=="EM"){
        fit = em_fit(Tr_data$X_u, Xaug, Tr_data$W_u, Waug, Tr_data$time, Tr_data$delta, nIter=100, n_folds=5, n_mu=50)
        Z_b = abs(fit$b_p)
        Z_beta = abs(fit$beta_p)
      } else if(method=="SCinCRM"){
        setwd("~/result210210_SCinCRM/SCinCRM/")
        source("highCoxCure.r")
        fit = highCoxCure(Tr_data$time, Tr_data$delta, X, fold=3,
                          nlambda1 = 15, nlambda2 = 15, lambda.min = 0.01,
                          lambda3 = c(0, 0.001, 0.01, 0.05, 0.1),	
                          nu=4, xi=1/sqrt(length(Tr_data$time)), eps=10^-4, maxitCD=1, maxitEM=10^3)
        Z_b = abs(fit$gam.hat.gcv[-(1:(ncol(Tr_data$X_u)+1))])
        Z_beta = abs(fit$bet.hat.gcv[-(1:ncol(Tr_data$X_u))])
      } else if(method=="Weibull"){
        fit = Weibull_fit(Tr_data$W_u, Waug, Tr_data$time, Tr_data$delta, nIter=1e4, n_folds=5)
        Z_b = rep(0, 2*J)
        Z_beta = abs(fit$beta_p)
      } else if(method=="Cox_EM"){
        fit = cox_em_fit(Tr_data$X_u, Xaug, Tr_data$W_u, Waug, Tr_data$time, Tr_data$delta, 
                         nIter=50, n_folds=5, n_mu=10)
        Z_b = abs(fit$b_p)
        Z_beta = abs(fit$beta_p)
      }
      orig = 1:J
      W_b = Z_b[orig] - Z_b[orig+J]
      W_beta = Z_beta[orig] - Z_beta[orig+J]
      T_b = knockoff.threshold(W_b, fdr = 0.2, offset = 1)  # data-dependent threshold
      T_beta = knockoff.threshold(W_beta, fdr = 0.2, offset = 1)  # data-dependent threshold
      selected_b = (1:J)[W_b > T_b]
      selected_beta = (1:J)[W_beta > T_beta]
      output = matrix(c(fdp(selected_b, Tr_data$nonzero_b),
                        power(selected_b, Tr_data$nonzero_b),
                        fdp(selected_beta, Tr_data$nonzero_beta),
                        power(selected_beta, Tr_data$nonzero_beta)), 
                      ncol = 1)
      rownames(output) = c("FDP_b","Power_b", "FDP_beta","Power_beta")
      return(output)
      
    } else{ # without knockoff
      
      if(method=="gmifs"){
        output = gmifs_fit(Tr_data$X_u, Tr_data$X_p, Tr_data$W_u, Tr_data$W_p, 
                           Tr_data$time, Tr_data$delta, 
                           cure_cutoff=cure_cutoff, #### temporary
                           nIter=1e4, n_folds=5)
      } else if(method=="EM"){
        output = em_fit(Tr_data$X_u, Tr_data$X_p, Tr_data$W_u, Tr_data$W_p, 
                        Tr_data$time, Tr_data$delta, nIter=100, n_folds=5, n_mu=50)
      } else if(method=="SCinCRM"){
        setwd("~/result210210_SCinCRM/SCinCRM/")
        source("highCoxCure.r")
        output = highCoxCure(Tr_data$time, Tr_data$delta, cbind(Tr_data$X_u, Tr_data$X_p), fold=3,
                             nlambda1 = 15, nlambda2 = 15, lambda.min = 0.01,
                             lambda3 = c(0, 0.001, 0.01, 0.05, 0.1),	
                             nu=4, xi=1/sqrt(length(Tr_data$time)), eps=10^-4, maxitCD=1, maxitEM=10^3)
      } else if(method=="Weibull"){
        output = Weibull_fit(Tr_data$W_u, Tr_data$W_p, Tr_data$time, Tr_data$delta, nIter=1e4, n_folds=5)
      } else if(method=="Cox_EM"){
        output = cox_em_fit(Tr_data$X_u, Tr_data$X_p, Tr_data$W_u, Tr_data$W_p, 
                        Tr_data$time, Tr_data$delta, nIter=50, n_folds=5, n_mu=10)
      }
      return(list(data = data, output = output))
      
    }
    
  } 
}
end.time = Sys.time()
end.time - start.time

ko_flag = ifelse(knockoff, "_ko", "")
setwd(paste0("~/result210901_data_generation/",method,ko_flag))
sign_flag = ifelse(same_signs, "_sign","")
model_flag = switch(model_ind, "", "_GG", "_Gompertz", "_np", "_ggbase")
itct_flag = paste0("_itct", itct_mean*10)
cens_flag = paste0("_cens", cens_ub/10)
save(res, file = paste0(method, ko_flag, sign_flag, model_flag, itct_flag, cens_flag,
                        "_A",A,"_",seed,".RData"))
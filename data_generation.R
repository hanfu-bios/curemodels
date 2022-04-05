
################## Data Generating Procedure ################## 
library(flexsurv)

data.gener <- function(N=400, J=500, nonp=2, nTrue=10, A=1, rho=0.5, itct_mean=0.5, cens_ub = 20,
                       alpha=1, lambda=2, same_signs=FALSE, model="Weibull"){
  # J: number of penalized covariates
  # nonp: number of non-penalized covariates
  tr_i = 1:round((3/4)*N) # training index
  te_i = (round((3/4)*N)+1):N # testing index  ##corrected Mar30
  
  ## covariance matrices
  sd = 0.5
  block_sz = round(J/nTrue)
  corr_X_p = matrix(0, J, J)
  for (i in 1:nTrue) 
    corr_X_p[(block_sz*(i-1)+1):(block_sz*i),(block_sz*(i-1)+1):(block_sz*i)] = 
    rho^abs(outer(1:block_sz, 1:block_sz, "-"))
  Sigma_X_p = sd^2*corr_X_p
  X_p = mvnfast::rmvn(N, mu=rep(0,J), sigma = Sigma_X_p)
  X_u = matrix(rnorm(N*nonp, sd=sd),ncol = nonp)
  W_p = X_p
  W_u = X_u
  
  if(!same_signs){ # true signals from two parts are randomly generated
    ## unpenalized 
    b_u = rnorm(nonp, mean = 0.3, sd = 0.1) 
    b_u = b_u * sample(c(1,-1), length(b_u), replace = T)
    beta_u = rnorm(nonp, mean = 0.3, sd = 0.1)
    beta_u = beta_u * sample(c(1,-1), length(beta_u), replace = T)
    
    ## penalized
    nonzero_b <- nonzero_beta <- rep(NA, nTrue)
    for(i in 1:nTrue){
      nonzero_b[i] = sample((block_sz*(i-1)+1):(block_sz*i), 1)
      nonzero_beta[i] = sample((block_sz*(i-1)+1):(block_sz*i), 1)
    } 
    b_p <- beta_p <- rep(0, J)
    b_p[nonzero_b] = A * sample(c(1,-1), nTrue, replace = T)
    beta_p[nonzero_beta] = A * sample(c(1,-1), nTrue, replace = T)
  } else { # true signals from two parts have same signs
    ## unpenalized 
    sign_u = sample(c(1,-1), nonp, replace = T)
    b_u = abs(rnorm(nonp, mean = 0.3, sd = 0.1)) * sign_u
    beta_u = abs(rnorm(nonp, mean = 0.3, sd = 0.1)) * sign_u
    ## penalized
    nonzero_b <- rep(NA, nTrue)
    for(i in 1:nTrue) nonzero_b[i] = sample((block_sz*(i-1)+1):(block_sz*i), 1) 
    b_p  <- rep(0, J)
    b_p[nonzero_b] = A * sample(c(1,-1), nTrue, replace = T)
    beta_p = b_p
    nonzero_beta = nonzero_b
  }
  
  ## cure or not
  itct = rnorm(1, mean = itct_mean, sd = 0.1)
  bx = itct + X_u %*% b_u + X_p %*% b_p
  p = 1/(1+exp(-bx))
  Y = rbinom(N, 1, p)
  
  # eps = rlogis(N, scale = 1)
  # itct = rnorm(1, mean = 0, sd = 0.1)
  # p = 1/(1+exp(- (itct + X_u %*% b_u + X_p %*% b_p) - eps))  # uncure probability
  # Y = ifelse(p > quantile(p,1-uncure), 1, 0) 
  
  ## survival
  beta_w = W_u%*%beta_u + W_p%*%beta_p
  if(model=="Weibull"){
    t = rweibull(N, shape = alpha, scale = 1/lambda * exp(- beta_w/alpha))
  }else if(model=="GG"){
    t = rgengamma(N, mu=-log(lambda)-beta_w/alpha, sigma=1/alpha, 
                  Q=2)
  }else if(model=="Gompertz"){
    t = rgompertz(N, shape = 0.2, rate = exp(beta_w))
  }else if(model=="nonparametric"){
    t = nonparametric_time_generator(N, beta_w, maxT = cens_ub, knots = 8)
  }else if(model=="GG_baseline"){
    t = rep(NA, N)
    for(i in 1:N){
      u = runif(1)
      t[i] = uniroot(function(a) 
        pgengamma(a, mu=-log(2), sigma=1, Q=0.5, lower.tail = FALSE)^exp(beta_w[i])-u, 
        c(0,20), extendInt = "yes")$root
    }
  }
  u = runif(N, 0, cens_ub)
  delta = ifelse(t>u | Y==0, 0, 1)
  time = pmin(t, u)
  time[Y==0] = u[Y==0]
  
  ## training and test
  Tr_data = list(X_u=X_u[tr_i,],X_p=X_p[tr_i,],W_u=W_u[tr_i,],W_p=W_p[tr_i,],time=time[tr_i],Y=Y[tr_i],
                 delta=delta[tr_i],nonzero_b=nonzero_b, nonzero_beta=nonzero_beta, 
                 b_u = b_u, beta_u = beta_u, b_p_nz = b_p[nonzero_b], beta_p_nz = beta_p[nonzero_beta], 
                 itct=itct)
  Te_data = list(X_u=X_u[te_i,],X_p=X_p[te_i,],W_u=W_u[te_i,],W_p=W_p[te_i,],time=time[te_i],Y=Y[te_i],delta=delta[te_i])
  return(list(Tr_data=Tr_data, Te_data=Te_data))
}

nonparametric_time_generator <- function (N, beta_w, maxT = 20, knots = 8) 
{
  time <- 0:maxT
  k <- c(0, sort(sample(time[2:maxT], size = knots, replace = FALSE)), maxT)
  heights <- c(0, sort(1-pmin(rexp(knots, rate = 10),1-1e-15)), 1)
  # heights <- c(0, sort(runif(knots)), 1)
  tk <- merge(data.frame(time), data.frame(time = k, heights), 
              by = "time", all = FALSE)
  MonotonicSpline <- stats::splinefun(x = tk$time, y = tk$heights, 
                                      method = "hyman")
  t = rep(NA, N)
  for(i in 1:N){
    u = runif(1)
    t[i] = uniroot(function(a) (1-MonotonicSpline(a))^exp(beta_w[i])-u, c(0,maxT))$root
  }
  return(t)
}

# dev.new(width=10,height=4,noRStudioGD = TRUE)
# par(mfrow=c(1,2), mar=c(5,5,4,2)+0.1)
# time = seq(0,20,0.1)
# plot(time, 1-MonotonicSpline(time), type = 'l', xlab = "Time", ylab = "Survival", lwd = 2) # survival
# plot(time[-length(time)], diff(MonotonicSpline(time))/(1-MonotonicSpline(time))[-length(time)], type = 'l',
#      xlab = "Time", ylab = "Hazard", lwd = 2) # hazard
# setwd("~/Documents/Research/CoxMCM/plot/")
# dev2bitmap("nonparametric_Surv_Haz_exp20_J100.jpeg", res = 750, height = 4, width=10)
# plot(diff(MonotonicSpline(time)))

# dev.new(width=10,height=4,noRStudioGD = TRUE)
# par(mfrow=c(1,2), mar=c(5,5,4,2)+0.1)
# time = seq(0.01,10,0.01)
# plot(time, pgengamma(time, mu=-log(2), sigma = 1, Q=0.5, lower.tail = FALSE), type = 'l', xlab = "Time", ylab = "Survival", lwd = 2) # survival
# plot(time, dgengamma(time, mu=-log(2), sigma = 1, Q=0.5)/
#                 pgengamma(time, mu=-log(2), sigma = 1, Q=0.5, lower.tail = FALSE), type = 'l',
#      xlab = "Time", ylab = "Hazard", lwd = 2) # hazard
# setwd("~/Documents/Research/CoxMCM/plot/")
# dev2bitmap("ggbase_Surv_Haz.jpeg", res = 750, height = 4, width=10)

# plot_gg <- function(mu, sigma, Q){
#   time = seq(0.01,10,0.01)
#   plot(time, dgengamma(time, mu = mu, sigma = sigma, Q=Q)/
#          pgengamma(time, mu=mu, sigma = sigma, Q=Q, lower.tail = FALSE), 
#        type="l",ylab="Hazard", main = paste0("sigma=",sigma,", Q=", Q))
# }
# 
# par(mfrow=c(2,2), mar=c(2,4,4,2)+0.1)
# for(Q in 1/(2:5)) plot_gg(mu = -log(2), sigma = 1, Q)

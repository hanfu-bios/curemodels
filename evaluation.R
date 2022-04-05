
################## Evaluation Metrics ################## 
library(survival)

C.stat <- function(cure_cutoff = 5, b_u_hat, itct_hat, b_p_hat, beta_u_hat, beta_p_hat, testing_delta, 
                   testing_time, X_u, X_p, W_u, W_p){
  C_csw_num = 0
  C_csw_denom = 0
  testing_n = length(testing_time)
  v = rep(0, testing_n)
  y = rep(999, testing_n)
  y[testing_time>cure_cutoff] = 0
  y[testing_time<=cure_cutoff & testing_delta==1] = 1
  v[y<2] = 1
  if(all(b_p_hat==0)) {
    p_hat = 1/(1+exp(-itct_hat-X_u %*% b_u_hat))
  } else{
    p_hat = 1/(1+exp(-itct_hat- X_u %*% b_u_hat - X_p[,b_p_hat!=0,drop=FALSE] %*% b_p_hat[b_p_hat!=0])) # uncure probability
  }
  temp = v*y + (1-v)*as.vector(p_hat)
  if(all(beta_p_hat==0)) {
    W_beta = W_u %*% beta_u_hat 
  } else{
    W_beta = W_u %*% beta_u_hat + W_p[,beta_p_hat!=0,drop=FALSE] %*% beta_p_hat[beta_p_hat!=0] 
  }
  for(i in 1:testing_n)
    for(j in 1:testing_n){
      if (j==i | !testing_delta[i] | testing_time[i]>testing_time[j]) next
      I_ij = testing_time[i]<testing_time[j] | (testing_time[i]==testing_time[j] & !testing_delta[j])
      if (!I_ij) next
      if (W_beta[i]>W_beta[j]) C_csw_num = C_csw_num + temp[j]
      C_csw_denom = C_csw_denom + temp[j]
    }
  return(C_csw_num / C_csw_denom)
}

c_stat_beta <- function(beta_u_hat,  beta_p_hat, testing_delta, testing_time, W_u, W_p){
  C_csw_num = 0
  C_csw_denom = 0
  testing_n = length(testing_time)
  if(all(beta_p_hat==0)) {
    W_beta = W_u %*% beta_u_hat 
  } else{
    W_beta = W_u %*% beta_u_hat + W_p[,beta_p_hat!=0,drop=FALSE] %*% beta_p_hat[beta_p_hat!=0] 
  }
  for(i in 1:testing_n)
    for(j in 1:testing_n){
      if (j==i | !testing_delta[i] | testing_time[i]>testing_time[j]) next
      I_ij = testing_time[i]<testing_time[j] | (testing_time[i]==testing_time[j] & !testing_delta[j])
      if (!I_ij) next
      if (W_beta[i]>W_beta[j]) C_csw_num = C_csw_num + 1
      C_csw_denom = C_csw_denom + 1
    }
  return(C_csw_num / C_csw_denom)
}

c_stat_b <- function(b_u_hat, itct_hat, b_p_hat, testing_delta, testing_time, X_u, X_p){
  C_csw_num = 0
  C_csw_denom = 0
  testing_n = length(testing_time)
  if(all(b_p_hat==0)) {
    xb = itct_hat+X_u %*% b_u_hat
  } else{
    xb = itct_hat+X_u %*% b_u_hat + X_p[,b_p_hat!=0,drop=FALSE] %*% b_p_hat[b_p_hat!=0]
  }
  for(i in 1:testing_n)
    for(j in 1:testing_n){
      if (j==i | !testing_delta[i] | testing_time[i]>testing_time[j]) next
      I_ij = testing_time[i]<testing_time[j] | (testing_time[i]==testing_time[j] & !testing_delta[j])
      if (!I_ij) next
      if (xb[i]>xb[j]) C_csw_num = C_csw_num + 1
      C_csw_denom = C_csw_denom + 1
    }
  return(C_csw_num / C_csw_denom)
}

# cure prediction accuracy
cure_accu <- function(b_u_hat, itct_hat, b_p_hat, X_u, X_p, delta, time, Y){
  true_p = mean(Y)
  if(all(b_p_hat==0)) {
    p_hat = 1/(1+exp(-itct_hat-X_u %*% b_u_hat))
  } else{
    p_hat = 1/(1+exp(-itct_hat- X_u %*% b_u_hat - X_p[,b_p_hat!=0,drop=FALSE] %*% b_p_hat[b_p_hat!=0])) # uncure probability
  }
  cure_out = summary(survfit(Surv(time, delta)~1), times = max(time[delta==1]))$surv
  Y_out = p_hat > quantile(p_hat,cure_out) + 0
  accu = mean(Y==Y_out)
  return(accu)
}

RME <- function(hat, true, star){
  num = matrix(hat-true,1)%*%Sigma_X_p%*%matrix(hat-true,ncol=1)
  denom = matrix(star-true,1)%*%Sigma_X_p%*%matrix(star-true,ncol=1)
  return(num/denom)
}
ERR <- function(hat, true, star){
  num = sum((hat-true)^2)
  denom = sum((star-true)^2)
  return(num/denom) 
}
RME_ERR <- function(b_p_hat, beta_p_hat, Tr_data, A){
  J = length(b_p_hat)
  Xoracle = Tr_data$X_p[,Tr_data$nonzero_b]
  Woracle = Tr_data$W_p[,Tr_data$nonzero_beta]
  out = weibull.cure(scale(Tr_data$X_u), scale(Xoracle),
                     scale(Tr_data$W_u), scale(Woracle),
                     Tr_data$time, Tr_data$delta,
                     nIter=1e4)
  star_b <- star_beta <- rep(0,J)
  nstep = length(out$logLikelihood)
  star_b[Tr_data$nonzero_b] = out$b_p_path[nstep,]
  star_beta[Tr_data$nonzero_beta] = out$beta_p_path[nstep,]
  
  true_b <- true_beta <- rep(0,J)
  true_b[Tr_data$nonzero_b] = Tr_data$b_p_nz * A
  true_beta[Tr_data$nonzero_beta] = Tr_data$beta_p_nz * A
  
  return(c(RME(b_p_hat, true_b, star_b),
           ERR(b_p_hat, true_b, star_b),
           RME(beta_p_hat, true_beta, star_beta),
           ERR(beta_p_hat, true_beta, star_beta)))
}

fdp <- function(selected, nonzero) {
  if (length(selected)==0) return(0)
  else return(sum(!selected %in% nonzero) / max(1, length(selected)))
}
power <- function(selected, nonzero) {
  if (length(selected)==0) return(0)
  else return(sum(selected %in% nonzero) / length(nonzero))
}
fdp_power <- function(nonzero_b, nonzero_beta, b_p_hat, beta_p_hat){
  selected_b = which(b_p_hat!=0)
  selected_beta = which(beta_p_hat!=0)
  return(c(fdp(selected_b, nonzero_b),
           power(selected_b, nonzero_b),
           fdp(selected_beta, nonzero_beta),
           power(selected_beta, nonzero_beta)))
}

AUC_msi <- function(cure_cutoff = 5, b_u_hat, itct_hat, b_p_hat, testing_delta, testing_time, X_u, X_p){
  testing_n = length(testing_time)
  v = rep(0, testing_n)
  y = rep(999, testing_n)
  y[testing_time>cure_cutoff] = 0
  y[testing_time<=cure_cutoff & testing_delta==1] = 1
  v[y<2] = 1
  if(all(b_p_hat==0)) {
    xb = itct_hat+X_u %*% b_u_hat
  } else{
    xb = itct_hat+X_u %*% b_u_hat + X_p[,b_p_hat!=0,drop=FALSE] %*% b_p_hat[b_p_hat!=0]
  }
  p_hat = 1/(1+exp(-xb))
  temp = v*y + (1-v)*p_hat
  temp1 = temp[order(p_hat, decreasing = T)]
  temp_f = v*(1-y) + (1-v) * (1-p_hat)
  temp1f = temp_f[order(p_hat, decreasing = T)]
  TPR = c(0, cumsum(temp1)/cumsum(temp1)[testing_n])
  FPR = c(0, cumsum(temp1f)/cumsum(temp1f)[testing_n])
  height = (TPR[-1]+TPR[-length(TPR)])/2
  width = diff(FPR)
  auc_msi = sum(height*width)
  return(auc_msi)
}

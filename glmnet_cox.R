
################## Cox - glmnet ################## 
library(glmnet)
library(survival)

cox_glmnet <- function(W_u, W_p, time, delta){
  fit_lat = cv.glmnet(x = cbind(W_u, W_p), 
                      y = Surv(time = time, event = delta),
                      family = "cox",
                      penalty.factor = c(rep(0, ncol(W_u)), rep(1, ncol(W_p))),
                      type.measure = "C")
  coef_lat = fit_lat$glmnet.fit$beta[,which.max(fit_lat$cvm)]
  beta_u = coef_lat[1:ncol(W_u)]
  beta_p = coef_lat[-(1:ncol(W_u))]
  output = list(beta_u = beta_u, beta_p = beta_p)
}
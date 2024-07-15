library(dplyr)
library(quantreg)
library(lqa)
library(BMisc)
library(glmnet)
library(Formula)
library(ggplot2)
library(MASS)
library(qte)
library(ggsci)
library(reshape2)

# Creation of simulated data
dat1 <- function(n, p, gamma, beta_tre, beta_untre, eta, alpha, senario_num) {
  x <- matrix(rnorm(p * n, 0, 1), byrow = FALSE, nrow = n)
  ga_z <- rowSums((x * matrix(alpha, nrow = n, ncol = p, byrow = TRUE)))
  pa <- expit(ga_z)
  z <- as.numeric(runif(n = length(pa)) < pa)
  
  if (senario_num == 1) {
    y <- eta * z + x[, 1:p] %*% beta_tre[1:p] + rnorm(nrow(x), 0, 1)
    data <- data.frame(x = x, y = y, z = z)
  } else if (senario_num == 2) {
    y <- eta * z + x[, 1:p] %*% beta_tre[1:p] + rnorm(nrow(x), 0, 1) * (1 + gamma * (x[, 2] + x[, 10]))
    data <- data.frame(x = x, y = y, z = z)
  } else if (senario_num == 3) {
    y <- eta * z * (1 + x[, 2]) + x[, 1:p] %*% beta_tre[1:p] + rnorm(nrow(x), 0, 1) * (1 + gamma * (x[, 2] + x[, 10]))
    data <- data.frame(x = x, y = y, z = z)
  } else {
    x_z <- cbind(x, z)
    x1 <- as.matrix(subset(x_z, z == 1))
    x0 <- as.matrix(subset(x_z, z == 0))
    y1 <- eta + x1[, 1:p] %*% beta_tre[1:p] + rnorm(nrow(x1), 0, 1) * (1 + gamma * (x1[, 2] + x1[, 10]))
    y0 <- x0[, 1:p] %*% beta_untre[1:p] + rnorm(nrow(x0), 0, 1) * (1 + gamma * (x0[, 2] + x0[, 10]))
    data <- data.frame(rbind(cbind(x1[, 1:p], y1, x1[, p + 1]), cbind(x0[, 1:p], y0, x0[, p + 1])))
  }
  
  return(data)
}

#Creation of true data
dat2 <- function(n, p, gamma, beta_tre, beta_untre, eta, senario_num) {
  x <- matrix(rnorm(p * n, 0, 1), byrow = FALSE, nrow = n)
  
  if (senario_num == 1) {
    y1 <- eta + x[, 1:p] %*% beta_tre[1:p] + rnorm(nrow(x), 0, 1)
    y0 <- x[, 1:p] %*% beta_tre[1:p] + rnorm(nrow(x), 0, 1)
    true_data <- cbind(y1, y0)
  } else if (senario_num == 2) {
    y1 <- eta + x[, 1:p] %*% beta_tre[1:p] + rnorm(nrow(x), 0, 1) * (1 + gamma * (x[, 2] + x[, 10]))
    y0 <- x[, 1:p] %*% beta_tre[1:p] + rnorm(nrow(x), 0, 1) * (1 + gamma * (x[, 2] + x[, 10]))
    true_data <- cbind(y1, y0)
  } else if (senario_num == 3) {
    y1 <- eta * (1 + x[, 2]) + x[, 1:p] %*% beta_tre[1:p] + rnorm(nrow(x), 0, 1) * (1 + gamma * (x[, 2] + x[, 10]))
    y0 <- x[, 1:p] %*% beta_tre[1:p] + rnorm(nrow(x), 0, 1) * (1 + gamma * (x[, 2] + x[, 10]))
    true_data <- cbind(y1, y0)
  } else {
    y1 <- eta + x[, 1:p] %*% beta_tre[1:p] + rnorm(nrow(x), 0, 1) * (1 + gamma * (x[, 2] + x[, 10]))
    y0 <- x[, 1:p] %*% beta_untre[1:p] + rnorm(nrow(x), 0, 1) * (1 + gamma * (x[, 2] + x[, 10]))
    true_data <- cbind(y1, y0)
  }
  
  return(true_data)
}

# Calculation of QTE
calcu_qte = function(y, z, ps, tau) {
  n <- length(z)
  treated.weights <- z / (n * ps)
  untreated.weights <- (1 - z) / (n * (1 - ps))
  
  treated.firpo.quantiles <- getWeightedQuantiles(tau, y, treated.weights, norm = TRUE)
  untreated.firpo.quantiles <- getWeightedQuantiles(tau, y, untreated.weights, norm = TRUE)
  est_qte <- treated.firpo.quantiles - untreated.firpo.quantiles
  return(est_qte)
}

# logit function
expit = function(x) {
  pr = (exp(x) / (1 + exp(x)))
  return(pr)
}

create_weights = function(fp, fA, fw) {
  fw = (fp)^(-1)
  fw[fA == 0] = (1 - fp[fA == 0])^(-1)
  return(fw)
}

# caluculate wAMD
wAMD_func = function(DataM, varlist, trt.var, wgt, beta) {
  trt = untrt = diff_vec = rep(NA, length(beta))
  names(trt) = names(untrt) = names(diff_vec) = varlist
  
  for (j in 1:length(varlist)) {
    this.var = paste("w", varlist[j], sep = "")
    DataM[, this.var] = DataM[, varlist[j]] * DataM[, wgt]
    trt[j] = sum(DataM[DataM[, trt.var] == 1, this.var]) / sum(DataM[DataM[, trt.var] == 1, wgt])
    untrt[j] = sum(DataM[DataM[, trt.var] == 0, this.var]) / sum(DataM[DataM[, trt.var] == 0, wgt])
    diff_vec[j] = abs(trt[j] - untrt[j])
  }
  
  wdiff_vec = diff_vec * abs(beta)
  wAMD = c(sum(wdiff_vec))
  ret = list(diff_vec = diff_vec, wdiff_vec = wdiff_vec, wAMD = wAMD)
  
  return(ret)
}

# Assume that the dataset has covariates, treatment variables and outcome.
# Columns 1 to p are covariates, p+1 columns are outcomes, p+2 columns are treatment variables
qoal_func = function(data, tau, lambda_vec, gamma_conv){
  # data formatting
  p = ncol(data) - 2
  n = nrow(data)
  colnames(data) = c(paste("X",1:p,sep=""),"Y","Z")
  
  # Outcome regression 
  form_y <- formula(paste0(paste("Y",paste(paste0("X", 1:p), collapse = "+"), sep = "~"),"+Z"))
  betaXY <- coef(rq(form_y, data = data, tau))[2:(p + 1)]

  
  # Calculation of regularization parameters
  gamma_vals = 2*( gamma_conv - lambda_vec + 1 )
  names(lambda_vec) = as.character(lambda_vec)
  names(gamma_vals) = names(lambda_vec)
  
  # Empty array creation
  wAMD_vec = rep(NA, length(lambda_vec))
  PS_pro = as.data.frame(matrix(NA,nrow=n,ncol=length(lambda_vec)))
  coeff_XA = as.data.frame(matrix(NA,nrow = p,ncol=length(lambda_vec)))
  names(wAMD_vec) = names(coeff_XA) = names(PS_pro)= names(lambda_vec)
  rownames(coeff_XA) = paste("X",1:p,sep="")
  
  w.full.form = formula(paste("Z~",paste(paste("X",1:p,sep=""),collapse="+")))
  for( lil in names(lambda_vec) ){
    il = lambda_vec[lil]
    ig = gamma_vals[lil]
    
    # Setting Adaptive Lasso weights
    oal_pen = adaptive.lasso(lambda = n^(il),al.weights = abs(betaXY)^(-ig) )
    # Running adaptive lasso
    logit_oal = lqa.formula(w.full.form, data = data, penalty=oal_pen, family=binomial(logit))
    
    data[,paste("f.pA",lil,sep="")] = predict.lqa(logit_oal)$mu.new
    PS_pro[,lil]= predict.lqa(logit_oal)$mu.new
    coeff_XA[,lil] = coef(logit_oal)[2:(p+1)]
    data[,paste("w",lil,sep="")] = create_weights(fp=data[,paste("f.pA",lil,sep="")],fA=data$Z)
    wAMD_vec[lil] = wAMD_func(DataM=data,varlist=paste("X",1:p,sep=""),trt.var="Z", wgt=paste("w",lil,sep=""),beta=betaXY)$wAMD
  } 
  
  tt = which.min(wAMD_vec)
  # print out the coefficients for the propensity score that corresponds with smalles wAMD value
  ps = PS_pro[,tt]
  coef = coeff_XA[,tt]
  qoal = list(coef = coef, best_lambda = lambda_vec[tt], ps = ps)
  
  return(qoal)
}



# Assume that the dataset has covariates, treatment variables and outcome
# Columns 1 to p are covariates, p+1 columns are outcomes, p+2 columns are treatment variables
weight_lasso = function(data, lambda_vec, gamma_conv, method){
  # data formatting
  p = ncol(data) - 2
  n = nrow(data)
  colnames(data) = c(paste("X",1:p,sep=""),"Y","Z")
  
  # Outcome regression 
  if (tolower(method) %in% c("oal", "adl")) {
    if (tolower(method) == "oal") {
      form_y <- formula(paste0(paste("Y",paste(paste0("X", 1:p), collapse = "+"), sep = "~"),"+Z"))
    } else {  # methods == "adl"
      form_y <- formula(paste0(paste("Z",paste(paste0("X", 1:p), collapse = "+"), sep = "~")))
    }
    
    betaXY <- coef(lm(form_y, data = data))[2:(p + 1)]
  } else {
    print("The method needs to include 'adl' or 'oal'.")
  }

  # Calculation of regularization parameters
  gamma_vals = 2*( gamma_conv - lambda_vec + 1 )
  names(lambda_vec) = as.character(lambda_vec)
  names(gamma_vals) = names(lambda_vec)
  
  # Empty array creation
  wAMD_vec = rep(NA, length(lambda_vec))
  PS_pro = as.data.frame(matrix(NA,nrow=n,ncol=length(lambda_vec)))
  coeff_XA = as.data.frame(matrix(NA,nrow = p,ncol=length(lambda_vec)))
  names(wAMD_vec) = names(coeff_XA) = names(PS_pro)= names(lambda_vec)
  rownames(coeff_XA) = paste("X",1:p,sep="")
  
  w.full.form = formula(paste("Z~",paste(paste("X",1:p,sep=""),collapse="+")))
  for( lil in names(lambda_vec) ){
    il = lambda_vec[lil]
    ig = gamma_vals[lil]
    
    # Setting Adaptive Lasso weights
    oal_pen = adaptive.lasso(lambda = n^(il),al.weights = abs(betaXY)^(-ig) )
    # Running adaptive lasso
    logit_oal = lqa.formula(w.full.form, data = data, penalty=oal_pen, family=binomial(logit))
    
    data[,paste("f.pA",lil,sep="")] = predict.lqa(logit_oal)$mu.new
    PS_pro[,lil]= predict.lqa(logit_oal)$mu.new
    coeff_XA[,lil] = coef(logit_oal)[2:(p+1)]
    data[,paste("w",lil,sep="")] = create_weights(fp=data[,paste("f.pA",lil,sep="")],fA=data$Z)
    wAMD_vec[lil] = wAMD_func(DataM=data,varlist=paste("X",1:p,sep=""),trt.var="Z", wgt=paste("w",lil,sep=""),beta=betaXY)$wAMD
  } 
  
  tt = which.min(wAMD_vec)
  # print out the coefficients for the propensity score that corresponds with smalles wAMD value
  ps = PS_pro[,tt]
  coef = coeff_XA[,tt]
  res_wlasso = list(coef = coef, best_lambda = lambda_vec[tt], ps = ps)
  return(res_wlasso)
}

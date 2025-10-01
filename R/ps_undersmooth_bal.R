
#' Undersmoothing Lasso PS Models using balance diagnostics
#'
#' @param data a dataset or matrix containing baseline covariates
#' @param treatment a binary vector for treatment
#' @param ps_dat a matrix of fitted propensity scores
#' @param normalized a boolean TRUE/FALSE to use normalized weighting
#' @param standardize a boolean TRUE/FALSE to use standardized differences
#'
#' @returns A list containing: 1) predicted values from the selected lasso model; 2) a vector with the value for the smallest average standardized difference among all models after PS weighting and the value of the maximum standardized difference among all models after PS weighting; 3) the index (column) of the selected model/predicted values for each balance criteria.
#' @details The ps_undersmooth_bal() function takes as input predicted values from several propensity score models. The function calculates the weighted standardized difference for each variable and returns the index (column) for the models where 1) the average standardized difference is the smallest, and 2) the maximum weighted standardized difference is the smallest. 
#' @export
#' @examples
#' #loading package
#' library(PSLassoSynthNC)
#' 
#' #creating some simulated data for testing
#' nstudy<- 2000
#' nvars<- 500
#' nc<- 100
#' ns<- nvars-(nc)
#' alpha_temp<- runif(nc, 0.0, 0.4)
#' beta_temp<- runif(nc, 0.0, 0.4)
#' random_neg<- sample(1:length(alpha_temp), 0.5*length(alpha_temp), replace=FALSE)
#' alpha_temp[random_neg]<- -1*alpha_temp[random_neg]
#' beta_temp[random_neg]<-  -1*beta_temp[random_neg]
#' alpha<-  matrix(c(alpha_temp, rep(0, ns)), ncol=1)
#' beta<-   matrix(c(beta_temp, rep(0, ns)), ncol=1)
#' betaE<- 0
#' cprev<- runif(nvars, 0, 0.3)
#' cprev<- sample(cprev)
#' oprev<- 0.05
#' tprev<- 0.4
#' Xcovs_sim<- matrix(rnorm((nstudy*nvars), 0, 1), nrow=nstudy, ncol=nvars)
#' Xcovs_sim<- as.data.frame(Xcovs_sim)  
#' names(Xcovs_sim)<- c(paste0('x', 1:nvars))
#' W<- as.matrix(Xcovs_sim)
#' colnames(W)<- c(paste0('x', 1:nvars))
#' linear_pred_e<- W %*% alpha
#' linear_pred_y<- W %*% beta
#' treatment_inc<- tprev
#' fn <- function(c) mean(plogis(c + linear_pred_e)) - treatment_inc
#' alpha0 <- uniroot(fn, lower = -20, upper = 20)$root
#' Ee <- (1 + exp( -(alpha0 + linear_pred_e) ))^-1
#' e<- rbinom(nstudy, 1, Ee)
#' outcome_inc<- oprev
#' fn <- function(c) mean(plogis(c + betaE*e + linear_pred_y  )) - outcome_inc
#' beta0 <- uniroot(fn, lower = -20, upper = 20)$root
#' Ey <- (1 + exp( -( beta0 + betaE*e + linear_pred_y )))^-1
#' y<- rbinom(nstudy, 1, Ey)
#' simdat <- as.data.frame(cbind(y, e, Ee, Xcovs_sim))
#'
#' #creating folid vector for testing
#' N <- length(e)
#' V=10
#' n<- 1:length(e)
#' cvfolds<- stratifyCVFoldsByYandID(V=V, Y = e)
#' folds <- cvfolds$validRows
#' foldid <- cvfolds$fold_id
#'
#' #running treatment_model function
#' trt_out<- treatment_model(data=Xcovs_sim, treatment=e, foldid=foldid, alpha=1, lambda_ratio=NULL, nlambda=100, nmodels=NULL, maxit=5000, penalty=NULL, par=FALSE)
#' 
#' #extracting the out-of-fold predictions from trt_out
#' #' ps_dat_out<- trt_out[[2]] 
#'
#' #running ps_undersmooth_bal() function
#' ps_undersmooth_bal(data=Xcovs_sim, treatment=e, ps_dat=ps_dat_out, method='ow')
ps_undersmooth_bal<- function(data,
                              treatment,
                              ps_dat,
                              method,
                              normalized=TRUE,
                              standardize=TRUE){
  
  ## note: balance_select calculates balance with normalized weighted averages
  cov_diff<- balance_weighted_diff(data=data,
                                   treatment=treatment,
                                   ps_dat=ps_dat,
                                   method=method,
                                   normalized=normalized,
                                   standardize=standardize)
  
  ## exclude first column which is unadjusted (crude) differences
  cov_diff<- cov_diff[,-1]
  
  ## standardized absolute differences
  cov_diff_abs<- apply(cov_diff, 2, abs)
  
  ## balance metric 1: minimum ASAMD
  cov_diff_abs_avg<- apply(cov_diff_abs, 2, mean)
  bal_m1_index<- which.min(cov_diff_abs_avg)
  bal_m1_value<- cov_diff_abs_avg[bal_m1_index]
  
  ## balance metric 2: model with the smallest max standardized difference
  cov_diff_max<- apply(cov_diff_abs, 2, max)
  bal_m2_index<- which.min(cov_diff_max)
  bal_m2_value<- cov_diff_max[bal_m2_index]
  select_value<- c(bal_m1_value, bal_m2_value)
  select_index<- c(bal_m1_index, bal_m2_index)
  select_preds<- ps_dat[,select_index]
  names(select_preds)<- c("model1", "model2")
  names(select_value)<- c("asamd", "max_diff")
  names(select_index)<- c("index1", "index2")
  results<- list(select_preds,
                 select_value,
                 select_index)
  names(results)<- c("predictions", "balance", "index")
  
  return(results)
}
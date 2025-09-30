##################################################################################################
##
## ps_undersmooth_bal: function to choose lambda value based on minimizing balance criteria
##
##################################################################################################

#’ Undersmoothing Lasso PS Models using balance diagnostics
#’
#’ @param data a dataset or matrix containing baseline covariates
#’ @param treatment a binary vector for treatment
#’ @param ps_dat a matrix of fitted propensity scores
#’ @param normalized a boolean TRUE/FALSE to use normalized weighting
#’ @param standardize a boolean TRUE/FALSE to use standardized differences
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
##########################################################################################
##
## Functions to calculate covariate balance:
## weighted.var() and mw_fun() are helper functions used within weighted_diff()
## weighted_diff() is a helper function used within balance_weighted_diff()
## balance_weighted_diff() calculates covariate balance for each covariate
##
##########################################################################################
#Helper function used within weighted.diff
#weighted.var function from Gavin Simpson. URL: https://stat.ethz.ch/pipermail/r-help/2008-July/168762.html
weighted.var <- function(x, w, na.rm = FALSE) {
  if (na.rm) {
    w <- w[i <- !is.na(x)]
    x <- x[i]
  }
  sum.w <- sum(w)
  sum.w2 <- sum(w^2)
  mean.w <- sum(x * w) / sum(w)
  (sum.w / (sum.w^2 - sum.w2)) * sum(w * (x - mean.w)^2, na.rm = na.rm)
}
#Helper function used within weighted.diff
mw_fun<- function(score, treatment){
  score_data<- as.data.frame(cbind(score, (1-score)))
  numerator<- apply(score_data, 1, min)
  weight<- numerator / (treatment*score + (1-treatment)*(1-score))
  return(weight)
}


#’ Helper function used within ’balance_weighted_diff()’ to calculate weighted standardized differences
#’
#’ @param data a dataset or matrix containing baseline covariates
#’ @param data0 a dataset or matrix containing baseline covariates for unexposed group
#’ @param data1 a dataset or matrix containing baseline covariates for exposed group
#’ @param score a dataset or matrix of fitted propensity score values (each column corresponds to different model)
#’ @param treatment a vector of binary indicators indicating treatment status
#’ @param method weighting method used to calculate weighted standardized differences
#’ @param normalized boolean TRUE/FALSE to indicate use of normalized weights (default is TRUE)
weighted_diff<- function(data,
                         data0,
                         data1,
                         score,
                         treatment,
                         method,
                         normalized,
                         standardize){
  ## IPTW weights
  if(method=="iptw"){
    e1_mean<- mean(treatment)
    e0_mean<- 1-e1_mean
    #weight<- (treatment*e1_mean)/score + ((1-treatment)*e0_mean)/(1-score) #stabilized weight
    weight<- treatment/score + (1-treatment)/(1-score) #unstabilized weight
  }
  ## matching weights
  if(method=="mw"){
    weight<- mw_fun(score=score, treatment=treatment)
  }
  ## overlap weights
  if(method=="ow"){
    weight<- treatment*(1-score) + (1-treatment)*score
  }
  ## creating weighted cohorts by treatment group
  weight0<- weight[treatment==0]
  weight1<- weight[treatment==1]
  
  #############################################################
  ## normalized weighted average (weighted.mean normalizes)
  if(normalized==TRUE){
    fun0<- function(x){weighted.mean(x,weight0)}
    fun1<- function(x){weighted.mean(x,weight1)}
    fun0.sd<- function(x){sqrt(weighted.var(x, weight0))}
    fun1.sd<- function(x){sqrt(weighted.var(x, weight1))}
    #fun0.sd<- function(x){sqrt(wtd.var(x, weight0, na.rm=TRUE))}
    #fun1.sd<- function(x){sqrt(wtd.var(x, weight1, na.rm=TRUE))}
    mean0.w<- apply(data0, 2, fun0)
    mean1.w<- apply(data1, 2, fun1)
    ## same as above but calculated manually
    #mean0.w<- apply(data0, 2, function(x) sum(x*weight0)/sum(weight0))
    #mean1.w<- apply(data1, 2, function(x) sum(x*weight1)/sum(weight1))
    sd0.w<- apply(data0, 2, fun0.sd)
    sd1.w<- apply(data1, 2, fun1.sd)
  }
  
  ############################################################################
  ## unnormalized weighted balance (Austin does not recommend this approach)
  if(normalized != TRUE){
    mean0.w<- apply(data0, 2, function(x) mean(x*weight0))
    mean1.w<- apply(data1, 2, function(x) mean(x*weight1))
    sd0.w<- apply(data0, 2, fun0.sd)
    sd1.w<- apply(data1, 2, fun1.sd)
    #sd0.w<- apply(data1, 2, function(x) sd(x*weight0))
    #sd1.w<- apply(data0, 2, function(x) sd(x*weight1))
  }
  
  #######################################
  ## Calculating mean difference
  if(standardize==TRUE){diff_weight<- (mean0.w-mean1.w) / sqrt((sd0.w^2 + sd1.w^2)/2)}
  if(standardize==FALSE){diff_weight<- (mean0.w-mean1.w)}
  abs_diff_weight<- abs(diff_weight)
  return(diff_weight)
}


#’ Calculates unadjusted and weighted standardized differences for each covariate
#’
#’ @param data a dataset or matrix containing baseline covariates
#’ @param treatment a vector of binary indicators indicating treatment status
#’ @param ps_dat a datset or matrix of fitted propensity score values (each column corresponds to a different model)
#’ @param method weighting method used to calculate weighted standardized differences
#’ @param normalized boolean TRUE/FALSE to indicate use of normalized weights (default is TRUE)
balance_weighted_diff<- function(data,
                                 treatment,
                                 ps_dat,
                                 method,
                                 normalized,
                                 standardize){
  treatment=treatment
  data=data
  ps_dat=as.data.frame(ps_dat)
  method=method
  
  ###############################
  ## Unadjusted balance
  data0<- data[treatment==0,]
  data1<- data[treatment==1,]
  data_test0<- data0
  data_test1<- data1
  mean0<- apply(data_test0, 2, mean)
  mean1<- apply(data_test1, 2, mean)
  sd0<- apply(data_test0, 2, sd)
  sd1<- apply(data_test1, 2, sd)
  sdf<- apply(data, 2, sd)
  if(standardize==TRUE){diff_crude<- (mean0-mean1) / sqrt((sd0^2 + sd1^2)/2)}
  if(standardize==FALSE){diff_crude<- (mean0-mean1)}
  abs_diff_crude<- abs(diff_crude)
  
  #####################################
  ## PS weighted balance
  bal_avg<- NULL
  bal_max<- NULL
  dat_diff_weight<- NULL
  for(iii in 1:ncol(ps_dat)){
    ps_select<- ps_dat[,iii]
    ## weighted balance (weighted_diff is defined above)
    weighted_differences<- weighted_diff(data=data,
                                         data0=data0,
                                         data1=data1,
                                         score=ps_select,
                                         treatment=treatment,
                                         method=method,
                                         normalized=normalized,
                                         standardize=standardize)
    diff_weight<- weighted_differences
    dat_diff_weight<- cbind(dat_diff_weight, diff_weight)
  }
  dat_diff_weight2<- cbind(diff_crude, dat_diff_weight)
  return(dat_diff_weight2)
}
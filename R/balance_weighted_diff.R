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


##########################################################################################
##
##    Balance Function: Calculates weighted standardized difference for each covariate 
##
##########################################################################################


#' Helper function used within ’balance_weighted_diff()’ to calculate weighted standardized differences after PS adjustment
#'
#' @param data A dataset or matrix containing baseline covariates
#' @param data0 A dataset or matrix containing baseline covariates for unexposed group
#' @param data1 A dataset or matrix containing baseline covariates for exposed group
#' @param score A dataset or matrix of fitted propensity score values (each column corresponds to predicted values from a different model) 
#' @param treatment A vector of binary indicators indicating treatment status
#' @param method weighting method used to calculate weighted standardized differences
#' @param normalized boolean TRUE/FALSE to indicate use of normalized weights (default is TRUE)
#' @returns A numeric vector containing the standardized differences for each covariate after PS adjustment.
#' @details The weighted_diff() function is used within the balance_weighted_diff() function to calculate the adjusted standardized differences for each covariate. 
#' @export
#' @examples
#' #load library
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
#' #running treatment_model() function
#' trt_out<- treatment_model(data=Xcovs_sim, treatment=e, foldid=foldid, alpha=1,lambda_ratio=.01, maxit=5000, nmodels=9)
#' 
#' #running weighted_diff function
#' ps_dat_crossfit<- trt_out[[2]][,1]
#' data0<- Xcovs_sim[e==0,]
#' data1<- Xcovs_sim[e==1,]
#' weighted_diff(data=Xcovs_sim, data0=data0, data1=data1, score=ps_dat_crossfit, treatment=e, method='ow', normalized=TRUE)
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

<!-- README.md is generated from README.Rmd. Please edit that file -->

# PSLassoSynthNC

<!-- badges: start -->

<!-- badges: end -->

## Background

The least absolute shrinkage and selection operator (LASSO) has become
the most widely used tool for fitting large-scale propensity score (PS)
models. LASSO uses $L_1$ regularized regression to prevent overfitting
by shrinking coefficients toward zero (setting some exactly to zero).
The degree of regularization is typically selected using
cross-validation to minimize out-of-sample prediction error. Both theory
and simulations have shown, however, that when using LASSO models for PS
weighting, less regularization is needed to minimize bias in PS weighted
estimators. This is referred to as undersmoothing the LASSO model, where
the optimal degree of undersmoothing can be derived from the target
causal parameter’s efficient influence function. In many settings,
however, the efficient influence function is unknown or difficult to
derive.

In recent work by Wyss et al. (2025), the authors consider the use of
balance metrics as a simple and generally applicable approach to select
the degree of undersmoothing when the efficient influence function is
unknown. Because LASSO models that are tuned using balance metrics alone
are not assured to minimize bias in PS weighted estimators—as such
metrics are blind to the efficient influence function—Wyss et al further
propose a framework to generate synthetic negative control exposures for
bias detection. This package provides R code that uses balance criteria
to undersmooth LASSO PS-weighted estimators, and the use of synthetic
negative control exposures for bias detection.

This package is based on work in the following paper:

- Wyss R, Hansen BB, Hahn G, van der Laan L, Lin KJ. Undersmoothed LASSO
  models for propensity score weighting and synthetic negative control
  exposures for bias detection. *arXiv:2506.17760\[stat.ME\]*. 2025.

Additional references on undersmoothing propensity score models are
provided in the following:

- Ju C, Wyss R, Franklin JM, et al. Collaborative-controlled LASSO for
  constructing propensity score-based estimators in high-dimensional
  data. *Statistical Methods in Medical Research*. 2019;28(4):1044-1063.

- Ertefaie A, Hejazi NS, van der Laan MJ. Nonparametric
  inverse-probability-weighted estimators based on the highly adaptive
  lasso. *Biometrics*. 2023;79(2):1029-1041.

- Wyss R, van der Laan M, Gruber S, et al. Targeted learning with an
  undersmoothed Lasso propensity score model for large-scale covariate
  adjustment in healthcare database studies. *American Journal of
  Epidemiology*. 2024. Online ahead of print.

## Description

This package provides functions to fit undersmoothed Lasso models for
propensity score weighting. The degree of undersmoothing is determined
by balance criteria applied to the propensity score weighted cohorts.
Balance criteria for selecting the regularization tuning parameter
include:

1)  *Largest standardized absolute mean difference*: Selects the
    penalization tuning parameter value that minimizes the maximal
    absolute value of standardized differences after PS weighting.
2)  *Average standardized absolute mean difference*: Selects the
    penalization tuning parameter value that minimizes the L1 norm of
    standardized differences after PS weighting.

Options for propensity score weighting include

1)  Inverse probability weighting
2)  Overlap weighting
3)  Matching weights

The package further provides diagnostics for evaluating model
performance. In addition to balance and prediction diagnostics, we
propose using synthetically generated negative control exposures for
bias detection. The package provides a function for the process of
generating synthetic negative control cohorts.

## Overview of Files in R Folder

A brief description of each R file in the package is provided below. A
more detailed description for each function (including a description of
input parameters) is provided in by rendering documentation for the
specific function through the help function in R (i.e.,
*help(name_of_function)*).

- **treatment_model.R**:
  - contains a function called *treatment_model()* that calls glmnet to
    fit several Lasso models for treatment assignment with varying
    degrees of undersmoothing and returns a dataframe of the out-of-fold
    predictions for each fitted lasso model.
- **balance_weighted_diff.R**:
  - contains a function called *balance_weighted_diff()* that calculates
    weighted standardized differences across treatment groups (called by
    models after PS weighting
- **ps_undersmooth_bal.R**:
  - contains a function called *ps_undersmooth_bal()* that takes as
    input a matrix of fitted PS values and selects the degree of
    undersmoothing using selected balance criteria.
- **utils.R**:
  - contains helper functions that are used within the simulation code.
    Helper functions include functions to calculate prediction
    diagonstics (auc and negative log-likelihood) along with a function
    to create rando folds for cross-validation that can be stratified
    evenly across a variable (e.g., exposure).
- **dat_gen.R**:
  - contains a helper function to generate data used for the simulation
    study in the paper by Wyss et al. (2025).
- **hal_model.R**:
  - simply calls the HAL function from the package hal9001 to generate
    the design matrix representing the expanded set of binary indicator
    functions of the covariates.
  - Coyle J, Hejazi N, Phillips R, van der Laan L, van der Laan M
    (2025). *hal9001: The scalable highly adaptive lasso*.
    <doi:10.5281/zenodo.3558313>
    <https://doi.org/10.5281/zenodo.3558313>, R package version 0.4.6,
    <https://github.com/tlverse/hal9001>.
  - Hejazi N, Coyle J, van der Laan M (2020). “hal9001: Scalable highly
    adaptive lasso regression in R.” *Journal of Open Source Software*.
    <doi:10.21105/joss.02526> <https://doi.org/10.21105/joss.02526>,
    <https://doi.org/10.21105/joss.02526>.

## Installation

You can install the development version of PSLassoSynthNC from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("richiewyss/PSLassoSynthNC")
```

## Example

Below, we provide an example for running the full analytic pipeline from
PS estimation, PS model selection/undersmoothing, treatment effect
estimation, and diagnostic assessment for one simulation run. This
pipeline was used in the publication by Wyss et al. (2025).

Note: The code below is just for one simulation run. When running a
simulation study (as in Wyss et al. (2025), the code needs to be run in
a 'for' loop (for i in 1:nsim) to get multiple simulation runs. To
reduce computation time for simulation studies, each simluation run
could be run in parallel with minor edits (e.g., using a parallel
processing package).

### R Code for Running Analysis on One Generated Cohort (one simulation run)

``` r
library(PSLassoSynthNC)
library(devtools)
devtools::load_all()

###########################
## Generate data

library(hal9001)
library(Matrix)
library(glmnet)
library(dplyr)

scenario<- 1
seed1<- 100       ## seed for simulating random data (when running in loop, should be different for each run)
seed2<- 110101    ## seed for setting global parameters for scenario 1 (should stay the same across runs)

dat<- dat_gen(n=5000, scenario=scenario, seed1=seed1, seed2=seed2)

############################
## Creating Folds
nfolds = 10  
cvfolds <- stratifyCVFoldsByYandID(V=nfolds, Y = dat$A)
folds <- cvfolds$validRows
foldid <- cvfolds$fold_id

#################################################################################################################
## Getting baseline covariates (for Scenario 2, requires fitting HAL on full data to extract design matrix)
x_names<- names(dat)[which(substr(names(dat), 1, 1)=='x')]
X_dat<- dat[,names(dat) %in% x_names]

if(scenario == 1){ x_basis_full<- X_dat }
if(scenario == 2){ x_basis_full<- hal_model(X=X_dat, Y=dat$A, max_degree=2, nfolds=NULL, foldid=foldid, num_knots=c(100, 25)) } 

## Reducing dimension of basis matrix (cleaning based on prevalence)
temp1<- apply(x_basis_full, 2, mean)
temp2<- (temp1 >= .01 & temp1 <= .99)
      
x_basis_full<- x_basis_full[,temp2]

################################################################################################
## Fitting glmnet on x_basis from above and getting predicted values for multiple lambda values
lasso_object<- treatment_model(data=x_basis_full, 
                               treatment=dat$A,
                               foldid=foldid,
                               alpha=1,
                               lambda_ratio=0.01, #ifelse(nrow(x_basis_full) < ncol(x_basis_full), 0.01, 1e-04)
                               nlambda=200,
                               nmodels=NULL,
                               maxit=5000,
                               penalty=NULL,
                               par=FALSE)

preds<- NULL
preds_cv<- NULL
preds_all<- NULL
coef_mat<- NULL
lambdas<- NULL
nvars<- NULL
      
preds<- lasso_object[[1]]       #same-sample predictions
preds_cv<- lasso_object[[2]]    #cross-validated (out-of-fold) predictions
coef_mat<- lasso_object[[3]]    #coefficient matrix (coefficients for each lambda)
lambdas<- lasso_object[[4]]     #lambda values
nvars<- lasso_object[[5]]       #number of variables selected by each lambda

#select a range of lambda values to reduce computation time instead of selecting all values
steps<- NULL
steps<- floor(quantile(1:ncol(preds), seq(0, 1, .02))) 
    
preds_all<- preds_cv[,steps]
coef_mat<- coef_mat[,steps]
lambdas<- lambdas[steps]
nvars<- nvars[steps]

######################################
## Prediction Diagnostics
cstat<- NULL
nll<- NULL
cstat<- apply(preds_all, 2, function(x) auc(dat$A, x)) #auc
nll<- apply(preds_all, 2, function(x) nloglik(dat$A, x)) #negative log-likelihood

##########################################################
## Undersmoothing using balance diagnostics
bal_results1<- NULL
bal_results2<- NULL
bal_results3<- NULL

## undersmoothing based on iptw balance
bal_results1<- ps_undersmooth_bal(data=x_basis_full,
                                  treatment=dat$A,
                                  ps_dat=preds_all,
                                  method=’iptw’,
                                  normalized=TRUE,
                                  standardize=TRUE)
                                  
## undersmoothing based on matching weight balance
bal_results2<- ps_undersmooth_bal(data=x_basis_full,
                                  treatment=dat$A,
                                  ps_dat=preds_all,
                                  method=’mw’,
                                  normalized=TRUE,
                                  standardize=TRUE)
                                  
## undersmoothing based on overlap weight balance
bal_results3<- ps_undersmooth_bal(data=x_basis_full,
                                  treatment=dat$A,
                                  ps_dat=preds_all,
                                  method=’ow’,
                                  normalized=TRUE,
                                  standardize=TRUE)
                                  
## cell with minimum ASAMD
bal_avg_counts1<- bal_results1$index[1]
bal_avg_counts2<- bal_results2$index[1]
bal_avg_counts3<- bal_results3$index[1]

## cell with smallest maximum standardized difference
bal_max_counts1<- bal_results1$index[2]
bal_max_counts2<- bal_results2$index[2]
bal_max_counts3<- bal_results3$index[2]

## asamd value
bal_avg1<- bal_results1$balance[1]
bal_avg2<- bal_results2$balance[1]
bal_avg3<- bal_results3$balance[1]

## value of maximum standardized difference
bal_max1<- bal_results1$balance[2]
bal_max2<- bal_results2$balance[2]
bal_max3<- bal_results3$balance[2]

bal_counts_avg<- c(bal_avg_counts1,
                   bal_avg_counts2,
                   bal_avg_counts3)
                   
bal_counts_max<- c(bal_max_counts1,
                   bal_max_counts2,
                   bal_max_counts3)
                   
bal_avg<- c(bal_avg1,
            bal_avg2,
            bal_avg3)
            
bal_max<- c(bal_max1,
            bal_max2,
            bal_max3)
            
names(bal_counts_avg)<- paste0(’count’, 1:length(bal_counts_avg))
names(bal_counts_max)<- paste0(’count’, 1:length(bal_counts_max))
names(bal_avg)<- paste0(’avg’, 1:length(bal_avg))
names(bal_max)<- paste0(’max’, 1:length(bal_max))


########################################
## Estimating Treatment Effects

## unadjusted estimate
est_crude <- NULL
est_crude <- mean(dat$Y[dat$A==1]) - mean(dat$Y[dat$A==0])

## adjusted estimates
est_ipw<- est_mw<- est_ow<- NULL
for(l in 1:ncol(preds_all)){
  A_hat<- preds_all[,l]
  
  ## IPW Estimation using Hajek estimator (normalized average, same as MLE)
  weight<- NULL
  weight<- dat$A * 1/A_hat + (1-dat$A) * 1/(1-A_hat)
  y1_ipw_est <- sum(dat$A * dat$Y * weight) / sum(dat$A * weight)
  y0_ipw_est <- sum((1-dat$A) * dat$Y * weight) / sum((1-dat$A) * weight)
  est_ipw[l] <- y1_ipw_est - y0_ipw_est
  
  ## Matching Weights using Hajek estimator (normalized average, same as MLE)
  score_data<- as.data.frame(cbind(A_hat, (1-A_hat)))
  numerator<- apply(score_data, 1, min)
  weight<- NULL
  weight<- numerator / (dat$A*A_hat + (1-dat$A)*(1-A_hat))
  y1_mw_est<- sum(dat$A * dat$Y * weight) / sum(dat$A * weight)
  y0_mw_est<- sum((1-dat$A) * dat$Y * weight) / sum((1-dat$A) * weight)
  est_mw[l]<- y1_mw_est - y0_mw_est
  
  ## Overlap Weights using Hajek estimator (normalized average, same as MLE)
  weight<- NULL
  weight<- dat$A * (1-A_hat) + (1-dat$A) * A_hat
  y1_ow_est<- sum(dat$A * dat$Y * weight) / sum(dat$A * weight)
  y0_ow_est<- sum((1-dat$A) * dat$Y * weight) / sum((1-dat$A) * weight)
  est_ow[l]<- y1_ow_est - y0_ow_est
}

## Writing estimates to data frame
est_ipw_all<- as.data.frame(rbind(est_ipw))
est_mw_all<- as.data.frame(rbind(est_mw))
est_ow_all<- as.data.frame(rbind(est_ow))
colnames(est_ipw_all)<- paste0(’ipw_est’, 1:ncol(est_ipw_all))
colnames(est_mw_all)<- paste0(’mw_est’, 1:ncol(est_mw_all))
colnames(est_ow_all)<- paste0(’ow_est’, 1:ncol(est_ow_all))

## effect estimate for ipw, mw, and ow from cross-validated Lasso model
est_ipw_all[1]
est_mw_all[1]
est_ow_all[1]

## effect estimate for ipw, mw, and ow from undersmoothed Lasso model minimizing ASAMD
est_ipw_all[bal_avg_counts1]
est_mw_all[bal_avg_counts2]
est_ow_all[bal_avg_counts3]

## effect estimate for ipw, mw and ow from undersmoothed Lasso model with smallest maximum standardized difference
est_ipw_all[bal_max_counts1]
est_mw_all[bal_max_counts2]
est_ow_all[bal_max_counts3]
```

### Generating Synthetic Negative Control Exposures and Running Analysis on Synthetic Cohorts

``` r

##predicted values from CV Lasso fitted to full data (model that minimizes CV prediction error)
psS<- preds_all[,1] 

## creating sampling probabilities
Xcovs_unexp <- x_basis_full[dat$A==0,]                   # subsetting cohort to unexposed
prop_unexp <- psS[dat$A==0]                              # subsetting predicted PS data to unexposed
theta<- log(prop_unexp/(1-prop_unexp))                   # log odds of the propensity score
exposureRate<- mean(dat$A)                               # prevalence of exposure in the full population
fn <- function(c) mean(plogis(c + theta)) - exposureRate # function to find intercept (i.e., finds value for c so that function is 0)
delta <- uniroot(fn, lower = -100, upper = 100)$root     # delta is intercept value
pi<- plogis(delta + theta)                               # pi are selection probabilities
y0<- dat$Y[dat$A==0]                                     # outcome in unexposed

## Loop to generate and run analyses on many synthetic negative control exposure cohorts
nruns<- 5 #number of synthetic datasets to generate and run analyses on
est_ipw_synth_cv<- est_mw_synth_cv<- est_ow_synth_cv<- NULL #objects to store estimates from CV model
est_ipw_synth_b1<- est_mw_synth_b1<- est_ow_synth_b1<- NULL #objects to store estimates from undersmoothed model 1
est_ipw_synth_b2<- est_mw_synth_b2<- est_ow_synth_b2<- NULL #objects to store estimates from undersmoothed model 2

for(j in 1:nruns){
  ## bootstrap oversampling and assigning synthetic exposure
  sample_index<- sample(1:nrow(Xcovs_unexp), nrow(x_basis_full), replace=TRUE) # sampling with replacement
  Xcovs_boot<- Xcovs_unexp[sample_index,]                                      # sample_index is the index for sampled individuals
  y0_boot<- y0[sample_index]                                                   # outcomes corresponding to sampled individuals
  pi_boot<- pi[sample_index]                                                   # probabilities for synthetic exposure
  zz<- rbinom(dim(Xcovs_boot)[1], 1, pi_boot)                                  # assignment of synthetic exposure status
  
  ## cleaning data to remove sparse variables
  temp1<- NULL
  temp2<- NULL
  temp1<- apply(Xcovs_boot, 2, mean)
  temp2<- (temp1 >= .01 & temp1 <= .99)
  Xcovs_boot<- Xcovs_boot[,temp2]
  
  ## refitting Lasso PS models in pseudo-population
  cvfolds<- NULL
  folds<- NULL
  foldid<- NULL
  
  cvfolds <- stratifyCVFoldsByYandID(V=nfolds, Y = zz)
  folds <- cvfolds$validRows
  foldid <- cvfolds$fold_id
  
  Lasso_object<- treatment_model(data=Xcovs_boot,
                                 treatment=zz,
                                 foldid=foldid,
                                 alpha=1,
                                 lambda_ratio=0.01, #ifelse(nrow(x_basis_full) < ncol(x_basis_full), 0.01, 1e-04)
                                 nlambda=200,
                                 nmodels=NULL,
                                 maxit=5000,
                                 penalty=NULL,
                                 par=FALSE)
                                 
  preds<- NULL
  preds_cv<- NULL
  preds_all<- NULL
  coef_mat<- NULL
  lambdas<- NULL
  nvars<- NULL
  
  preds<- Lasso_object[[1]]    #same-sample predictions
  preds_cv<- Lasso_object[[2]] #cross-validated (out-of-fold) predictions
  coef_mat<- Lasso_object[[3]] #coefficient matrix (coefficients for each lambda)
  lambdas<- Lasso_object[[4]]  #lambda values
  nvars<- Lasso_object[[5]]    #number of variables selected by each lambda
  
  #index values for subset of predictions and lambdas (select a range to reduce computation time instead of selecting all values)
  steps<- NULL
  steps<- floor(quantile(1:ncol(preds), seq(0, 1, .02)))
  preds_all<- preds_cv[,steps]
  coef_mat<- coef_mat[,steps]
  lambdas<- lambdas[steps]
  nvars<- nvars[steps]
  
  ## prediction diagnostics
  cstat<- NULL
  nll<- NULL
  cstat<- apply(preds_all, 2, function(x) auc(zz, x))
  nll<- apply(preds_all, 2, function(x) nloglik(zz, x))
  
  ## Undersmoothing using balance diagnostics
  bal_results1<- NULL
  bal_results2<- NULL
  bal_results3<- NULL
  bal_results1<- ps_undersmooth_bal(data=Xcovs_boot,
                                    treatment=zz,
                                    ps_dat=preds_all,
                                    method=’iptw’,
                                    normalized=TRUE,
                                    standardize=TRUE)
                                    
  bal_results2<- ps_undersmooth_bal(data=Xcovs_boot,
                                    treatment=zz,
                                    ps_dat=preds_all,
                                    method=’mw’,
                                    normalized=TRUE,
                                    standardize=TRUE)
                                    
  bal_results3<- ps_undersmooth_bal(data=Xcovs_boot,
                                    treatment=zz,
                                    ps_dat=preds_all,
                                    method=’ow’,
                                    normalized=TRUE,
                                    standardize=TRUE)
                                    
  ## cell with minimum ASAMD
  bal_avg_counts1<- bal_results1$index[1]
  bal_avg_counts2<- bal_results2$index[1]
  bal_avg_counts3<- bal_results3$index[1]
  
  ## cell with smallest maximum standardized difference
  bal_max_counts1<- bal_results1$index[2]
  bal_max_counts2<- bal_results2$index[2]
  bal_max_counts3<- bal_results3$index[2]
  
  ## asamd
  bal_avg1<- bal_results1$balance[1]
  bal_avg2<- bal_results2$balance[1]
  bal_avg3<- bal_results3$balance[1]
  
  ## maximum standardized difference
  bal_max1<- bal_results1$balance[2]
  bal_max2<- bal_results2$balance[2]
  bal_max3<- bal_results3$balance[2]
  
  bal_counts_avg<- c(bal_avg_counts1,
                     bal_avg_counts2,
                     bal_avg_counts3)
                     
  bal_counts_max<- c(bal_max_counts1,
                     bal_max_counts2,
                     bal_max_counts3)
                     
  bal_avg<- c(bal_avg1,
              bal_avg2,
              bal_avg3)
              
  bal_max<- c(bal_max1,
              bal_max2,
              bal_max3)
  
  names(bal_counts_avg)<- paste0(’count’, 1:length(bal_counts_avg))
  names(bal_counts_max)<- paste0(’count’, 1:length(bal_counts_max))
  names(bal_avg)<- paste0(’avg’, 1:length(bal_avg))
  names(bal_max)<- paste0(’max’, 1:length(bal_max))
  
  ############################################
  ## Estimating Synthetic Exposure Effects
  
  ## unadjusted synthetic estimate
  est_crude_synth<- NULL
  est_crude_synth<- mean(y0_boot[zz==1]) - mean(y0_boot[zz==0])
  
  ## adjusted synthetic estimates
  est_ipw_synth<- est_mw_synth<- est_ow_synth<- NULL
  for(l in 1:ncol(preds_all)){
    zz_hat<- preds_all[,l]
    
    ## IPW Estimation using Hajek estimator (normalized average, same as MLE)
    weight<- NULL
    weight<- zz * 1/zz_hat + (1-zz) * 1/(1-zz_hat)
    y1_ipw_synth <- sum(zz * y0_boot* weight) / sum(zz * weight)
    y0_ipw_synth <- sum((1-zz) * y0_boot * weight) / sum((1-zz) * weight)
    est_ipw_synth[l] <- y1_ipw_synth - y0_ipw_synth
    
    ## Matching Weights using Hajek estimator (normalized average, same as MLE)
    score_data<- as.data.frame(cbind(zz_hat, (1-zz_hat)))
    numerator<- apply(score_data, 1, min)
    weight<- NULL
    weight<- numerator / (zz*zz_hat + (1-zz)*(1-zz_hat))
    y1_mw_synth<- sum(zz * y0_boot * weight) / sum(zz * weight)
    y0_mw_synth<- sum((1-zz) * y0_boot * weight) / sum((1-zz) * weight)
    est_mw_synth[l]<- y1_mw_synth - y0_mw_synth
    
    ## Overlap Weights using Hajek estimator (normalized average, same as MLE)
    weight<- NULL
    weight<- zz * (1-zz_hat) + (1-zz) * zz_hat
    y1_ow_synth<- sum(zz * y0_boot * weight) / sum(zz * weight)
    y0_ow_synth<- sum((1-zz) * y0_boot * weight) / sum((1-zz) * weight)
    est_ow_synth[l]<- y1_ow_synth - y0_ow_synth
  }
  
  ## Writing estimates to data frame
  est_ipw_synth_all<- as.data.frame(rbind(est_ipw_synth))
  est_mw_synth_all<- as.data.frame(rbind(est_mw_synth))
  est_ow_synth_all<- as.data.frame(rbind(est_ow_synth))
  colnames(est_ipw_synth_all)<- paste0(’ipw_synth’, 1:ncol(est_ipw_synth_all))
  colnames(est_mw_synth_all)<- paste0(’mw_synth’, 1:ncol(est_mw_synth_all))
  colnames(est_ow_synth_all)<- paste0(’ow_synth’, 1:ncol(est_ow_synth_all))
  
  ## appending results for synthetic effect estimate for ipw, mw, and ow from cross-validated Lasso model
  est_ipw_synth_cv <- c(est_ipw_synth_cv, est_ipw_synth_all[1])
  est_mw_synth_cv <- c(est_mw_synth_cv, est_mw_synth_all[1])
  est_ow_synth_cv <- c(est_ow_synth_cv, est_ow_synth_all[1])
  
  ## appending results for synthetic effect estimate for ipw, mw and ow from undersmoothed Lasso model minimizing ASAMD
  est_ipw_synth_b1 <- c(est_ipw_synth_b1, est_ipw_synth_all[bal_avg_counts1])
  est_mw_synth_b1 <- c(est_mw_synth_b1, est_mw_synth_all[bal_avg_counts2])
  est_ow_synth_b1 <- c(est_ow_synth_b1, est_ow_synth_all[bal_avg_counts3])
  
  ## appending results for synthetic effect estimate from undersmoothed Lasso model with smallest maximum standardized difference
  est_ipw_synth_b2 <- c(est_ipw_synth_b2, est_ipw_synth_all[bal_max_counts1])
  est_mw_synth_b2 <- c(est_mw_synth_b2, est_mw_synth_all[bal_max_counts2])
  est_ow_synth_b2 <- c(est_ow_synth_b2, est_ow_synth_all[bal_max_counts3])
}

## synthetic effect estimate for CV estimators averaged across all synthetic datasets
mean(as.numeric(est_ipw_synth_cv))
mean(as.numeric(est_mw_synth_cv))
mean(as.numeric(est_ow_synth_cv))

## synthetic effect estimate for undersmoothed estimators minimizing ASAMD averaged across all synthetic datasets
mean(as.numeric(est_ipw_synth_b1))
mean(as.numeric(est_mw_synth_b1))
mean(as.numeric(est_ow_synth_b1))

## synthetic effect estimate for undersmoothed estimators with smallest maximum st diff averaged across all synthetic datasets
mean(as.numeric(est_ipw_synth_b2))
mean(as.numeric(est_mw_synth_b2))
mean(as.numeric(est_ow_synth_b2))
```

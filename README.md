
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PSLassoSynthNC

<!-- badges: start -->

<!-- badges: end -->

# Background

This package fits undersmoothed LASSO models for propensity score
weighting. When fitting Lasso propensity score models, both theory and
simulations have shown that undersmoothing the Lasso model is beneficial
for reducint bias in estimated treatment effects. However, determining
the degree of undersmoothing can be challenging. Too much undersmoothing
can cause severe overfitting resulting in poorly fit PS models. However,
too little undersmoothing will result in fitted propensity score models
that exclude important confounder information.

This package uses balance diagnostics for undersmoothing
high-dimensional Lasso propensity score models. The package provides
diagnostics for evaluating model performance. In addition to balance and
prediciton diagnostics, we propose using synthetically generated
negative control exposures for bias detection. The package provides a
function to automate the process of synthetic negative control exposure
generation.

Options for propensity score weighting include 1) inverse probability
weighting; 2) overlap weighting; and 3) matching weights.

This package is based on work in the following paper:

- Wyss R, Hansen BB, Hahn G, van der Laan L, Lin KJ.
  Collaborative-controlled LASSO for constructing propensity score-based
  estimators in high-dimensional data. *arXiv:2506.17760\[stat.ME\]*.
  2025.

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

# Description

The function **main_fun()** runs the full analytic pipeline from PS
estimation, PS model selection/undersmoothing, treatment effect
estimation, and diagnostic assessment.

Within *main_fun()*, there are several helper functions that are called.
Each of these functions can be called outside of *main_fun()* to give
users flexibility to use specific procedures outside of *main_fun()*.
For example, if users are interested in fitting undersmoothed lasso PS
models, but want to apply the fitted propensity scores in a causal
inference framework that is not provided by *main_fun()*, users can call
the functions outside of *main_fun()* to run selected parts of the
analytic pipeline.

The list of helper functions within *main_fun()* are described below.

# Helper functions called by **main_fun()**

- **treatment_model()**:
  - calls glmnet to fit several Lasso models for treatment assignment
    with varying degrees of undersmoothing and returns a dataframe of
    propensity score values for each fitted lasso model.
- **outcome_model()**:
  - calls glmnet to fit an outcome lasso model that is tuned using
    cross-validation and returns a vector of the predicted values along
    with variables selected by the model to be used when running the
    outcome adaptive lasso for propensity score estimation.  
- **ps_dist_plot()**:
  - plots propensity score distributions for all fitted lasso propensity
    score models across treatment groups and calculates c-statistics.
- **balance_weighted_diff()**:
  - calculates weighted standardized differences across treatment groups
    (called by *cov_diff_plot()*).
- **cov_diff_plot()**:
  - plots standardized differences across treatment groups for all of
    the fitted lasso PS models after PS weighting or PS matching.
- **ps_undersmooth()**:
  - takes as input a matrix of fitted PS values and selects the degree
    of undersmoothing using selected criteria.
  - Options include:
    - collaborative learning
    - score equation
    - balance minimization
- **ps_weighting()**:
  - Uses PS weighting to estimate the treatment effect
  - Options include:
    - inverse probability weights (target population ATE)
    - standardized mortality ratio weighting (target population ATT)
    - overlap weights (depends on overlap in estimated PS)
    - matching weights (depends on overlap in estimated PS)
- **matchit()**:
  - calls the *matchit()* function within the *MatchIt* package for PS
    matching.
- **tmle()**:
  - calls the *tmle()* function within the *TMLE* package when using
    targeted minimum loss-based estimation for estimating treatment
    effects.

## Installation

You can install the development version of PSLassoSynthNC from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("richiewyss/PSLassoSynthNC")
```

## Example

This is a basic example for running **main()**. Examples for running
each of the helper functions listed above are provided in separate
vignettes.

``` r
library(PSLassoSynthNC)
library(devtools)
devtools::load_all()
```

## Simulate some example data

``` r
nstudy<- 2000
nvars<- 500

nc<- 100
ns<- nvars-(nc)

alpha_temp<- runif(nc, 0.0, 0.4)
beta_temp<- runif(nc, 0.0, 0.4)

random_neg<- sample(1:length(alpha_temp), 0.5*length(alpha_temp), replace=FALSE)
alpha_temp[random_neg]<- -1*alpha_temp[random_neg]
beta_temp[random_neg]<-  -1*beta_temp[random_neg]

alpha<-  matrix(c(alpha_temp, rep(0, ns)), ncol=1)
beta<-   matrix(c(beta_temp, rep(0, ns)), ncol=1)

betaE<- 0

cprev<- runif(nvars, 0, 0.3)
cprev<- sample(cprev)
oprev<- 0.05
tprev<- 0.4
    
Xcovs_sim<- matrix(rnorm((nstudy*nvars), 0, 1), nrow=nstudy, ncol=nvars)
Xcovs_sim<- as.data.frame(Xcovs_sim)  
names(Xcovs_sim)<- c(paste0('x', 1:nvars))
      
W<- as.matrix(Xcovs_sim)
colnames(W)<- c(paste0('x', 1:nvars))

linear_pred_e<- W %*% alpha
linear_pred_y<- W %*% beta
      
treatment_inc<- tprev
fn <- function(c) mean(plogis(c + linear_pred_e)) - treatment_inc
alpha0 <- uniroot(fn, lower = -20, upper = 20)$root
Ee <- (1 + exp( -(alpha0 + linear_pred_e) ))^-1
e<- rbinom(nstudy, 1, Ee)
      
outcome_inc<- oprev
fn <- function(c) mean(plogis(c + betaE*e + linear_pred_y  )) - outcome_inc
beta0 <- uniroot(fn, lower = -20, upper = 20)$root
Ey <- (1 + exp( -( beta0 + betaE*e + linear_pred_y )))^-1
y<- rbinom(nstudy, 1, Ey)

simdat<- NULL
simdat <- as.data.frame(cbind(y, e, Ee, Xcovs_sim))
```

## Running **main_fun()**

The function **main_fun()** takes as input:

- *outcome*:
  - a numeric vector representing the outcome.
- *treatment*:
  - a vector of binary indicators representing treatment.
- *Xcovs*:
  - a dataframe or matrix of baseline covariates.
- *model_type*:
  - the type of lasso PS model. Options include ‘lasso’ or ‘oal’ for
    outcome adaptive lasso.
- *undersmooth*:
  - the method for selecting the degree of undersmoothing. Options
    include ‘CV’, ‘balance’, ‘ctmle’, or ‘score’.
- *nmodels*
  - the number of lasso models to fit and consider for degrees of
    undersmoothing (must be between 1 and 12).
- *cross_fitting*
  - whether or not to use cross-fitted PS values (default is TRUE).
- *nfolds*
  - number of CV folds to use when fitting lasso models (default is 10).
- *trunc_type*
  - optional to truncate PS values or PS weights (default is NULL).
    Options include ‘symmetric’, ‘assymetric’, ‘psvalue’, ‘psweight’.
- *trunc_val_low*
  - if ‘trunc_type’ is not null, then need to specify the lower value on
    which to truncate PS (between 0 and 1)
- *trunc_val_high*
  - if ‘trunc_type’ is not null, then need to specify the upper value on
    which to truncate PS (between 0 and 1)
- *method*
  - the causal inference method used for estimating treatment effects.
    Options include ‘ow’, ‘mw’, ‘ipw’, ‘smrw’, ‘matching’, or ‘tmle’.

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.

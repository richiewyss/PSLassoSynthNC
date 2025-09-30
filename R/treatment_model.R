##################################################################################
##
## treatment_model: function to to fit lasso model for many values of lambda
##
##################################################################################
#’ Calls glmnet to fit propensity score models corresponding to different degrees of regularization
#’
#’ @param data a dataset or matrix containing baseline covariates
#’ @param treatment a binary vector for treatment
#’ @param foldid fold each subject belongs to
#’ @param alpha the elasticnet tuning parameter defined in glmnet (default is 1 for Lasso)
#’ @param lambda_ratio the ratio from the largest to smallest lambda to consider (constrains the range of lambda values)
#’ @param nlambda the number of lambda tuning parameters to consider when fitting glmnet
#’ @param maxit maximum number of iterations as defined in glmnet
#’ @param nmodels the number of undersmoothed models to return PS values for
#’ @param penalty optional vector to specify different penalties for variables like in adaptive lasso (default is NULL)
#’ @param par optional TRUE/FALSE to implement parallel computing (default is FALSE)
treatment_model<- function(data,
                           treatment,
                           foldid,
                           alpha=1,
                           lambda_ratio=NULL,
                           nlambda=100,
                           nmodels=NULL,
                           maxit=5000,
                           penalty=NULL,
                           par=FALSE){
  Wmat = as.matrix(data)
  sx = Matrix::Matrix(Wmat, sparse=TRUE)
  if(is.null(penalty)){
    penalty_factor<- rep(1, ncol(Wmat))
  }
  if(!is.null(penalty)){
    penalty_factor<- penalty
  }
  if(is.null(lambda_ratio)){
    lambda_ratio<- ifelse(nrow(sx) < ncol(sx), 0.01, 1e-04)
  }
  glmnet.e<- NULL
  glmnet.e<- glmnet::cv.glmnet(x = sx,
                               y = treatment,
                               family = "binomial",
                               type.measure = "deviance",
                               alpha = alpha, ##ridge regression alpha=0, lasso alpha=1
                               nlambda = nlambda,
                               lambda.min.ratio = lambda_ratio, #ifelse(nrow(sx) < ncol(sx), 0.01, 1e-04),
                               #nfolds = 10, ##don’t specify this when using foldid
                               #penalty.factor = c(rep(1, dim(sx)[2])),
                               penalty.factor = penalty_factor,
                               parallel = par,
                               standardize = FALSE, ###HAL does not standardize design matrix
                               maxit=maxit, #5000,
                               foldid = foldid,
                               keep=TRUE)
  
  ## extracting predicted values & lambda starting and stopping values
  gns1<- NULL
  gns2<- NULL
  gns1<- predict(glmnet.e, newx = Wmat, s=glmnet.e$lambda, type = "response") ##predicted values for each lambda
  gns2<- plogis(glmnet.e$fit.preval[, 1:ncol(gns1)]) ## Cross Validated (out-of-fold) predicted values for each lambda
  
  ## lambda value that optimizes CV prediction (minimizes CV prediction error)
  n.lambda.start<- which(glmnet.e$lambda == glmnet.e$lambda.min)
  lambda.start<- glmnet.e$lambda.min
  
  ## only keeping lambda values that are less than or equal to lambda that optimizes CV prediction
  preds_under1<- gns1[,n.lambda.start:ncol(gns1)] ##same-sample predicted values
  preds_under2<- gns2[,n.lambda.start:ncol(gns2)] ##out-of-fold predicted values
  lambda_vector<- glmnet.e$lambda[n.lambda.start:length(glmnet.e$lambda)]
  lassocoef = glmnet.e$glmnet.fit$beta[,n.lambda.start:length(glmnet.e$lambda)]
  coef_mat<- lassocoef
  
  ##number of selected variables for each lambda value
  n_selected_vars<- apply(lassocoef, 2, function(x){sum(x != 0)})
  
  results<- list(preds_under1,
                 preds_under2,
                 coef_mat,
                 lambda_vector,
                 n_selected_vars)
  
  names(results)<- c("preds",
                     "preds_cf",
                     "coef_mat",
                     "lambdas",
                     "n_vars")
  return(results)
}
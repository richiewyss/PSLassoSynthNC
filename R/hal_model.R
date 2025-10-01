###################################################################################################
##
## Function to fit HAL to generate matrix of indicator basis functions (used for Scenario 2) 
##      see hal9001 package in R for details
##
###################################################################################################

#' @export
hal_model<- function(X, Y, max_degree, num_knots, nfolds, foldid){
  # fitting HAL
  mod_full<- fit_hal(X=X,
                     Y=Y,
                     X_unpenalized = NULL,
                     max_degree = max_degree,
                     smoothness_orders = 0,
                     num_knots = num_knots,
                     reduce_basis = 0.01,
                     family = c("binomial"),
                     lambda = 10,
                     id = NULL,
                     offset = NULL,
                     fit_control = list(cv_select = FALSE, nfolds = nfolds, foldid = foldid, use_min = TRUE,
                                        lambda.min.ratio = .001, prediction_bounds = "default"),
                     basis_list = NULL,
                     return_lasso = TRUE,
                     return_x_basis = TRUE,
                     yolo = FALSE)
  ## design matrix
  x_basis_full<- mod_full$x_basis[,-1] ## first column is intercept term (all ones)
  return(x_basis_full)
}

#############################################################################################
##
## helper functions to calculate prediciton diagnostics: NLL and Cstat
##      1) nloglik(): calculates the negative log-likelihood for binary response variable
##      2) auc(): calculates the c-statistic for binary response variable
##
#############################################################################################

#' @export
nloglik <- function(y, pred, trunc = 0.001) {
  pred <- pmin(pmax(pred, trunc), 1-trunc)
  if (all(y == round(y))) {
    - mean(ifelse(y==1, log(pred), log(1-pred)))
  } else {
    - mean(y * log(pred) + (1-y) * log(1-pred))
  }
}

#' @export
auc <- function(y, pred) {
  require("ROCR")
  performance(prediction(pred, y), "auc")@y.values[[1]]
}




###################################################################################################################
##
## helper function to create folds stratified by variable: created by Susan Gruber (used in Wyss et al. 2024)
##
###################################################################################################################

#' @export
stratifyCVFoldsByYandID <- function (V, Y, id = NULL) {
  # 1. distribute the ids that have Y = 1 in any of the rows equally among all the folds,
  # 2. separately, distribute the ids that have Y = 0 for all rows equally among the folds
  # ensure that V is less than the number of cases
  if (is.null(id)) id <- 1:length(Y)
  case_status_by_id <- by(Y, id, sum) # this gives n.unique results, sorted by id #
  case_ids <- names(case_status_by_id)[ case_status_by_id > 0]
  noncase_ids <- names(case_status_by_id)[ case_status_by_id == 0]
  if (V > min(length(case_ids), length(noncase_ids))) {
    stop("number of observations in minority class is less than the number of folds")
  }
  valSet.case_ids <- split(sample(case_ids), rep(1:V, length = length(case_ids)))
  valSet.noncase_ids <- split(sample(noncase_ids), rep(1:V, length = length(noncase_ids)))
  validRows <- vector("list", length = V)
  names(validRows) <- paste(seq(V))
  fold_id <- rep(NA, length(Y))
  for (v in seq(V)){
    validRows[[v]] <- which(as.character(id) %in% c(valSet.case_ids[[v]], valSet.noncase_ids[[v]]))
    fold_id[validRows[[v]]] <- v
  }
  return(list(validRows = validRows, fold_id = fold_id))
}
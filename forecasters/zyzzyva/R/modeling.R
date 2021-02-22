#' @include calibration.R transformation.R
NULL

#' Fit a learner to the training data
#'
#' @param train_test list of training and testing data
#' @param modeling_options a named list, additional elements of which
#'     overrides learner-dependent options set within the code
ml.fit_model <- function(train_test,
                         modeling_options) {
  if (modeling_options$learner == "linear") {
    ml.linear(train_test, modeling_options)
  } else {
    stop("available learners are ('linear')")
  }
}

ml.extract_slope_and_lag_columns <- function(X, base_covariates, location_covariates) {
  # create new_train_X, new_test_X with slope variables
  slope_base_vars <- paste(sapply(base_covariates, function(x) x$name),
                           "lag", sep = "_")
  new_X <- matrix(NA, nrow=nrow(X), ncol=length(slope_base_vars))
  for (i in 1:length(slope_base_vars)) {
    var_cols <- which(startsWith(colnames(X), slope_base_vars[i]))
    n_var_cols <- length(var_cols)
    if (n_var_cols <= 1) {
      stop("Could not create slope vars in training X, not enough columns for ", var)
    }
    new_X[,i] <- stats::lsfit(-(1:n_var_cols), t(X[, var_cols]))$coef[2,]
  }
  # add first lags of handpicked variables
  location_covariate_names <- sapply(location_covariates, function(x) x$name)
  first_cases <- sapply(c(slope_base_vars, location_covariate_names),
                        function(var) which(startsWith(colnames(X), var))[1],
                        USE.NAMES = F)
  new_X <- cbind(new_X, X[, first_cases])
  
  covariate_names <- c(paste(slope_base_vars, "slope", sep="_"),
                       colnames(X)[first_cases])
  
  colnames(new_X) <- covariate_names
  return(new_X)
}


#' Fit a simple linear model.
#'
#' The following model is fit:
#'    response ~ sum(lag(c) for c in covariates) + log(pop) +
#'            sum(slope(c) for c in covariates).
#' The slope is computed using least squares over past lags of each signal,
#' for each location individually.
#' Calibration is not currently performed.
#' 
#' @param train_test list of training and testing data
#' @param modeling_options a named list, additional elements of which
#'     overrides learner-dependent options set within the code
#' @importFrom logger log_debug
#' @importFrom Iso pava
#' @importFrom stats lsfit
#' @importFrom quantgen quantile_lasso
ml.linear <- function(train_test,
                      modeling_options) {
  # extract items from train_test
  train_X <- train_test$train_X
  train_y <- train_test$train_y
  test_X <- train_test$test_X
  train_row_locations <- train_test$train_row_locations[[1]]
  test_row_locations <- train_test$test_row_locations[[1]]

  # create new_train_X, new_test_X with slope variables
  new_train_X <- ml.extract_slope_and_lag_columns(train_X,
                                                  modeling_options$base_covariates,
                                                  modeling_options$location_covariates)
  new_test_X <- ml.extract_slope_and_lag_columns(test_X,
                                                 modeling_options$base_covariates,
                                                 modeling_options$location_covariates)
  logger::log_debug("Input colnames:\n",
                    paste(colnames(new_train_X), collapse='\n'))
  
  # finally, get predictions
  out <- matrix(NA, nrow=nrow(test_X), ncol=length(modeling_options$quantiles))
  # get quantiles
  fit_quantiles <- quantgen::quantile_lasso(new_train_X,
                                            train_y,
                                            tau = modeling_options$quantiles,
                                            standardize = FALSE,
                                            lambda = 0)
  preds <-  predict(fit_quantiles, new_test_X, sort = TRUE)

  if (FALSE) {
  # if (modeling_options$recalibrate) {
    # isotonic regression to ensure quantiles ordering
    for (i in 1:nrow(preds)) {
      preds[i, ] <- Iso::pava(preds[i, ])
    }

    tau_calibrated <- ca.get_empirical_quantiles(modeling_options$quantiles,
                                                 preds,
                                                 train_y)
    fit_quantiles <- quantgen::quantile_lasso(new_train_X,
                                          train_y,
                                          tau = tau_calibrated,
                                          standardize = FALSE,
                                          lambda = 0)
    preds <- predict(fit_quantiles, new_test_X, sort = TRUE)
  }

  colnames(preds) <- modeling_options$quantiles
  out <- preds
  
  # isotonic regression to ensure quantiles ordering
  for (i in 1:nrow(out)) {
    out[i, ] <- Iso::pava(out[i, ])
  }

  # transform back to correct scale for clipping
  if (modeling_options$log_response) {
    out <- tr.inv_log_pad(out)
  }
  
  # ensure negatives are set to 0
  out[out < 0] <- 0
  
  # format for output
  final_out <- NULL
  for (i in 1:length(modeling_options$quantiles)) {
    final_out <- rbind(final_out, cbind(modeling_options$ahead,
                                        train_test$test_row_locations, 
                                        modeling_options$quantiles[i],
                                        out[, i]))
  }
  names(final_out) <- c("ahead", "geo_value", "quantile", "value")
  
  final_out
}

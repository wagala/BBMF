#' Compute Binary Matrix Factorization Reconstruction Metrics
#'
#' Computes reconstruction-performance measures comparing an observed or
#' reference binary matrix with a reconstructed binary matrix.
#'
#' @param X A binary reference matrix containing values in \code{0, 1}.
#' @param X_hat A binary reconstructed matrix with the same dimensions as
#'   \code{X}.
#'
#' @return A named numeric vector containing:
#' \describe{
#'   \item{\code{Specificity}}{
#'     The true-negative rate,
#'     \eqn{TN / (TN + FP)}.
#'   }
#'   \item{\code{F1}}{
#'     The harmonic mean of precision and sensitivity.
#'   }
#'   \item{\code{MCC}}{
#'     The Matthews correlation coefficient.
#'   }
#'   \item{\code{Recon_Error_Rate}}{
#'     The proportion of matrix entries incorrectly reconstructed,
#'     \eqn{(FN + FP) / (KG)}.
#'   }
#' }
#'
#' @details
#' The function first computes the confusion-matrix counts:
#' true positives (TP), false negatives (FN), false positives (FP),
#' and true negatives (TN).
#'
#' Sensitivity, precision, specificity, F1-score, F2-score, and the
#' Matthews correlation coefficient are calculated internally.
#' The current return value includes specificity, F1-score, MCC, and
#' reconstruction error rate.
#'
#' Undefined metrics are returned as \code{NA}.
#'
#' The confusion counts are coerced to double precision before computing
#' the Matthews correlation coefficient to reduce the risk of integer
#' overflow for large matrices.
#'
#' @examples
#' X <- matrix(
#'   c(1, 0, 1, 0,
#'     1, 1, 0, 0),
#'   nrow = 2,
#'   byrow = TRUE
#' )
#'
#' X_hat <- matrix(
#'   c(1, 0, 0, 0,
#'     1, 1, 1, 0),
#'   nrow = 2,
#'   byrow = TRUE
#' )
#'
#' bmf_metrics(X, X_hat)
#'
#' @export
#'
bmf_metrics <- function(X, X_hat) {
  
  if (!all(dim(X) == dim(X_hat))) {
    stop("X and X_hat must have the same dimensions.")
  }
  
  # Confusion counts
  TP <- as.numeric(
    sum(X == 1 & X_hat == 1)
  )
  
  FN <- as.numeric(
    sum(X == 1 & X_hat == 0)
  )
  
  FP <- as.numeric(
    sum(X == 0 & X_hat == 1)
  )
  
  TN <- as.numeric(
    sum(X == 0 & X_hat == 0)
  )
  
  # Sensitivity (recall)
  sensitivity <- ifelse(
    (TP + FN) == 0,
    NA,
    TP / (TP + FN)
  )
  
  # Precision
  precision <- ifelse(
    (TP + FP) == 0,
    NA,
    TP / (TP + FP)
  )
  
  # Specificity
  specificity <- ifelse(
    (TN + FP) == 0,
    NA,
    TN / (TN + FP)
  )
  
  # F1-score
  f1 <- ifelse(
    is.na(precision) |
      is.na(sensitivity) |
      (precision + sensitivity) == 0,
    NA,
    2 * precision * sensitivity /
      (precision + sensitivity)
  )
  
  # F2-score
  f2 <- ifelse(
    is.na(precision) |
      is.na(sensitivity) |
      (4 * precision + sensitivity) == 0,
    NA,
    5 * precision * sensitivity /
      (4 * precision + sensitivity)
  )
  
  # Matthews correlation coefficient
  denom_terms <- c(
    TP + FP,
    TP + FN,
    TN + FP,
    TN + FN
  )
  
  if (any(denom_terms == 0)) {
    
    mcc <- NA
    
  } else {
    
    mcc <- (
      TP * TN - FP * FN
    ) / sqrt(
      prod(denom_terms)
    )
  }
  
  # Reconstruction error
  recon_error <- FN + FP
  
  recon_error_rate <- recon_error / (
    nrow(X) * ncol(X)
  )
  
  return(
    c(
      Specificity = specificity,
      F1 = f1,
      MCC = mcc,
      Recon_Error_Rate = recon_error_rate
    )
  )
}
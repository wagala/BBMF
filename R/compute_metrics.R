#' Reconstruction Accuracy Metrics
#'
#' Computes specificity, F1 score, Matthews Correlation Coefficient (MCC),
#' and reconstruction error rate for an estimated reconstruction against a
#' known ground-truth binary matrix.
#'
#' @param Z_est A binary (0/1) matrix of dimension \eqn{K \times G}: the
#'   estimated reconstruction (e.g. \code{bool_prod(W_hat, H_hat)} from a
#'   fitted \code{\link{run_BBMF}}, or the reconstruction produced by
#'   a competing method such as Asso or GreConD+).
#' @param X_truth A binary (0/1) matrix of the same dimension as
#'   \code{Z_est}: the true, noise-free matrix used to generate the
#'   simulated data.
#'
#' @return A one-row \code{data.frame} with columns \code{Specificity},
#'   \code{F1}, \code{MCC}, and \code{ReconError}.
#'
#' @details
#' Let \eqn{TP}, \eqn{TN}, \eqn{FP}, \eqn{FN} count entries where
#' \code{Z_est} and \code{X_truth} are, respectively, both 1, both 0,
#' 1-but-truth-0, and 0-but-truth-1. Then
#' \deqn{\mathrm{Specificity} = \frac{TN}{TN+FP}, \qquad
#'  \mathrm{F1} = \frac{2\,\mathrm{precision}\cdot\mathrm{recall}}
#'  {\mathrm{precision}+\mathrm{recall}},}
#' \deqn{\mathrm{MCC} = \frac{TP\cdot TN - FP\cdot FN}
#'  {\sqrt{(TP+FP)(TP+FN)(TN+FP)(TN+FN)}}, \qquad
#'  \mathrm{ReconError} = \frac{FP+FN}{TP+TN+FP+FN}.}
#' \code{MCC} is returned as \code{NA} if its denominator is zero (e.g.\
#' when \code{Z_est} or \code{X_truth} is constant).
#'
#' \strong{Numerical note:} \code{sum()} applied to a logical vector
#' returns an R integer, and at realistic problem sizes (\eqn{K \times G}
#' in the thousands or more) the product of four such counts in the MCC
#' denominator can exceed the 32-bit integer limit
#' (\code{.Machine$integer.max}), silently returning \code{NA} via an
#' integer-overflow warning rather than a numeric answer. Every count is
#' therefore coerced to \code{double} with \code{as.numeric()} before
#' multiplication, which is not optional at scale.
#'
#' @examples
#' \dontrun{
#' metrics_table <- rbind(
#'   data.frame(Method = "Asso",     compute_metrics(X_asso,       X_true)),
#'   data.frame(Method = "GreConD+", compute_metrics(X_GrenC_plus, X_true)),
#'   data.frame(Method = "BBMF",     compute_metrics(Z_hat,        X_true))
#' )
#' print(metrics_table)
#' }
#'
#' @export
compute_metrics <- function(Z_est, X_truth) {

  TP <- sum(Z_est == 1 & X_truth == 1)
  TN <- sum(Z_est == 0 & X_truth == 0)
  FP <- sum(Z_est == 1 & X_truth == 0)
  FN <- sum(Z_est == 0 & X_truth == 1)

  specificity <- TN / (TN + FP)
  precision   <- TP / (TP + FP)
  recall      <- TP / (TP + FN)
  F1          <- 2 * precision * recall / (precision + recall)

  # Coerce to double BEFORE multiplying: sum() of a logical vector is an
  # R integer, and (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) can overflow the 32-bit
  # integer range well before any of the four factors looks large on its
  # own -- confirmed to occur in practice at K=70, G=44 (K*G=3080).
  mcc_num <- as.numeric(TP) * as.numeric(TN) - as.numeric(FP) * as.numeric(FN)
  mcc_den <- sqrt(as.numeric(TP + FP) * as.numeric(TP + FN) *
                  as.numeric(TN + FP) * as.numeric(TN + FN))
  MCC <- ifelse(mcc_den == 0, NA, mcc_num / mcc_den)

  recon_error <- (FP + FN) / (TP + TN + FP + FN)

  data.frame(Specificity = specificity, F1 = F1, MCC = MCC, ReconError = recon_error)
}

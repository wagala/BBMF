#' Add Asymmetric Noise to a Binary Matrix
#'
#' Introduces false-negative and false-positive errors into a binary
#' ground-truth matrix. Entries equal to 1 are independently changed to 0
#' with probability `fnr`, while entries equal to 0 are independently changed
#' to 1 with probability `fpr`.
#'
#' The corresponding BBMF observation probabilities are
#' \eqn{p_{11} = 1 - \mathrm{fnr}} and \eqn{p_{10} = \mathrm{fpr}}.
#'
#' @param X_truth A binary matrix containing only 0 and 1. This represents
#'   the noise-free Boolean matrix.
#' @param fnr Numeric scalar in `[0, 1]`. The false-negative rate,
#'   \eqn{P(X = 0 \mid Z = 1)}.
#' @param fpr Numeric scalar in `[0, 1]`. The false-positive rate,
#'   \eqn{P(X = 1 \mid Z = 0)}.
#' @param seed Optional integer random seed for reproducibility. If `NULL`,
#'   the current random-number-generator state is used.
#'
#' @return A list containing:
#' \describe{
#'   \item{X_sim}{The noisy binary matrix.}
#'   \item{fnr}{The requested false-negative rate.}
#'   \item{fpr}{The requested false-positive rate.}
#'   \item{p11}{The implied true-positive probability, `1 - fnr`.}
#'   \item{p10}{The implied false-positive probability, `fpr`.}
#'   \item{realized_fnr}{The observed proportion of true 1s changed to 0.}
#'   \item{realized_fpr}{The observed proportion of true 0s changed to 1.}
#'   \item{realized_noise}{The overall proportion of matrix entries changed.}
#' }
#'
#' @examples
#' X_truth <- matrix(c(1, 0, 1, 0, 1, 0), nrow = 2)
#'
#' noisy <- add_asymmetric_noise(
#'   X_truth,
#'   fnr = 0.10,
#'   fpr = 0.05,
#'   seed = 1
#' )
#'
#' noisy$X_sim
#' noisy$realized_noise
#'
#' @export
add_asymmetric_noise <- function(X_truth, fnr, fpr, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  stopifnot(
    all(X_truth %in% c(0, 1)),
    fnr >= 0, fnr <= 1,
    fpr >= 0, fpr <= 1
  )
  
  X_sim <- X_truth
  
  # Positions of true ones and zeros
  idx1 <- which(X_truth == 1)
  idx0 <- which(X_truth == 0)
  
  # 1 -> 0 flips: false negatives
  flip1 <- stats::runif(length(idx1)) < fnr
  
  # 0 -> 1 flips: false positives
  flip0 <- stats::runif(length(idx0)) < fpr
  
  X_sim[idx1[flip1]] <- 0L
  X_sim[idx0[flip0]] <- 1L
  
  list(
    X_sim = X_sim,
    fnr = fnr,
    fpr = fpr,
    p11 = 1 - fnr,
    p10 = fpr,
    realized_fnr =
      if (length(idx1) > 0) mean(X_sim[idx1] == 0) else NA_real_,
    realized_fpr =
      if (length(idx0) > 0) mean(X_sim[idx0] == 1) else NA_real_,
    realized_noise = mean(X_sim != X_truth)
  )
}
#' Update Rows of W by Block Sampling
#'
#' Performs a block Gibbs update of the rows of the binary matrix
#' \code{W}. For each row, the function enumerates all possible binary
#' configurations across the \eqn{R} latent factors, evaluates their
#' log-posterior probabilities, and samples one configuration from the
#' resulting conditional distribution.
#'
#' @param W A binary matrix of dimension \eqn{K \times R} containing
#'   the current factor activations.
#' @param H A binary matrix containing the current factor-specific
#'   feature patterns.
#' @param X An observed binary data matrix.
#' @param alpha Prior probability parameter(s) associated with the entries
#'   of \code{W}.
#' @param p11 Probability of observing an entry equal to 1 when the
#'   corresponding reconstructed entry is 1.
#' @param p10 Probability of observing an entry equal to 1 when the
#'   corresponding reconstructed entry is 0.
#'
#' @return An updated binary matrix \code{W} with the same dimensions
#'   as the input matrix.
#'
#' @details
#' For each row \eqn{k} of \code{W}, the function constructs all
#' \eqn{2^R} possible binary configurations. Each candidate configuration
#' is temporarily assigned to row \eqn{k}, and its log-posterior
#' probability is evaluated using \code{compute_log_postW()}.
#'
#' The resulting log-posterior values are normalized using the
#' log-sum-exp calculation implemented in \code{sumLog()}. One candidate
#' configuration is then sampled using \code{stats::rmultinom()} and assigned
#' to \code{W[k, ]}.
#'
#' Rows are updated sequentially, so each subsequent row update
#' conditions on the most recently updated values of \code{W}.
#'
update_W_Block_rows <- function(W, H, X, alpha, p11, p10) {
  K <- nrow(W)  # Number of rows
  R <- ncol(W)  # Number of columns
  
  for (k in 1:K) {
    # Generate all possible configurations for the current row
    configurations <- expand.grid(rep(list(c(0, 1)), R))
    
    log_posteriors <- apply(configurations, 1, function(row_config) {
      W[k, ] <- as.numeric(row_config)
      compute_log_postW(W, H, X, alpha, p11, p10)
    })
    
    # Normalize the log posteriors to probabilities
    log_post = sumLog(log_posteriors)
    probs = exp(log_posteriors - log_post)
    
    # Use a multinomial approach to sample a configuration based on computed probabilities
    sampled_config_index <- which(stats::rmultinom(1, 1, probs) == 1)
    W[k, ] <- as.numeric(configurations[sampled_config_index, ])
  }
  
  return(W)
}
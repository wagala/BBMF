#' Update Columns of H by Block Sampling
#'
#' Performs a block Gibbs update of the columns of the binary matrix
#' \code{H}. For each column, the function enumerates all possible binary
#' configurations across the \eqn{R} latent factors, evaluates their
#' log-posterior probabilities, and samples one configuration from the
#' resulting conditional distribution.
#'
#' @param W A binary matrix containing the current factor activations.
#' @param H A binary matrix of dimension \eqn{R \times G} containing
#'   the current factor-specific feature patterns.
#' @param X An observed binary data matrix.
#' @param beta Prior probability parameter(s) associated with the entries
#'   of \code{H}.
#' @param p11 Probability of observing an entry equal to 1 when the
#'   corresponding reconstructed entry is 1.
#' @param p10 Probability of observing an entry equal to 1 when the
#'   corresponding reconstructed entry is 0.
#'
#' @return An updated binary matrix \code{H} with the same dimensions
#'   as the input matrix.
#'
#' @details
#' For each column \eqn{g} of \code{H}, the function constructs all
#' \eqn{2^R} possible binary configurations. Each candidate configuration
#' is temporarily assigned to column \eqn{g}, and its log-posterior
#' probability is evaluated using \code{compute_log_posteriorH()}.
#'
#' The resulting log-posterior values are normalized using the
#' log-sum-exp calculation implemented in \code{sumLog()}. One candidate
#' configuration is then sampled using \code{stats::rmultinom()} and assigned
#' to \code{H[, g]}.
#'
#' Columns are updated sequentially, so each subsequent column update
#' conditions on the most recently updated values of \code{H}.
#'
update_H_Block_columns<- function(W, H, X, beta, p11, p10) {
  R <- nrow(H)  # Number of rows
  G <- ncol(H)  # Number of columns
  
  for (g in 1:G) {
    # Generate all possible configurations for the current column
    configurations <- expand.grid(rep(list(c(0, 1)), R))
    
    log_posteriors <- apply(configurations, 1, function(col_config) {
      H[, g] <- as.numeric(col_config)
      # compute log posterior
      compute_log_posteriorH(W, H, X, beta, p11, p10)
    })
    
    # Convert log posteriors to probabilities using log-sum-exp for stability
    log_post = sumLog(log_posteriors)
    probs = exp(log_posteriors -log_post)
    
    # Use stats::rmultinom to sample a configuration based on computed probabilities
    sampled_config_index <- which(stats::rmultinom(1, 1, probs) == 1)
    H[, g] <- as.numeric(configurations[sampled_config_index, ])
  }
  
  return(H)
}
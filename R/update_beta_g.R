#' Update the Feature-Specific Probabilities Beta
#'
#' Samples the feature-specific probabilities \eqn{\beta_g} from their
#' full conditional Beta distributions given the current binary matrix
#' \code{H} and the sparsity indicators \eqn{\psi_g}.
#'
#' @param H A binary matrix of dimension \eqn{R \times G} containing
#'   the current factor-specific feature patterns.
#' @param psi_g A binary vector of length \eqn{G} containing the
#'   feature-specific sparsity indicators \eqn{\psi_g}.
#' @param b1 First shape parameter of the Beta prior used when
#'   \eqn{\psi_g = 1}.
#' @param b2 Second shape parameter of the Beta prior used when
#'   \eqn{\psi_g = 1}.
#' @param c1 First shape parameter of the Beta prior used when
#'   \eqn{\psi_g = 0}.
#' @param c2 Second shape parameter of the Beta prior used when
#'   \eqn{\psi_g = 0}.
#'
#' @return A numeric vector of length \eqn{G} containing the updated
#'   feature-specific probabilities \eqn{\beta_g}.
#'
#' @details
#' For each feature \eqn{g}, the function computes the number of active
#' entries in column \eqn{g} of \code{H},
#'
#' \deqn{
#' H_g = \sum_{r=1}^{R} H_{rg}.
#' }
#'
#' If \eqn{\psi_g = 1}, the full conditional distribution is
#'
#' \deqn{
#' \beta_g \mid - \sim
#' \mathrm{Beta}(H_g + b_1,\; R - H_g + b_2).
#' }
#'
#' If \eqn{\psi_g = 0}, the full conditional distribution is
#'
#' \deqn{
#' \beta_g \mid - \sim
#' \mathrm{Beta}(H_g + c_1,\; R - H_g + c_2).
#' }
#'
update_beta_g <- function(H, psi_g, b1, b2, c1, c2) {
  R <- nrow(H)  # Number of observations per group (number of rows)
  G <- ncol(H)  # Number of groups (number of columns)
  
  # Initialize a vector to store sampled beta_g values
  beta_g_samples <- numeric(G)
  
  # Calculate the number of successes for each group g
  H_g <- colSums(H)
  
  # Iterate over each group g
  for (g in 1:G) {
    if (psi_g[g] == 1) {
      # Sample from Beta(H_g + b1, R - H_g + b2) when psi_g = 1
      beta_g_samples[g] <- stats::rbeta(1, H_g[g] + b1, R - H_g[g] + b2)
    } else if (psi_g[g] == 0) {
      # Sample from Beta(H_g + c1, R - H_g + c2) when psi_g = 0
      beta_g_samples[g] <- stats::rbeta(1, H_g[g] + c1, R - H_g[g] + c2)
    } else {
      stop("psi_g must be either 0 or 1.")
    }
  }
  
  return(beta_g_samples)
}
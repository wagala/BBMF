#' Update the Global Sparsity Probability Pi
#'
#' Samples the global sparsity probability \code{pi} from its full
#' conditional Beta distribution given the current binary sparsity
#' indicators \eqn{\psi_g}.
#'
#' @param psi_g A binary vector of length \eqn{G} containing the
#'   feature-specific sparsity indicators.
#' @param d1 First shape parameter of the Beta prior for \code{pi}.
#' @param d2 Second shape parameter of the Beta prior for \code{pi}.
#'
#' @return A numeric scalar containing the updated value of \code{pi}.
#'
#' @details
#' Let \eqn{\boldsymbol{\psi} = (\psi_1, \ldots, \psi_G)} denote the
#' vector of binary sparsity indicators. Assuming
#'
#' \deqn{
#' \pi \sim \mathrm{Beta}(d_1, d_2),
#' }
#'
#' the full conditional distribution is
#'
#' \deqn{
#' \pi \mid \boldsymbol{\psi}
#' \sim
#' \mathrm{Beta}\left(
#' d_1 + \sum_{g=1}^{G}\psi_g,
#' d_2 + G - \sum_{g=1}^{G}\psi_g
#' \right).
#' }
#'
update_pi <- function(psi_g, d1, d2) {
  # Calculate the total number of successes (psi_g = 1) and failures (psi_g = 0)
  sum_psi <- sum(psi_g)  # Total number where psi_g = 1
  G <- length(psi_g)     # Total number of psi_g values
  
  # Calculate the parameters for the Beta distribution
  alpha_post <- sum_psi + d1
  beta_post <- G - sum_psi + d2
  
  # Sample pi from the Beta distribution
  pi_s <- stats::rbeta(1, alpha_post, beta_post)
  
  return(pi_s)
}
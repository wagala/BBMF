#' Log-Prior Density for the Factor-Pattern Matrix H
#'
#' Computes the log-density of the factor-pattern matrix \code{H}
#' under its Bernoulli prior, with success probability
#' \eqn{\beta_g \phi_r} for entry \eqn{H_{rg}}.
#'
#' @param H A binary matrix of dimension \eqn{R \times G} containing
#'   the current factor-specific feature patterns.
#' @param beta A numeric vector of length \eqn{G} containing the
#'   per-feature activation probabilities \eqn{\beta_g}.
#' @param phi A numeric vector of length \eqn{R} containing the
#'   per-factor inclusion probabilities \eqn{\phi_r}.
#'
#' @return A numeric scalar containing the log-prior density of
#'   \code{H}.
#'
#' @details
#' Each entry is assumed independent with
#' \eqn{H_{rg} \mid \beta_g, \phi_r \sim \mathrm{Bernoulli}(\beta_g \phi_r)}.
#'
#' @export
compute_log_prior_H <- function(H, beta, phi) {

  R <- nrow(H)
  G <- ncol(H)

  p_mat <- outer(phi, beta)  # entry [r,g] = phi_r * beta_g

  logprior <- sum(H * log(p_mat) + (1 - H) * log(1 - p_mat))

  return(logprior)
}

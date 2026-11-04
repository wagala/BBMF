#' Log-Prior Density for the Factor-Inclusion Probabilities phi
#'
#' Computes the log-density of the factor-level inclusion
#' probabilities \code{phi} under independent Beta priors.
#'
#' @param phi A numeric vector of length \eqn{R} containing the
#'   per-factor inclusion probabilities \eqn{\phi_r}.
#' @param s1 Shape-1 hyperparameter(s) of the Beta prior on
#'   \code{phi}: either a scalar (shared across all \eqn{r}) or a
#'   numeric vector of length \eqn{R}.
#' @param s2 Shape-2 hyperparameter(s) of the Beta prior on
#'   \code{phi}: either a scalar or a numeric vector of length
#'   \eqn{R}.
#'
#' @return A numeric scalar containing the log-prior density of
#'   \code{phi}.
#'
#' @details
#' Each \eqn{\phi_r \sim \mathrm{Beta}(s_{1,r}, s_{2,r})}
#' independently. The Beta normalizing constant is omitted, since
#' this is a single (non-mixture) prior compared against nothing else
#' in any downstream full conditional -- the omitted constant is an
#' overall additive term that never affects any decision made from
#' this density.
#'
#' A numerically safe convention is used at the boundary: whenever a
#' shape parameter equals 1, the corresponding exponent is exactly
#' zero and the term is treated as identically zero, even if
#' \code{phi} lies at 0 or 1 (where the naive computation
#' \code{0 * log(0)} would otherwise evaluate to \code{NaN} rather
#' than its correct limiting value of zero).
#'
#' @export
compute_log_prior_phi <- function(phi, s1, s2) {

  R <- length(phi)
  if (length(s1) == 1) s1 <- rep(s1, R)
  if (length(s2) == 1) s2 <- rep(s2, R)

  term1 <- ifelse(s1 == 1, 0, (s1 - 1) * log(phi))
  term2 <- ifelse(s2 == 1, 0, (s2 - 1) * log(1 - phi))

  logprior <- sum(term1 + term2)

  return(logprior)
}

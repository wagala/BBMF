#' Sample the Factor-Inclusion Probabilities phi from Their Prior
#'
#' Draws a fresh vector of factor-level inclusion probabilities
#' \code{phi} from independent Beta priors.
#'
#' @param s1 Shape-1 hyperparameter(s) of the Beta prior: either a
#'   scalar (shared across all \eqn{r}) or a numeric vector of length
#'   \eqn{R}.
#' @param s2 Shape-2 hyperparameter(s) of the Beta prior: either a
#'   scalar or a numeric vector of length \eqn{R}.
#' @param R The number of latent factors.
#'
#' @return A numeric vector of length \eqn{R}, with each entry drawn
#'   from \eqn{\mathrm{Beta}(s_{1,r}, s_{2,r})}.
#'
#' @export
sample_prior_phi <- function(s1, s2, R) {

  if (length(s1) == 1) s1 <- rep(s1, R)
  if (length(s2) == 1) s2 <- rep(s2, R)

  phi <- stats::rbeta(R, s1, s2)

  return(phi)
}

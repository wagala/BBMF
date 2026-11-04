#' Update the Slab-Prior Probability pi (Full Conditional)
#'
#' Draws a new value of \code{pi}, the prior probability that a
#' feature belongs to the slab component, from its full conditional
#' distribution.
#'
#' @param psi A binary vector of length \eqn{G}, the current
#'   sparsity indicators.
#' @param d1 Shape-1 hyperparameter of the Beta prior on \code{pi}.
#' @param d2 Shape-2 hyperparameter of the Beta prior on \code{pi}.
#'
#' @return A numeric scalar, the updated \code{pi}.
#'
#' @details
#' An ordinary Beta-Bernoulli conjugate update,
#' \deqn{\pi \mid \mathrm{rest} \sim \mathrm{Beta}\!\left(d_1 +
#'  \textstyle\sum_g \psi_g,\ \ d_2 + G - \sum_g \psi_g\right).}
#'
#' @export
update_pi <- function(psi, d1, d2) {

  G <- length(psi)
  pi_g <- stats::rbeta(1, d1 + sum(psi), d2 + G - sum(psi))

  return(pi_g)
}

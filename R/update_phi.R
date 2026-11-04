#' Update the Factor-Inclusion Probabilities phi (Full Conditional)
#'
#' Draws a new vector of factor-level inclusion probabilities
#' \code{phi} from its full conditional distribution, using the
#' auxiliary variables \code{B} and \code{D}.
#'
#' @param B A binary matrix of dimension \eqn{K \times R}, from
#'   \code{\link{update_AB}}.
#' @param D A binary matrix of dimension \eqn{R \times G}, from
#'   \code{\link{update_CD}}.
#' @param s1 Shape-1 hyperparameter(s) of the Beta prior on
#'   \code{phi}: scalar or length-\eqn{R} vector.
#' @param s2 Shape-2 hyperparameter(s) of the Beta prior on
#'   \code{phi}: scalar or length-\eqn{R} vector.
#'
#' @return A numeric vector of length \eqn{R}, the updated \code{phi}.
#'
#' @details
#' Because \eqn{B_{kr} \mid \phi_r \sim \mathrm{Bernoulli}(\phi_r)}
#' for every sample \eqn{k}, and \eqn{D_{rg} \mid \phi_r \sim
#' \mathrm{Bernoulli}(\phi_r)} for every feature \eqn{g}, pooling all
#' \eqn{K+G} of these auxiliary draws with the Beta prior gives an
#' ordinary conjugate update,
#' \deqn{\phi_r \mid \mathrm{rest} \sim \mathrm{Beta}\!\left(s_{1,r} +
#'  \textstyle\sum_k B_{kr} + \sum_g D_{rg},\ \ s_{2,r} + (K -
#'  \textstyle\sum_k B_{kr}) + (G - \sum_g D_{rg})\right).}
#' This is the mechanism by which the effective number of active
#' factors can fall below the user-specified upper bound \eqn{R}: a
#' factor with \eqn{\phi_r \to 0} is switched off in both \code{W} and
#' \code{H} simultaneously.
#'
#' @export
update_phi <- function(B, D, s1, s2) {

  K <- nrow(B)
  R <- ncol(B)
  G <- ncol(D)

  if (length(s1) == 1) s1 <- rep(s1, R)
  if (length(s2) == 1) s2 <- rep(s2, R)

  sumB <- colSums(B)
  sumD <- rowSums(D)

  phi <- stats::rbeta(R,
                       s1 + sumB + sumD,
                       s2 + (K - sumB) + (G - sumD))

  return(phi)
}

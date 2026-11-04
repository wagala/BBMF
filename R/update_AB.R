#' Update the Auxiliary Pair (A, B) for the Factor-Activation Matrix
#'
#' Draws the auxiliary matrices \code{A} and \code{B} whose logical
#' AND reproduces the current \code{W}, given \eqn{A_{kr} \mid
#' \alpha_k \sim \mathrm{Bernoulli}(\alpha_k)} and \eqn{B_{kr} \mid
#' \phi_r \sim \mathrm{Bernoulli}(\phi_r)}.
#'
#' @param W A binary matrix of dimension \eqn{K \times R}.
#' @param alpha A numeric vector of length \eqn{K}.
#' @param phi A numeric vector of length \eqn{R}.
#'
#' @return A list with elements \code{A} and \code{B}, each a binary
#'   matrix of dimension \eqn{K \times R}.
#'
#' @details
#' \code{A} is the variable directly governed by \code{alpha} and is
#' used in place of \code{W} when updating \code{alpha} (see
#' \code{\link{update_alpha}}); \code{B} is directly governed by
#' \code{phi} and is pooled with \code{D} (from
#' \code{\link{update_CD}}) when updating \code{phi} (see
#' \code{\link{update_phi}}).
#'
#' @export
update_AB <- function(W, alpha, phi) {

  K <- nrow(W)
  R <- ncol(W)

  Pmat <- matrix(alpha, nrow = K, ncol = R)
  Qmat <- matrix(phi, nrow = K, ncol = R, byrow = TRUE)

  result <- sample_aux_pair(W, Pmat, Qmat)

  return(result)
}

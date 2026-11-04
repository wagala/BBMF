#' Update the Auxiliary Pair (C, D) for the Factor-Pattern Matrix
#'
#' Draws the auxiliary matrices \code{C} and \code{D} whose logical
#' AND reproduces the current \code{H}, given \eqn{C_{rg} \mid
#' \beta_g \sim \mathrm{Bernoulli}(\beta_g)} and \eqn{D_{rg} \mid
#' \phi_r \sim \mathrm{Bernoulli}(\phi_r)}.
#'
#' @param H A binary matrix of dimension \eqn{R \times G}.
#' @param beta A numeric vector of length \eqn{G}.
#' @param phi A numeric vector of length \eqn{R}.
#'
#' @return A list with elements \code{C} and \code{D}, each a binary
#'   matrix of dimension \eqn{R \times G}.
#'
#' @details
#' \code{C} is the variable directly governed by \code{beta} and is
#' used in place of \code{H} when updating \code{beta} (see
#' \code{\link{update_beta_psi}}); \code{D} is directly governed by
#' \code{phi} and is pooled with \code{B} (from
#' \code{\link{update_AB}}) when updating \code{phi} (see
#' \code{\link{update_phi}}).
#'
#' @export
update_CD <- function(H, beta, phi) {

  R <- nrow(H)
  G <- ncol(H)

  Pmat <- matrix(beta, nrow = R, ncol = G, byrow = TRUE)
  Qmat <- matrix(phi, nrow = R, ncol = G)

  result <- sample_aux_pair(H, Pmat, Qmat)
  names(result) <- c("C", "D")

  return(result)
}

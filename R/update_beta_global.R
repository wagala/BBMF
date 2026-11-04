#' Update a Single Pooled Activation Probability (BBMF-Global)
#'
#' Draws a single shared activation probability from every
#' \eqn{(r,g)} pair pooled together, in place of the per-feature
#' \eqn{\beta_g}. Used by the \code{"Global"} and \code{"Sym"}
#' ablation variants (see \code{\link{ablation_flags}}). No
#' spike-and-slab: \code{"Global"} sits below \code{"Flat"} in the
#' nested family, so the slab-only restriction is already in force.
#'
#' @param C A binary matrix (typically \eqn{R \times G}), from
#'   \code{\link{update_CD}} or equal to \code{H} directly if
#'   \code{phi} is fixed at 1.
#' @param b1,b2 Beta hyperparameters (the slab shape parameters).
#'
#' @return A numeric scalar.
#'
#' @seealso \code{\link{update_beta_psi}}, \code{\link{update_beta_flat}},
#'   \code{\link{update_alpha_global}}
#'
#' @export
update_beta_global <- function(C, b1, b2) {
  total <- sum(C)
  n <- length(C)
  stats::rbeta(1, b1 + total, b2 + n - total)
}

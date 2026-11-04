#' Update Per-Feature Activation Rates Without Spike-and-Slab (BBMF-Flat)
#'
#' Draws per-feature activation rates \code{beta}, retaining
#' \eqn{\beta_g}'s per-feature indexing but drawing every entry only
#' from the slab component -- \eqn{\psi_g} is never sampled and plays
#' no role. Used by the \code{"Flat"}, \code{"Global"}, and
#' \code{"Sym"} ablation variants (see \code{\link{ablation_flags}}).
#'
#' @param C A binary matrix of dimension \eqn{R \times G}, from
#'   \code{\link{update_CD}} or equal to \code{H} directly if
#'   \code{phi} is fixed at 1.
#' @param b1,b2 Beta hyperparameters of the slab component.
#'
#' @return A numeric vector of length \eqn{G}.
#'
#' @details
#' This is exactly the \eqn{\psi_g = 1} branch of
#' \code{\link{update_beta_psi}}, extracted on its own:
#' \deqn{\beta_g \mid \mathrm{rest} \sim \mathrm{Beta}(n_g + b_1,\
#'  R - n_g + b_2), \qquad n_g = \textstyle\sum_r C_{rg}.}
#'
#' @seealso \code{\link{update_beta_psi}}, \code{\link{update_beta_global}}
#'
#' @export
update_beta_flat <- function(C, b1, b2) {
  R <- nrow(C)
  G <- ncol(C)
  n_g <- colSums(C)
  stats::rbeta(G, n_g + b1, R - n_g + b2)
}

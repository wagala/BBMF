#' Update a Single Pooled Activation Probability (BBMF-Global)
#'
#' Draws a single shared activation probability from every
#' \eqn{(k,r)} pair pooled together, in place of the per-sample
#' \eqn{\alpha_k}. Used by the \code{"Global"} and \code{"Sym"}
#' ablation variants (see \code{\link{ablation_flags}}).
#'
#' @param A A binary matrix (typically \eqn{K \times R}), from
#'   \code{\link{update_AB}} or equal to \code{W} directly if
#'   \code{phi} is fixed at 1.
#' @param a1,a2 Beta hyperparameters.
#'
#' @return A numeric scalar.
#'
#' @details
#' \deqn{\alpha_0 \mid \mathrm{rest} \sim \mathrm{Beta}\!\left(a_1 +
#'  \textstyle\sum A,\ \ a_2 + n - \sum A\right),}
#' where \eqn{n} is the total number of entries in \code{A}. This is
#' the same Beta-Bernoulli conjugacy as \code{\link{update_alpha}},
#' pooled across the whole matrix rather than row by row.
#'
#' @seealso \code{\link{update_alpha}}, \code{\link{update_beta_global}}
#'
#' @export
update_alpha_global <- function(A, a1, a2) {
  total <- sum(A)
  n <- length(A)
  stats::rbeta(1, a1 + total, a2 + n - total)
}

#' Update the Per-Sample Activation Probabilities alpha (Full Conditional)
#'
#' Draws a new vector of per-sample activation probabilities
#' \code{alpha} from its full conditional distribution, using the
#' auxiliary matrix \code{A}.
#'
#' @param A A binary matrix of dimension \eqn{K \times R}, from
#'   \code{\link{update_AB}} (or equal to \code{W} directly if
#'   \code{phi} is fixed at 1; see Details).
#' @param a1 Shape-1 hyperparameter of the Beta prior on
#'   \code{alpha}.
#' @param a2 Shape-2 hyperparameter of the Beta prior on
#'   \code{alpha}.
#'
#' @return A numeric vector of length \eqn{K}, the updated
#'   \code{alpha}.
#'
#' @details
#' In the augmented model, \code{A} (not \code{W}) is the variable
#' directly governed by \code{alpha}:
#' \eqn{A_{kr} \mid \alpha_k \sim \mathrm{Bernoulli}(\alpha_k)}, giving
#' the conjugate update
#' \deqn{\alpha_k \mid \mathrm{rest} \sim \mathrm{Beta}\!\left(a_1 +
#'  \textstyle\sum_r A_{kr},\ \ a_2 + R - \sum_r A_{kr}\right).}
#' When \code{phi} is fixed at 1 for every factor (e.g. the
#' \code{BBMF-NoPhi} ablation variant), \eqn{A_{kr} = W_{kr}}
#' deterministically, and \code{W} may be passed directly in place of
#' \code{A}.
#'
#' @export
update_alpha <- function(A, a1, a2) {

  K <- nrow(A)
  R <- ncol(A)

  sumA <- rowSums(A)
  alpha <- stats::rbeta(K, a1 + sumA, a2 + R - sumA)

  return(alpha)
}

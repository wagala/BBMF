#' Update the Symmetric Noise Parameter rho (BBMF-Sym)
#'
#' Draws a single match probability \code{rho}, used in place of the
#' two separate noise parameters \code{p11}, \code{p10}, with
#' \eqn{p_{11} = \rho} and \eqn{p_{10} = 1-\rho}. Used by the
#' \code{"Sym"} ablation variant (see \code{\link{ablation_flags}}),
#' the restriction that coincides with OrMachine's noise model.
#'
#' @param X An observed binary matrix of dimension \eqn{K \times G}.
#' @param Z The current reconstructed matrix, \eqn{Z = W \circ H}
#'   (see \code{\link{bool_prod}}).
#' @param e1,e2 Beta hyperparameters for \code{rho}.
#'
#' @return A numeric scalar, \code{rho}, satisfying \code{rho > 0.5}.
#'
#' @details
#' Writing \eqn{n_{\mathrm{match}} = \sum_{k,g} \mathbb{1}(X_{kg} =
#' Z_{kg})} and \eqn{N = KG},
#' \deqn{\rho \mid \mathrm{rest} \sim \mathrm{Beta}(e_1 +
#'  n_{\mathrm{match}},\ e_2 + N - n_{\mathrm{match}})\ \text{truncated
#'  to}\ (\tfrac12, 1).}
#' The truncation to \eqn{(\tfrac12,1)} is the one-parameter
#' specialization of the general identifiability constraint
#' \eqn{p_{11} > p_{10}}, and coincides exactly with OrMachine's own
#' constraint \eqn{\lambda \in \mathbb{R}_+}.
#'
#' @seealso \code{\link{update_noise}}, \code{\link{rtrunc_beta}}
#'
#' @export
update_rho <- function(X, Z, e1, e2) {
  Znum <- Z * 1
  n_match <- sum(Znum * X) + sum((1 - Znum) * (1 - X))
  N <- length(X)
  rtrunc_beta(1, e1 + n_match, e2 + N - n_match, lower = 0.5, upper = 1)
}

#' Update the Noise Parameters p11 and p10 (Full Conditional)
#'
#' Draws new values of the sensitivity \code{p11} and false-positive
#' rate \code{p10} from their full conditional distributions, subject
#' to the identifiability constraint \eqn{p_{11} > p_{10}}.
#'
#' @param X An observed binary matrix of dimension \eqn{K \times G}.
#' @param Z The current reconstructed matrix,
#'   \eqn{Z = W \circ H} (see \code{\link{bool_prod}}), of the
#'   same dimension as \code{X}.
#' @param p11 The current sensitivity (used as the lower truncation
#'   bound for the new \code{p10} draw).
#' @param p10 The current false-positive rate (unused in the
#'   computation but retained for interface symmetry with other
#'   update functions).
#' @param e1,e2 Shape hyperparameters of the Beta prior on
#'   \code{p11}.
#' @param f1,f2 Shape hyperparameters of the Beta prior on
#'   \code{p10}.
#'
#' @return A list with elements \code{p11} and \code{p10}, the updated
#'   noise parameters, satisfying \code{p11 > p10} exactly.
#'
#' @details
#' The likelihood is invariant under the joint relabeling
#' \eqn{(Z, p_{11}, p_{10}) \mapsto (1-Z, p_{10},
#' p_{11})}; imposing \eqn{p_{11} > p_{10}} resolves this symmetry by
#' selecting the state \eqn{Z=1} as the one more likely to
#' produce an observed 1. The constraint is enforced exactly via
#' truncated-Beta draws (\code{\link{rtrunc_beta}}), not by rejection
#' or approximation:
#' \deqn{p_{11} \mid \mathrm{rest} \sim \mathrm{Beta}(n_1+e_1,\,
#'  n-n_1+e_2) \ \text{truncated to}\ (p_{10}, 1),}
#' \deqn{p_{10} \mid \mathrm{rest} \sim \mathrm{Beta}(m_1+f_1,\,
#'  m-m_1+f_2) \ \text{truncated to}\ (0, p_{11}),}
#' where \eqn{n_1, n} count entries with \eqn{Z_{kg}=1}
#' (respectively \eqn{X_{kg}=1} among them), and \eqn{m_1, m} count
#' entries with \eqn{Z_{kg}=0}. \code{p10} is drawn using the
#' newly-updated \code{p11}, in that order.
#'
#' @export
update_noise <- function(X, Z, p11, p10, e1, e2, f1, f2) {

  n1 <- sum(Z * X)
  n  <- sum(Z)
  m1 <- sum((1 - Z) * X)
  m  <- sum(1 - Z)

  p11_new <- rtrunc_beta(1, n1 + e1, n - n1 + e2, lower = p10, upper = 1)
  p10_new <- rtrunc_beta(1, m1 + f1, m - m1 + f2, lower = 0, upper = p11_new)

  return(list(p11 = p11_new, p10 = p10_new))
}

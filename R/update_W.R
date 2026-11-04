#' Update the Factor-Activation Matrix W (Full Conditional)
#'
#' Draws a new factor-activation matrix \code{W} from its full
#' conditional distribution given the data, the current \code{H}, and
#' the current values of every other parameter.
#'
#' @param W A binary matrix of dimension \eqn{K \times R} containing
#'   the current factor activations (updated column by column).
#' @param H A binary matrix of dimension \eqn{R \times G} containing
#'   the current factor-specific feature patterns (held fixed during
#'   this call).
#' @param X An observed binary matrix of dimension \eqn{K \times G}.
#' @param alpha A numeric vector of length \eqn{K} containing the
#'   per-sample activation probabilities \eqn{\alpha_k}.
#' @param phi A numeric vector of length \eqn{R} containing the
#'   per-factor inclusion probabilities \eqn{\phi_r}.
#' @param p11 Sensitivity: \eqn{P(X_{kg}=1 \mid Z_{kg}=1)}.
#' @param p10 False-positive rate: \eqn{P(X_{kg}=1 \mid Z_{kg}=0)}.
#'
#' @return An updated binary matrix \code{W} of the same dimension as
#'   the input.
#'
#' @details
#' For each factor \eqn{r}, the coverage of every entry \eqn{(k,g)} by
#' the remaining factors, \eqn{Z_{kg}^{(-r)}}, is tracked via
#' the coverage matrix \eqn{S = WH}, updated incrementally by a
#' rank-one correction as each column of \code{W} is resampled (rather
#' than recomputed from scratch), giving \eqn{O(KGR)} total cost per
#' sweep. Only entries \eqn{g} with \eqn{H_{rg}=1} and
#' \eqn{Z_{kg}^{(-r)}=0} are informative about \eqn{W_{kr}};
#' for these entries \eqn{Z_{kg} = W_{kr}} exactly, giving a
#' closed-form Bernoulli full conditional with success probability
#' \deqn{\omega_{kr} = \frac{\alpha_k \phi_r \, p_{11}^{m^{(1)}_{kr}}
#'  (1-p_{11})^{m^{(0)}_{kr}}}{\alpha_k \phi_r \, p_{11}^{m^{(1)}_{kr}}
#'  (1-p_{11})^{m^{(0)}_{kr}} + (1-\alpha_k \phi_r) \, p_{10}^{m^{(1)}_{kr}}
#'  (1-p_{10})^{m^{(0)}_{kr}}},}
#' where \eqn{m^{(1)}_{kr}} and \eqn{m^{(0)}_{kr}} count the
#' informative entries with \eqn{X_{kg}=1} and \eqn{X_{kg}=0}
#' respectively.
#'
#' @export
update_W <- function(W, H, X, alpha, phi, p11, p10) {

  K <- nrow(W)
  R <- ncol(W)
  G <- ncol(H)

  S <- W %*% H  # coverage counts, computed once per call

  for (r in seq_len(R)) {

    cross_r <- outer(W[, r], H[r, ])
    not_covered <- (S - cross_r) <= 0

    informative <- not_covered * matrix(H[r, ], nrow = K, ncol = G, byrow = TRUE)
    m1 <- rowSums(informative * X)
    m0 <- rowSums(informative * (1 - X))

    log_num <- log(alpha * phi[r]) + m1 * log(p11) + m0 * log(1 - p11)
    log_den <- log(1 - alpha * phi[r]) + m1 * log(p10) + m0 * log(1 - p10)
    omega <- 1 / (1 + exp(log_den - log_num))

    W_new_r <- stats::rbinom(K, 1, omega)
    S <- S - cross_r + outer(W_new_r, H[r, ])  # incremental update, O(KG)
    W[, r] <- W_new_r
  }

  return(W)
}

#' Update the Factor-Pattern Matrix H (Full Conditional)
#'
#' Draws a new factor-pattern matrix \code{H} from its full
#' conditional distribution given the data, the current \code{W}, and
#' the current values of every other parameter.
#'
#' @param W A binary matrix of dimension \eqn{K \times R} containing
#'   the current factor activations (held fixed during this call; is
#'   normally the just-updated output of \code{\link{update_W}}).
#' @param H A binary matrix of dimension \eqn{R \times G} containing
#'   the current factor-specific feature patterns (updated row by
#'   row).
#' @param X An observed binary matrix of dimension \eqn{K \times G}.
#' @param beta A numeric vector of length \eqn{G} containing the
#'   per-feature activation probabilities \eqn{\beta_g}.
#' @param phi A numeric vector of length \eqn{R} containing the
#'   per-factor inclusion probabilities \eqn{\phi_r}.
#' @param p11 Sensitivity: \eqn{P(X_{kg}=1 \mid Z_{kg}=1)}.
#' @param p10 False-positive rate: \eqn{P(X_{kg}=1 \mid Z_{kg}=0)}.
#'
#' @return An updated binary matrix \code{H} of the same dimension as
#'   the input.
#'
#' @details
#' Exact mirror of \code{\link{update_W}}, with the roles of samples
#' and features exchanged: for each factor \eqn{r}, only patients
#' \eqn{k} with \eqn{W_{kr}=1} and no other factor already covering
#' \eqn{(k,g)} are informative about \eqn{H_{rg}}, giving a
#' closed-form Bernoulli full conditional with success probability
#' \deqn{\eta_{rg} = \frac{\beta_g \phi_r \, p_{11}^{n^{(1)}_{rg}}
#'  (1-p_{11})^{n^{(0)}_{rg}}}{\beta_g \phi_r \, p_{11}^{n^{(1)}_{rg}}
#'  (1-p_{11})^{n^{(0)}_{rg}} + (1-\beta_g \phi_r) \, p_{10}^{n^{(1)}_{rg}}
#'  (1-p_{10})^{n^{(0)}_{rg}}}.}
#' The coverage matrix is maintained incrementally, as in
#' \code{\link{update_W}}, giving \eqn{O(KGR)} total cost per call.
#'
#' @export
update_H <- function(W, H, X, beta, phi, p11, p10) {

  K <- nrow(W)
  R <- ncol(W)
  G <- ncol(H)

  S <- W %*% H

  for (r in seq_len(R)) {

    cross_r <- outer(W[, r], H[r, ])
    not_covered <- (S - cross_r) <= 0

    informative <- not_covered * matrix(W[, r], nrow = K, ncol = G, byrow = FALSE)
    n1 <- colSums(informative * X)
    n0 <- colSums(informative * (1 - X))

    log_num <- log(beta * phi[r]) + n1 * log(p11) + n0 * log(1 - p11)
    log_den <- log(1 - beta * phi[r]) + n1 * log(p10) + n0 * log(1 - p10)
    eta <- 1 / (1 + exp(log_den - log_num))

    H_new_r <- stats::rbinom(G, 1, eta)
    S <- S - cross_r + outer(W[, r], H_new_r)
    H[r, ] <- H_new_r
  }

  return(H)
}

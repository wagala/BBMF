#' Update the Spike-and-Slab Indicator psi and Activation Rates beta
#'
#' Draws a new sparsity indicator \code{psi} and per-feature
#' activation rate \code{beta} from their full conditional
#' distributions, using the auxiliary matrix \code{C}.
#'
#' @param C A binary matrix of dimension \eqn{R \times G}, from
#'   \code{\link{update_CD}} (or equal to \code{H} directly if
#'   \code{phi} is fixed at 1).
#' @param beta A numeric vector of length \eqn{G}, the current
#'   per-feature activation rates (used to draw \code{psi} before
#'   \code{beta} itself is updated).
#' @param pi_g The current prior probability that a feature belongs
#'   to the slab component.
#' @param b1,b2 Shape hyperparameters of the slab component,
#'   \eqn{\mathrm{Beta}(b_1, b_2)}.
#' @param c1,c2 Shape hyperparameters of the spike component,
#'   \eqn{\mathrm{Beta}(c_1, c_2)}.
#'
#' @return A list with elements \code{beta} (the updated activation
#'   rates), \code{psi} (the updated sparsity indicators), and
#'   \code{zeta} (the posterior probability that each feature belongs
#'   to the slab, prior to sampling \code{psi}).
#'
#' @details
#' \code{psi_g} is drawn first, from the \emph{current} \code{beta_g},
#' with
#' \deqn{\zeta_g = \frac{\pi/B(b_1,b_2) \cdot \beta_g^{b_1-1}
#'  (1-\beta_g)^{b_2-1}}{\pi/B(b_1,b_2) \cdot \beta_g^{b_1-1}
#'  (1-\beta_g)^{b_2-1} + (1-\pi)/B(c_1,c_2) \cdot \beta_g^{c_1-1}
#'  (1-\beta_g)^{c_2-1}}.}
#' Unlike \code{beta_g}'s own full conditional (where a single
#' branch's Beta normalizing constant cancels out of the
#' proportionality), \code{psi_g}'s full conditional compares the slab
#' and spike branches directly against each other, so the normalizing
#' constants \eqn{B(b_1,b_2)} and \eqn{B(c_1,c_2)} must be retained
#' whenever \eqn{(b_1,b_2) \neq (c_1,c_2)} -- computed here via
#' \code{lbeta()} for numerical stability. \code{beta_g} is then
#' redrawn from whichever component \code{psi_g} selects, using
#' \eqn{n_g = \sum_r C_{rg}}:
#' \deqn{\beta_g \mid \mathrm{rest} \sim \mathrm{Beta}(n_g+b_1,\,
#'  R-n_g+b_2) \ \text{if}\ \psi_g=1, \quad \mathrm{Beta}(n_g+c_1,\,
#'  R-n_g+c_2) \ \text{if}\ \psi_g=0.}
#'
#' @export
update_beta_psi <- function(C, beta, pi_g, b1, b2, c1, c2) {

  R <- nrow(C)
  G <- ncol(C)
  n_g <- colSums(C)

  log_slab  <- log(pi_g)     - lbeta(b1, b2) + (b1 - 1) * log(beta) + (b2 - 1) * log(1 - beta)
  log_spike <- log(1 - pi_g) - lbeta(c1, c2) + (c1 - 1) * log(beta) + (c2 - 1) * log(1 - beta)
  zeta <- 1 / (1 + exp(log_spike - log_slab))

  psi <- stats::rbinom(G, 1, zeta)

  beta_new <- ifelse(
    psi == 1,
    stats::rbeta(G, n_g + b1, R - n_g + b2),
    stats::rbeta(G, n_g + c1, R - n_g + c2)
  )

  return(list(beta = beta_new, psi = psi, zeta = zeta))
}

#' Compute the Unnormalized Log-Posterior
#'
#' Computes the full unnormalized log-posterior density at one
#' complete set of parameter values, summing the data log-likelihood
#' and every prior log-density in the model.
#'
#' @param X An observed binary matrix of dimension \eqn{K \times G}.
#' @param W A binary matrix of dimension \eqn{K \times R}.
#' @param H A binary matrix of dimension \eqn{R \times G}.
#' @param alpha A numeric vector of length \eqn{K}.
#' @param beta A numeric vector of length \eqn{G}.
#' @param psi A binary vector of length \eqn{G}.
#' @param pi_g A numeric scalar.
#' @param phi A numeric vector of length \eqn{R}.
#' @param p11,p10 Numeric scalars, the noise parameters.
#' @param a1,a2 Beta hyperparameters for \code{alpha}.
#' @param b1,b2 Beta hyperparameters for the slab component of
#'   \code{beta}.
#' @param c1,c2 Beta hyperparameters for the spike component of
#'   \code{beta}.
#' @param d1,d2 Beta hyperparameters for \code{pi_g}.
#' @param e1,e2 Beta hyperparameters for \code{p11}.
#' @param f1,f2 Beta hyperparameters for \code{p10}.
#' @param s1,s2 Beta hyperparameters for \code{phi}.
#'
#' @return A numeric scalar, the unnormalized log-posterior at the
#'   supplied parameter values.
#'
#' @details
#' This is the quantity used to identify the Maximum A Posteriori
#' (MAP) sample across a Gibbs chain: the iteration attaining the
#' largest value of this function is reported as the point estimate
#' of the latent factor matrices, since the posterior mean would give
#' non-binary matrices. The same 0*(-Inf)=NaN floating-point guard
#' used in \code{\link{compute_log_prior_phi}} is applied to the
#' \code{alpha} and noise-parameter prior terms, which are otherwise
#' vulnerable to the same trap whenever a hyperparameter equals 1 and
#' the corresponding parameter lies exactly at a 0/1 boundary (as
#' happens by construction under several of the ablation restrictions
#' in \code{\link{run_BBMF_ablation}}).
#'
#' @export
compute_log_posterior <- function(X, W, H, alpha, beta, psi, pi_g, phi, p11, p10,
                                   a1, a2, b1, b2, c1, c2, d1, d2, e1, e2, f1, f2, s1, s2) {

  ll <- compute_likelihood(W, H, X, p11, p10)
  lp_W <- compute_log_prior_W(W, alpha, phi)
  lp_H <- compute_log_prior_H(H, beta, phi)

  lp_alpha <- sum(ifelse(a1 == 1, 0, (a1 - 1) * log(alpha)) +
                  ifelse(a2 == 1, 0, (a2 - 1) * log(1 - alpha)))

  lp_beta_psi <- sum(ifelse(psi == 1,
                             -lbeta(b1, b2) + (b1 - 1) * log(beta) + (b2 - 1) * log(1 - beta),
                             -lbeta(c1, c2) + (c1 - 1) * log(beta) + (c2 - 1) * log(1 - beta)))

  lp_psi <- sum(psi * log(pi_g) + (1 - psi) * log(1 - pi_g))
  lp_pi  <- (d1 - 1) * log(pi_g) + (d2 - 1) * log(1 - pi_g)
  lp_phi <- compute_log_prior_phi(phi, s1, s2)

  lp_noise <- ifelse(e1 == 1, 0, (e1 - 1) * log(p11)) + ifelse(e2 == 1, 0, (e2 - 1) * log(1 - p11)) +
              ifelse(f1 == 1, 0, (f1 - 1) * log(p10)) + ifelse(f2 == 1, 0, (f2 - 1) * log(1 - p10))

  logpost <- ll + lp_W + lp_H + lp_alpha + lp_beta_psi + lp_psi + lp_pi + lp_phi + lp_noise

  return(logpost)
}

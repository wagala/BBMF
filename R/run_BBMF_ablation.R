#' Run a BBMF Ablation Variant (Array-Based Output, with Burn-In)
#'
#' Runs the full Gibbs sampler for one variant of the nested ablation
#' family (Section 2.3: \code{"Full"}, \code{"NoPhi"}, \code{"Flat"},
#' \code{"Global"}, \code{"Sym"}) for a single chain, discarding an
#' initial burn-in period and returning posterior samples as
#' pre-allocated arrays/matrices/vectors, in the same format as
#' \code{\link{run_BBMF}}.
#'
#' @param X An observed binary matrix of dimension \eqn{K \times G}.
#' @param R The number of latent factors (upper bound).
#' @param n_iter Total number of Gibbs sweeps to run, including
#'   burn-in.
#' @param burn_in Number of initial sweeps to discard before storing
#'   samples. Default \code{0}.
#' @param variant One of \code{"Full"}, \code{"NoPhi"}, \code{"Flat"},
#'   \code{"Global"}, \code{"Sym"} (see \code{\link{ablation_flags}}).
#'   Default \code{"Full"}, in which case this function behaves
#'   identically to \code{\link{run_BBMF}}.
#' @param W_init,H_init,init_method,asso_seed See
#'   \code{\link{run_BBMF}}; identical initialization logic.
#' @param a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,f1,f2,s1,s2 Hyperparameters;
#'   see \code{\link{run_BBMF}} for defaults and meaning. Not all are
#'   used by every variant (e.g.\ \code{c1,c2,d1,d2} are irrelevant
#'   once the spike-and-slab is disabled), but all are accepted for a
#'   uniform interface across variants.
#' @param show_progress Display a \code{progress} package progress
#'   bar. See \code{\link{run_BBMF}} for details.
#'
#' @return A list with the same elements as \code{\link{run_BBMF}}
#'   (\code{samples_W}, \code{samples_H}, \code{samples_alpha},
#'   \code{samples_beta}, \code{samples_psi}, \code{samples_phi},
#'   \code{samples_pi}, \code{samples_p11}, \code{samples_p10},
#'   \code{logpost_trace}), plus \code{variant} recording which
#'   restriction was used.
#'
#' @details
#' Under \code{"NoPhi"}, \code{"Flat"}, \code{"Global"}, and
#' \code{"Sym"}, \eqn{\phi_r} is fixed at 1 throughout (never sampled),
#' and by Corollary 1 the auxiliary variables \eqn{A=W}, \eqn{C=H}
#' deterministically -- so no auxiliary-variable sampling is performed
#' for these variants, only for \code{"Full"}.
#'
#' Under \code{"Global"} and \code{"Sym"}, \code{alpha} and
#' \code{beta} are single pooled scalars each sweep
#' (\code{\link{update_alpha_global}}, \code{\link{update_beta_global}}),
#' but are broadcast to full length (\eqn{K} and \eqn{G} respectively)
#' when stored in \code{samples_alpha}/\code{samples_beta}, so that
#' every downstream function (metrics, plotting, trace extraction)
#' works identically regardless of which variant produced the fit,
#' with no special-casing required. The same applies to \code{phi}
#' (always stored as a length-\eqn{R} vector, even when fixed at 1 for
#' every entry).
#'
#' Under \code{"Sym"}, the noise parameters are updated via
#' \code{\link{update_rho}} and \eqn{p_{10}} is set to \eqn{1-p_{11}}
#' exactly; the stored \code{samples_p11}/\code{samples_p10} therefore
#' satisfy this relationship exactly at every retained iteration.
#'
#' @seealso \code{\link{run_BBMF}} for the unablated model,
#'   \code{\link{run_multiple_BBMF_ablation}} for running several
#'   chains in parallel, \code{\link{ablation_flags}} for the
#'   variant-to-restriction mapping.
#'
#' @export
run_BBMF_ablation <- function(X, R, n_iter, burn_in = 0,
                               variant = c("Full", "NoPhi", "Flat", "Global", "Sym"),
                               W_init = NULL, H_init = NULL,
                               init_method = c("asso", "random"), asso_seed = NULL,
                               a1 = 1, a2 = 1, b1 = 1, b2 = 1, c1 = 1, c2 = 20, d1 = 1, d2 = 1,
                               e1 = 1, e2 = 1, f1 = 1, f2 = 1, s1 = 1, s2 = 1,
                               show_progress = TRUE) {

  variant <- match.arg(variant)
  flags <- ablation_flags(variant)
  init_method <- match.arg(init_method)

  K <- nrow(X)
  G <- ncol(X)

  if (burn_in >= n_iter) {
    stop("burn_in must be strictly smaller than n_iter.")
  }
  if (xor(is.null(W_init), is.null(H_init))) {
    stop("W_init and H_init must be supplied together, or not at all.")
  }

  if (!is.null(W_init)) {
    stopifnot(nrow(W_init) == K, ncol(W_init) == R)
    stopifnot(nrow(H_init) == R, ncol(H_init) == G)
    W <- W_init
    H <- H_init
  } else if (init_method == "asso") {
    if (!requireNamespace("rBMF", quietly = TRUE)) {
      stop(
        "init_method = \"asso\" (the default) requires the 'rBMF' package, ",
        "which is not installed. Either install it (install.packages(\"rBMF\")), ",
        "supply W_init and H_init directly, or set init_method = \"random\".",
        call. = FALSE
      )
    }
    asso_fit <- Asso_BBMF(X, R, seed = asso_seed)
    W <- asso_fit$W_asso
    H <- asso_fit$H_asso
  } else {
    W <- matrix(stats::rbinom(K * R, 1, 0.2), K, R)
    H <- matrix(stats::rbinom(R * G, 1, 0.2), R, G)
  }

  phi <- if (flags$fix_phi) rep(1, R) else rep(0.8, R)
  if (flags$global_priors) {
    alpha <- 0.3
    beta  <- 0.3
  } else {
    alpha <- rep(0.3, K)
    beta  <- rep(0.3, G)
  }
  psi  <- rep(1L, G)
  pi_g <- 0.5
  if (flags$symmetric_noise) {
    rho <- 0.9
    p11 <- rho
    p10 <- 1 - rho
  } else {
    p11 <- 0.9
    p10 <- 0.1
  }

  n_keep <- n_iter - burn_in

  samples_W     <- array(0L, dim = c(K, R, n_keep))
  samples_H     <- array(0L, dim = c(R, G, n_keep))
  samples_alpha <- matrix(0, nrow = K, ncol = n_keep)
  samples_beta  <- matrix(0, nrow = G, ncol = n_keep)
  samples_psi   <- matrix(0L, nrow = G, ncol = n_keep)
  samples_phi   <- matrix(0, nrow = R, ncol = n_keep)
  samples_pi    <- numeric(n_keep)
  samples_p11   <- numeric(n_keep)
  samples_p10   <- numeric(n_keep)
  logpost_trace <- numeric(n_keep)
  loglik_trace     <- numeric(n_keep)
  logprior_W_trace <- numeric(n_keep)
  logprior_H_trace <- numeric(n_keep)

  progress_available <- requireNamespace("progress", quietly = TRUE)
  if (show_progress && !progress_available) {
    message(
      "Package 'progress' not installed; proceeding without a progress bar. ",
      "Install it with install.packages(\"progress\") for progress display."
    )
  }
  use_progress <- show_progress && interactive() && progress_available
  if (use_progress) {
    pb_chain <- progress::progress_bar$new(
      format = paste0("  [", variant, "] sampling [:bar] :percent eta: :eta"),
      total = n_iter, clear = FALSE, width = 60
    )
  }

  for (iter in seq_len(n_iter)) {

    W <- update_W(W, H, X, alpha, phi, p11, p10)
    H <- update_H(W, H, X, beta, phi, p11, p10)

    if (!flags$fix_phi) {
      AB <- update_AB(W, alpha, phi)
      CD <- update_CD(H, beta, phi)
      A <- AB$A; B <- AB$B
      C <- CD$C; D <- CD$D
      phi <- update_phi(B, D, s1, s2)
    } else {
      # Corollary 1: phi = 1 makes A = W, C = H deterministically --
      # no auxiliary-variable sampling needed for these variants.
      A <- W
      C <- H
    }

    if (flags$global_priors) {
      alpha <- update_alpha_global(A, a1, a2)
      beta  <- update_beta_global(C, b1, b2)
    } else if (flags$disable_spike_slab) {
      alpha <- update_alpha(A, a1, a2)
      beta  <- update_beta_flat(C, b1, b2)
    } else {
      alpha <- update_alpha(A, a1, a2)
      bp <- update_beta_psi(C, beta, pi_g, b1, b2, c1, c2)
      beta <- bp$beta
      psi <- bp$psi
      pi_g <- update_pi(psi, d1, d2)
    }

    Z <- bool_prod(W, H)
    if (flags$symmetric_noise) {
      rho <- update_rho(X, Z, e1, e2)
      p11 <- rho
      p10 <- 1 - rho
    } else {
      np <- update_noise(X, Z, p11, p10, e1, e2, f1, f2)
      p11 <- np$p11
      p10 <- np$p10
    }

    if (iter > burn_in) {
      idx <- iter - burn_in
      alpha_full <- if (length(alpha) == 1) rep(alpha, K) else alpha
      beta_full  <- if (length(beta)  == 1) rep(beta,  G) else beta

      samples_W[, , idx]   <- W
      samples_H[, , idx]   <- H
      samples_alpha[, idx] <- alpha_full
      samples_beta[, idx]  <- beta_full
      samples_psi[, idx]   <- psi
      samples_phi[, idx]   <- phi
      samples_pi[idx]      <- pi_g
      samples_p11[idx]     <- p11
      samples_p10[idx]     <- p10
      logpost_trace[idx]   <- compute_log_posterior(
        X, W, H, alpha_full, beta_full, psi, pi_g, phi, p11, p10,
        a1, a2, b1, b2, c1, c2, d1, d2, e1, e2, f1, f2, s1, s2
      )
      # Decomposed pieces: computed fresh from the CURRENT W, H (already in
      # memory this iteration) -- no extra storage needed beyond one scalar
      # each, discarded from memory immediately, unlike storing full W/H.
      loglik_trace[idx]     <- compute_likelihood(W, H, X, p11, p10)
      logprior_W_trace[idx] <- compute_log_prior_W(W, alpha_full, phi)
      logprior_H_trace[idx] <- compute_log_prior_H(H, beta_full, phi)
    }

    if (use_progress) pb_chain$tick()
  }

  if (use_progress) pb_chain$terminate()

  return(list(
    variant = variant,
    samples_W = samples_W, samples_H = samples_H,
    samples_alpha = samples_alpha, samples_beta = samples_beta,
    samples_psi = samples_psi, samples_phi = samples_phi,
    samples_pi = samples_pi, samples_p11 = samples_p11, samples_p10 = samples_p10,
    logpost_trace = logpost_trace,
    loglik_trace = loglik_trace,
    logprior_W_trace = logprior_W_trace,
    logprior_H_trace = logprior_H_trace
  ))
}

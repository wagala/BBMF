#' Run the BBMF Gibbs Sampler (Array-Based Output, with Burn-In)
#'
#' Runs the full Gibbs sampler for the Bayesian Boolean Matrix
#' Factorization model for a single chain, discarding an initial
#' burn-in period and returning posterior samples as pre-allocated
#' arrays/matrices/vectors rather than a list of per-iteration
#' snapshots.
#'
#' @param X An observed binary matrix of dimension \eqn{K \times G}.
#' @param R The number of latent factors (upper bound; the effective
#'   number in use can fall below \code{R} via the factor-level
#'   sparsity mechanism, see \code{\link{update_phi}}).
#' @param n_iter Total number of Gibbs sweeps to run, including
#'   burn-in.
#' @param burn_in Number of initial sweeps to discard before storing
#'   samples. Must be strictly less than \code{n_iter}. Default
#'   \code{0} (no burn-in discarded).
#' @param W_init,H_init,init_method,asso_seed Initialization controls for the factor matrices and optional Asso warm start.
#' @param a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,f1,f2,s1,s2 Beta-prior hyperparameters; defaults are given in the function signature.
#' @param show_progress Display a \code{progress} package progress
#'   bar tracking iteration count (including burn-in). Default
#'   \code{TRUE}. Requires the \code{progress} package; if it is not
#'   installed, sampling proceeds without a progress bar and a
#'   one-time message is printed. Automatically ignored (treated as
#'   \code{FALSE}) in non-interactive sessions.
#'
#' @return A list with elements \code{samples_W} (array,
#'   \eqn{K \times R \times (n\_iter - burn\_in)}), \code{samples_H}
#'   (array, \eqn{R \times G \times (n\_iter - burn\_in)}),
#'   \code{samples_alpha} (matrix, \eqn{K \times (n\_iter - burn\_in)}),
#'   \code{samples_beta} (matrix, \eqn{G \times (n\_iter - burn\_in)}),
#'   \code{samples_psi} (matrix, \eqn{G \times (n\_iter - burn\_in)}),
#'   \code{samples_phi} (matrix, \eqn{R \times (n\_iter - burn\_in)}),
#'   \code{samples_pi}, \code{samples_p11}, \code{samples_p10}, and
#'   \code{logpost_trace} (each a numeric vector of length
#'   \code{n_iter - burn_in}).
#'
#' @details
#' Posterior samples are stored in pre-allocated arrays using the same
#' underlying model and identical per-sweep update logic (the same
#' \code{\link{update_W}}, \code{\link{update_H}}, etc.\ functions are
#' used internally, in the same order), differing only in how results
#' are returned: as pre-allocated arrays indexed by post-burn-in
#' iteration, rather than a list of one snapshot per iteration. The
#' array layout keeps memory contiguous per parameter (e.g.\ all of
#' \code{samples_W} in one block) rather than per iteration, which can
#' be more convenient for downstream array-based post-processing (e.g.
#' \code{samples_W[,, idx]} for the \eqn{K \times R} draw at retained
#' iteration \code{idx}, or \code{apply(samples_phi, 1, mean)} for the
#' posterior mean of each \eqn{\phi_r}).
#'
#' The Maximum A Posteriori (MAP) sample is identified after the fact
#' via \code{which.max(logpost_trace)}, e.g.:
#' \preformatted{
#' idx   <- which.max(fit$logpost_trace)
#' W_hat <- fit$samples_W[,, idx]
#' H_hat <- fit$samples_H[,, idx]
#' }
#'
#' @seealso \code{\link{run_BBMF_ablation}} for fitting a specified nested-family variant,
#'   \code{\link{run_multiple_BBMF}} for running several chains in
#'   parallel.
#'
#' @export
run_BBMF <- function(X, R, n_iter, burn_in = 0, W_init = NULL, H_init = NULL,
                      init_method = c("asso", "random"), asso_seed = NULL,
                      a1 = 1, a2 = 1, b1 = 1, b2 = 1, c1 = 1, c2 = 20, d1 = 1, d2 = 1,
                      e1 = 1, e2 = 1, f1 = 1, f2 = 1, s1 = 1, s2 = 1,
                      show_progress = TRUE) {

  K <- nrow(X)
  G <- ncol(X)
  init_method <- match.arg(init_method)

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

  alpha <- rep(0.3, K)
  beta  <- rep(0.3, G)
  psi   <- rep(1L, G)
  pi_g  <- 0.5
  phi   <- rep(0.8, R)
  p11   <- 0.9
  p10   <- 0.1

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
      format = "  sampling [:bar] :percent eta: :eta",
      total = n_iter, clear = FALSE, width = 60
    )
  }

  for (iter in seq_len(n_iter)) {

    W <- update_W(W, H, X, alpha, phi, p11, p10)
    H <- update_H(W, H, X, beta, phi, p11, p10)

    AB <- update_AB(W, alpha, phi)
    CD <- update_CD(H, beta, phi)

    alpha <- update_alpha(AB$A, a1, a2)
    bp <- update_beta_psi(CD$C, beta, pi_g, b1, b2, c1, c2)
    beta <- bp$beta
    psi <- bp$psi
    phi <- update_phi(AB$B, CD$D, s1, s2)
    pi_g <- update_pi(psi, d1, d2)

    Z <- bool_prod(W, H)
    np <- update_noise(X, Z, p11, p10, e1, e2, f1, f2)
    p11 <- np$p11
    p10 <- np$p10

    if (iter > burn_in) {
      idx <- iter - burn_in
      samples_W[, , idx]  <- W
      samples_H[, , idx]  <- H
      samples_alpha[, idx] <- alpha
      samples_beta[, idx]  <- beta
      samples_psi[, idx]   <- psi
      samples_phi[, idx]   <- phi
      samples_pi[idx]      <- pi_g
      samples_p11[idx]     <- p11
      samples_p10[idx]     <- p10
      logpost_trace[idx]   <- compute_log_posterior(
        X, W, H, alpha, beta, psi, pi_g, phi, p11, p10,
        a1, a2, b1, b2, c1, c2, d1, d2, e1, e2, f1, f2, s1, s2
      )
    }

    if (use_progress) pb_chain$tick()
  }

  if (use_progress) pb_chain$terminate()

  return(list(
    samples_W = samples_W, samples_H = samples_H,
    samples_alpha = samples_alpha, samples_beta = samples_beta,
    samples_psi = samples_psi, samples_phi = samples_phi,
    samples_pi = samples_pi, samples_p11 = samples_p11, samples_p10 = samples_p10,
    logpost_trace = logpost_trace
  ))
}

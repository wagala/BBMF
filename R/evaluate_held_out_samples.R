#' Held-Out Sample Prediction Evaluation
#'
#' Evaluates predictive accuracy by holding out entire samples (rows),
#' fitting the model on the remaining samples only, then inferring the
#' held-out samples' factor activations from the fitted \eqn{\mathbf{H}}
#' and predicting their reconstruction.
#'
#' @param X An observed binary matrix of dimension \eqn{K \times G}.
#' @param R The number of latent factors.
#' @param test_frac Fraction of samples (rows) to hold out. Default
#'   \code{0.2}.
#' @param n_iter_fit Number of Gibbs sweeps used to fit the model on
#'   the training samples.
#' @param n_iter_infer Number of Gibbs sweeps used to infer each
#'   held-out sample's own \eqn{W_{k\cdot}}, holding the fitted
#'   \eqn{\mathbf{H}} (and other fitted parameters) fixed.
#' @param seed Optional seed for the train/test split.
#' @param ... Additional arguments passed to \code{\link{run_BBMF}}
#'   for the training fit (e.g.\ \code{burn_in}, \code{init_method}).
#'
#' @return A list with elements \code{metrics} (the output of
#'   \code{\link{compute_metrics}} comparing the predicted held-out
#'   reconstruction against the actual held-out data), \code{Z_test_pred}
#'   (the predicted reconstruction for the held-out samples),
#'   \code{train_idx}, \code{test_idx}, and \code{fit} (the full
#'   training-set model fit, for further inspection).
#'
#' @details
#' This evaluates a genuinely different question than
#' \code{\link{compute_metrics}} applied to \code{mc_fit$best_sample}
#' on the full dataset: it asks whether the \emph{factor patterns}
#' learned from one set of samples generalize to predicting entirely
#' unseen samples, rather than how well the model fits the same data
#' it was trained on. It does not require any missing-data support in
#' the core sampler (unlike a per-cell masking evaluation, which this
#' package's \code{\link{update_W}}/\code{\link{update_H}} do not
#' currently support): held-out samples are simply excluded from the
#' training fit entirely, then their own \eqn{W_{k\cdot}} is inferred
#' afterward via ordinary Gibbs sweeps of \code{\link{update_W}},
#' holding the fitted \eqn{\mathbf{H}}, \eqn{\phi}, \eqn{p_{11}},
#' \eqn{p_{10}} fixed and using the mean of the training samples'
#' fitted \eqn{\alpha_k} as a shared prior rate for the new samples
#' (since no sample-specific \eqn{\alpha_k} exists for a sample that
#' was never part of the fit).
#'
#' @examples
#' \dontrun{
#' held_out <- evaluate_held_out_samples(X_noisy, R = 4, test_frac = 0.2,
#'                                        n_iter_fit = 50000, burn_in = 10000,
#'                                        seed = 1)
#' print(held_out$metrics)
#' }
#'
#' @seealso \code{\link{compute_metrics}}, \code{\link{run_BBMF}},
#'   \code{\link{update_W}}
#'
#' @export
evaluate_held_out_samples <- function(X, R, test_frac = 0.2,
                                       n_iter_fit = 10000, n_iter_infer = 500,
                                       seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)

  K <- nrow(X); G <- ncol(X)
  n_test <- round(test_frac * K)
  test_idx  <- sample(seq_len(K), size = n_test)
  train_idx <- setdiff(seq_len(K), test_idx)

  fit <- run_BBMF(X[train_idx, , drop = FALSE], R = R, n_iter = n_iter_fit, ...)
  best_idx <- which.max(fit$logpost_trace)

  H_hat   <- fit$samples_H[, , best_idx]
  p11_hat <- fit$samples_p11[best_idx]
  p10_hat <- fit$samples_p10[best_idx]
  phi_hat <- fit$samples_phi[, best_idx]
  alpha_new <- rep(mean(fit$samples_alpha[, best_idx]), n_test)

  X_test <- X[test_idx, , drop = FALSE]
  W_test <- matrix(stats::rbinom(n_test * R, 1, 0.2), n_test, R)
  for (it in seq_len(n_iter_infer)) {
    W_test <- update_W(W_test, H_hat, X_test, alpha_new, phi_hat, p11_hat, p10_hat)
  }

  Z_test_pred <- bool_prod(W_test, H_hat)
  metrics <- compute_metrics(Z_test_pred, X_test)

  list(metrics = metrics, Z_test_pred = Z_test_pred,
       train_idx = train_idx, test_idx = test_idx, fit = fit)
}

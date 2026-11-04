#' Gelman-Rubin Convergence Diagnostic (R-hat)
#'
#' Computes the potential scale reduction factor (\eqn{\hat{R}}) for a
#' scalar quantity tracked across multiple MCMC chains, comparing
#' between-chain variance to within-chain variance.
#'
#' @param chain_list A list of numeric vectors, one per chain, all of
#'   the same length (e.g.\ each chain's \code{logpost_trace}, or each
#'   chain's draws of a single scalar parameter such as \code{p11}).
#'
#' @return A numeric scalar. Values close to \code{1} (conventionally
#'   \code{< 1.01}) indicate the chains agree with each other; larger
#'   values indicate the chains have not yet converged to the same
#'   distribution.
#'
#' @details
#' Following the standard Gelman-Rubin construction: with \eqn{M}
#' chains of length \eqn{N}, between-chain variance
#' \eqn{B = \frac{N}{M-1}\sum_m (\bar\theta_m-\bar\theta)^2} and
#' within-chain variance \eqn{W = \frac{1}{M}\sum_m s_m^2} combine into
#' \deqn{\hat{R} = \sqrt{\frac{\frac{N-1}{N}W + \frac{1}{N}B}{W}}.}
#'
#' \strong{Caveat specific to warm-started chains:} \eqn{\hat{R}}'s
#' interpretation classically assumes chains start from
#' \emph{deliberately dispersed} initial points, so agreement across
#' chains is informative about convergence from different starting
#' regions. Chains produced by \code{\link{run_multiple_BBMF}} /
#' \code{\link{run_multiple_BBMF_ablation}} all warm-start from the
#' \emph{same} Asso solution by default, differing only in their
#' subsequent random draws -- so a low \eqn{\hat{R}} here is weaker
#' evidence of true convergence than it would be for chains started
#' from genuinely different regions of the parameter space, and should
#' be reported with this caveat rather than as an unqualified pass.
#'
#' @seealso \code{\link{compute_ess}} for the companion effective sample size diagnostic.
#'
#' @export
compute_rhat <- function(chain_list) {
  M <- length(chain_list)
  N <- length(chain_list[[1]])
  stopifnot(all(sapply(chain_list, length) == N))

  chain_means <- sapply(chain_list, mean)
  chain_vars  <- sapply(chain_list, stats::var)

  grand_mean <- mean(chain_means)
  B <- (N / (M - 1)) * sum((chain_means - grand_mean)^2)
  W <- mean(chain_vars)

  var_hat_plus <- ((N - 1) / N) * W + (1 / N) * B
  sqrt(var_hat_plus / W)
}

#' Effective Sample Size (ESS)
#'
#' Estimates the effective sample size of a scalar quantity tracked
#' across multiple (autocorrelated) MCMC chains, i.e.\ the number of
#' effectively independent draws the combined chains actually carry.
#'
#' @param chain_list A list of numeric vectors, one per chain, all of
#'   the same length; see \code{\link{compute_rhat}}.
#' @param max_lag Maximum autocorrelation lag to consider. Default
#'   \code{NULL}, which uses \code{length(chain_list[[1]]) - 1}.
#'
#' @return A numeric scalar, the estimated effective sample size,
#'   between \code{1} (a chain providing essentially one independent
#'   draw) and \code{M * N} (the chains behaving as if fully
#'   independent draws, the theoretical maximum).
#'
#' @details
#' Each chain is first centered by subtracting its own mean (removing
#' between-chain mean differences, so only within-chain
#' autocorrelation structure remains), then all centered chains are
#' pooled into one combined autocorrelation function. Following
#' Geyer's initial positive sequence estimator, consecutive pairs of
#' lag autocorrelations are summed only while their sum remains
#' positive (stopping at the first non-positive pair, which controls
#' the well-known noisiness of high-lag sample autocorrelation
#' estimates), giving
#' \deqn{\mathrm{ESS} = \frac{MN}{1 + 2\sum_{k} \hat\rho_k}.}
#' This is a standard, widely-used estimator (in the spirit of the one
#' implemented by the \code{coda} package), not a novel method -- but
#' is implemented here directly to avoid adding \code{coda} as a
#' dependency for two functions.
#'
#' @seealso \code{\link{compute_rhat}}
#'
#' @export
compute_ess <- function(chain_list, max_lag = NULL) {
  M <- length(chain_list)
  N <- length(chain_list[[1]])
  stopifnot(all(sapply(chain_list, length) == N))
  if (is.null(max_lag)) max_lag <- N - 1

  centered <- lapply(chain_list, function(x) x - mean(x))
  pooled <- unlist(centered)

  acf_vals <- stats::acf(pooled, lag.max = max_lag, plot = FALSE)$acf[, 1, 1]
  rho <- acf_vals[-1]  # drop lag 0 (always 1)

  # Geyer's initial positive sequence: sum consecutive PAIRS while their
  # sum stays positive, stopping at the first non-positive pair.
  sum_rho <- 0
  k <- 1
  while (k < length(rho)) {
    pair_sum <- rho[k] + rho[k + 1]
    if (is.na(pair_sum) || pair_sum <= 0) break
    sum_rho <- sum_rho + pair_sum
    k <- k + 2
  }

  ess <- (M * N) / (1 + 2 * sum_rho)
  max(1, min(ess, M * N))  # clamp to the theoretically valid range
}

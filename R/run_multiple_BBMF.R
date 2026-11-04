#' Run Multiple BBMF Chains in Parallel
#'
#' Runs several independent \code{\link{run_BBMF}} chains in parallel,
#' each with its own seed, using \code{pbmcapply::pbmclapply} for
#' fork-based parallelism with an aggregate progress display, and
#' identifies the single best sample (highest unnormalized
#' log-posterior) across all chains and iterations combined.
#'
#' @param X An observed binary matrix of dimension \eqn{K \times G}.
#' @param R The number of latent factors.
#' @param n_iter Total number of Gibbs sweeps per chain, including
#'   burn-in.
#' @param burn_in Number of initial sweeps to discard per chain before
#'   storing samples. Default \code{0}.
#' @param n_chains Number of independent chains to run.
#' @param base_seed Integer added to each chain index to form that
#'   chain's seed (\code{set.seed(base_seed + chain_id)}). Default
#'   \code{1001}, matching \code{chain_id} running from \code{1} to
#'   \code{n_chains}.
#' @param mc.cores Number of cores to use. If \code{NULL} (the
#'   default), uses \code{min(n_chains, parallel::detectCores())}.
#' @param ... Additional arguments passed through to
#'   \code{\link{run_BBMF}} (e.g. hyperparameters, \code{init_method}).
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{chains}}{A named list of length \code{n_chains}
#'     (\code{"Chain_1"}, ..., \code{"Chain_n"}), each element the
#'     full array-based output of \code{\link{run_BBMF}} for that
#'     chain.}
#'   \item{\code{best_chain}}{Integer index (1 to \code{n_chains}) of
#'     the chain containing the single best sample.}
#'   \item{\code{best_chain_name}}{The corresponding \code{"Chain_k"}
#'     name.}
#'   \item{\code{best_iter}}{Index (within the post-burn-in samples of
#'     \code{best_chain}) of the best sample.}
#'   \item{\code{best_logpost}}{The log-posterior at that sample.}
#'   \item{\code{best_sample}}{A list with elements \code{W}, \code{H},
#'     \code{alpha}, \code{beta}, \code{psi}, \code{phi}, \code{pi},
#'     \code{p11}, \code{p10}, giving the parameter values at the best
#'     sample -- i.e.\ the relevant slice already extracted from each
#'     array/matrix/vector in \code{chains}, so no manual indexing
#'     into \code{samples_W[,, idx]} etc.\ is required.}
#' }
#'
#' @details
#' Each chain's own internal progress bar (see \code{\link{run_BBMF}}'s
#' \code{show_progress} argument) is forced off here, regardless of
#' any value passed via \code{...}: with \code{n_chains} chains running
#' in separate forked processes, several independent per-iteration
#' progress bars writing to the same console would interleave into
#' unreadable output. Progress is instead shown once, in aggregate,
#' across chains, via \code{pbmcapply}'s own \code{mc.style = "ETA"}
#' display.
#'
#' Fork-based parallelism (the mechanism underlying
#' \code{pbmcapply::pbmclapply}, and \code{parallel::mclapply} beneath
#' it) is not available on Windows; there, chains run sequentially on
#' a single core regardless of \code{mc.cores}, with no error but also
#' no speedup.
#'
#' The function also provides a convenience summary across all fitted chains:
#' rather than only returning each chain's raw arrays, the single best
#' sample across every chain and every retained iteration is
#' identified and extracted automatically.
#'
#' @seealso \code{\link{run_BBMF}} for fitting a single chain.
#'
#' @export
run_multiple_BBMF <- function(X, R, n_iter, burn_in = 0, n_chains,
                               base_seed = 1001, mc.cores = NULL, ...) {

  if (!requireNamespace("pbmcapply", quietly = TRUE)) {
    stop(
      "Package 'pbmcapply' is required for run_multiple_BBMF(). ",
      "Install it with install.packages(\"pbmcapply\"), or use ",
      "Alternatively, run individual chains with run_BBMF().",
      call. = FALSE
    )
  }

  if (is.null(mc.cores)) {
    mc.cores <- min(n_chains, parallel::detectCores())
  }

  run_chain <- function(chain_id) {
    set.seed(base_seed + chain_id)
    run_BBMF(X = X, R = R, n_iter = n_iter, burn_in = burn_in,
             show_progress = FALSE, ...)
  }

  all_chain_results <- pbmcapply::pbmclapply(
    X        = seq_len(n_chains),
    FUN      = run_chain,
    mc.cores = mc.cores,
    mc.style = "ETA"
  )

  names(all_chain_results) <- paste0("Chain_", seq_len(n_chains))

  # --- Identify the single best sample across all chains and iterations,
  # summarizing the best sample across all chains and retained iterations. ---
  best_per_chain <- sapply(all_chain_results, function(ch) max(ch$logpost_trace))
  best_chain <- which.max(best_per_chain)
  best_chain_name <- names(all_chain_results)[best_chain]
  best_iter <- which.max(all_chain_results[[best_chain]]$logpost_trace)
  best_logpost <- all_chain_results[[best_chain]]$logpost_trace[best_iter]

  bc <- all_chain_results[[best_chain]]
  best_sample <- list(
    W     = bc$samples_W[, , best_iter],
    H     = bc$samples_H[, , best_iter],
    alpha = bc$samples_alpha[, best_iter],
    beta  = bc$samples_beta[, best_iter],
    psi   = bc$samples_psi[, best_iter],
    phi   = bc$samples_phi[, best_iter],
    pi    = bc$samples_pi[best_iter],
    p11   = bc$samples_p11[best_iter],
    p10   = bc$samples_p10[best_iter]
  )

  return(list(
    chains          = all_chain_results,
    best_chain      = best_chain,
    best_chain_name = best_chain_name,
    best_iter       = best_iter,
    best_logpost    = best_logpost,
    best_sample     = best_sample
  ))
}

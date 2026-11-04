#' Run Multiple Chains of a BBMF Ablation Variant in Parallel
#'
#' Runs several independent \code{\link{run_BBMF_ablation}} chains in
#' parallel for a single specified ablation variant, each with its own
#' seed, using \code{pbmcapply::pbmclapply}, and identifies the single
#' best sample across all chains and iterations combined.
#'
#' @param X An observed binary matrix of dimension \eqn{K \times G}.
#' @param R The number of latent factors.
#' @param n_iter Total number of Gibbs sweeps per chain, including
#'   burn-in.
#' @param burn_in Number of initial sweeps to discard per chain before
#'   storing samples. Default \code{0}.
#' @param n_chains Number of independent chains to run.
#' @param variant One of \code{"Full"}, \code{"NoPhi"}, \code{"Flat"},
#'   \code{"Global"}, \code{"Sym"} (see \code{\link{ablation_flags}}).
#' @param base_seed Integer added to each chain index to form that
#'   chain's seed. Default \code{1001}.
#' @param mc.cores Number of cores to use. If \code{NULL} (the
#'   default), uses \code{min(n_chains, parallel::detectCores())}.
#' @param ... Additional arguments passed through to
#'   \code{\link{run_BBMF_ablation}} (e.g. hyperparameters,
#'   \code{init_method}, \code{W_init}, \code{H_init}).
#'
#' @return A list with elements \code{variant}, \code{chains} (a named
#'   list, \code{"Chain_1"}, ..., \code{"Chain_n"}, each the full
#'   array-based output of \code{\link{run_BBMF_ablation}}),
#'   \code{best_chain}, \code{best_chain_name}, \code{best_iter},
#'   \code{best_logpost}, and \code{best_sample} (a list with elements
#'   \code{W}, \code{H}, \code{alpha}, \code{beta}, \code{psi},
#'   \code{phi}, \code{pi}, \code{p11}, \code{p10}, already extracted
#'   from the relevant arrays).
#'
#' @details
#' As with \code{\link{run_multiple_BBMF}}, each chain's own internal
#' progress bar is forced off, relying solely on \code{pbmcapply}'s
#' aggregate \code{mc.style = "ETA"} display; fork-based parallelism is
#' unavailable on Windows, where chains run sequentially with no error
#' but also no speedup.
#'
#' To reproduce the full ablation comparison (Table 1, Section 4.4),
#' call this function once per variant with matching \code{n_iter},
#' \code{burn_in}, \code{n_chains}, and (typically) the same
#' \code{W_init}/\code{H_init} warm start:
#' \preformatted{
#' fits <- lapply(c("Sym","Global","Flat","NoPhi","Full"), function(v) {
#'   run_multiple_BBMF_ablation(X, R = 4, n_iter = 100000, burn_in = 20000,
#'                               n_chains = 4, variant = v,
#'                               W_init = W_asso, H_init = H_asso)
#' })
#' names(fits) <- c("Sym","Global","Flat","NoPhi","Full")
#' }
#' \code{run_multiple_BBMF(...)} is equivalent to
#' \code{run_multiple_BBMF_ablation(..., variant = "Full")}.
#'
#' @seealso \code{\link{run_BBMF_ablation}} for the single-chain
#'   sampler, \code{\link{run_multiple_BBMF}} for the unablated model,
#'   \code{\link{ablation_flags}} for the variant definitions.
#'
#' @export
run_multiple_BBMF_ablation <- function(X, R, n_iter, burn_in = 0, n_chains,
                                        variant = c("Full", "NoPhi", "Flat", "Global", "Sym"),
                                        base_seed = 1001, mc.cores = NULL, ...) {

  variant <- match.arg(variant)

  if (!requireNamespace("pbmcapply", quietly = TRUE)) {
    stop(
      "Package 'pbmcapply' is required for run_multiple_BBMF_ablation(). ",
      "Install it with install.packages(\"pbmcapply\").",
      call. = FALSE
    )
  }

  if (is.null(mc.cores)) {
    mc.cores <- min(n_chains, parallel::detectCores())
  }

  run_chain <- function(chain_id) {
    set.seed(base_seed + chain_id)
    run_BBMF_ablation(X = X, R = R, n_iter = n_iter, burn_in = burn_in,
                       variant = variant, show_progress = FALSE, ...)
  }

  all_chain_results <- pbmcapply::pbmclapply(
    X        = seq_len(n_chains),
    FUN      = run_chain,
    mc.cores = mc.cores,
    mc.style = "ETA"
  )

  names(all_chain_results) <- paste0("Chain_", seq_len(n_chains))

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
    variant         = variant,
    chains          = all_chain_results,
    best_chain      = best_chain,
    best_chain_name = best_chain_name,
    best_iter       = best_iter,
    best_logpost    = best_logpost,
    best_sample     = best_sample
  ))
}

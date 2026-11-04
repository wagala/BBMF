#' Extract Log-Likelihood, Log-Posterior, and Log-Prior Traces per Chain
#'
#' Builds a long-format data frame of four diagnostic quantities --
#' log-likelihood, log-posterior, and the log-priors of \code{W} and
#' \code{H} separately -- at every retained iteration of every chain in
#' a fitted multi-chain object, suitable for faceted trace plots.
#'
#' @param mc_fit The output of \code{\link{run_multiple_BBMF}} or
#'   \code{\link{run_multiple_BBMF_ablation}}: a list with a
#'   \code{chains} element, each chain holding array-based samples
#'   (\code{samples_W}, \code{samples_H}, \code{samples_alpha},
#'   \code{samples_beta}, \code{samples_phi}, \code{samples_p11},
#'   \code{samples_p10}, \code{logpost_trace}).
#' @param X The observed binary data matrix used to fit \code{mc_fit},
#'   needed to recompute the log-likelihood at each retained iteration.
#'
#' @return A data frame with one row per (chain, iteration) pair and
#'   columns \code{chain} (factor, labelled \code{"Chain 1"}, ...),
#'   \code{iteration}, \code{LogLikelihood}, \code{LogPosterior},
#'   \code{LogPrior_W}, and \code{LogPrior_H}.
#'
#' @details
#' \code{LogPosterior} is read directly from \code{logpost_trace}
#' (already computed during sampling); \code{LogLikelihood},
#' \code{LogPrior_W}, and \code{LogPrior_H} are recomputed after the
#' fact via \code{\link{compute_likelihood}}, \code{\link{compute_log_prior_W}},
#' and \code{\link{compute_log_prior_H}} respectively, from the stored
#' arrays. These four quantities do \emph{not} sum to anything exact:
#' \code{LogPosterior} additionally includes the priors on
#' \code{alpha}, \code{beta}/\code{psi}, \code{pi}, \code{phi}, and the
#' noise parameters, none of which are broken out here. They are four
#' informative slices of the posterior, not a complete additive
#' decomposition.
#'
#' This function expects the \strong{array-based} output of
#' \code{\link{run_multiple_BBMF}}/\code{\link{run_multiple_BBMF_ablation}}
#' (indexing \code{samples_W[,, idx]} etc.), not the list-based output
#' of the older list-based sampler (which stores one snapshot per
#' iteration instead).
#'
#' @examples
#' \dontrun{
#' trace_df <- extract_traces_detailed(mc_fit, X_noisy)
#'
#' trace_long <- reshape2::melt(
#'   trace_df, id.vars = c("chain", "iteration"),
#'   measure.vars = c("LogLikelihood", "LogPosterior", "LogPrior_W", "LogPrior_H"),
#'   variable.name = "Quantity", value.name = "Value"
#' )
#'
#' ggplot2::ggplot(trace_long, ggplot2::aes(x = iteration, y = Value)) +
#'   ggplot2::geom_line(alpha = 0.8, color = "steelblue") +
#'   ggplot2::facet_grid(chain ~ Quantity, scales = "free_y") +
#'   ggplot2::theme_minimal()
#' }
#'
#' @seealso \code{\link{extract_phi_traces}} for the companion
#'   phi-inclusion-probability trace extractor.
#'
#' @export
extract_traces_detailed <- function(mc_fit, X) {

  n_chains <- length(mc_fit$chains)

  traces_list <- lapply(seq_len(n_chains), function(c) {
    ch <- mc_fit$chains[[c]]
    n_kept_c <- length(ch$logpost_trace)

    loglik <- sapply(seq_len(n_kept_c), function(idx) {
      compute_likelihood(ch$samples_W[, , idx], ch$samples_H[, , idx], X,
                          ch$samples_p11[idx], ch$samples_p10[idx])
    })
    logprior_W <- sapply(seq_len(n_kept_c), function(idx) {
      compute_log_prior_W(ch$samples_W[, , idx], ch$samples_alpha[, idx], ch$samples_phi[, idx])
    })
    logprior_H <- sapply(seq_len(n_kept_c), function(idx) {
      compute_log_prior_H(ch$samples_H[, , idx], ch$samples_beta[, idx], ch$samples_phi[, idx])
    })

    data.frame(
      chain = factor(c, levels = seq_len(n_chains), labels = paste("Chain", seq_len(n_chains))),
      iteration = seq_len(n_kept_c),
      LogLikelihood = loglik,
      LogPosterior  = ch$logpost_trace,
      LogPrior_W    = logprior_W,
      LogPrior_H    = logprior_H
    )
  })

  do.call(rbind, traces_list)
}

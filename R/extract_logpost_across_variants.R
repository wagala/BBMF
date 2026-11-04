#' Extract Log-Posterior Traces Across an Entire Ablation Family
#'
#' Builds a long-format data frame of the log-posterior trajectory at
#' every retained iteration of every chain, across every variant in a
#' fitted ablation family, for convergence checking across the whole
#' family at once.
#'
#' @param fits A named list of multi-chain fits, one per ablation
#'   variant, e.g.\ the result of calling
#'   \code{\link{run_multiple_BBMF_ablation}} once per variant and
#'   collecting the results with \code{names(fits) <- variants}. Each
#'   element must have the same structure as the output of
#'   \code{\link{run_multiple_BBMF}}/\code{\link{run_multiple_BBMF_ablation}}
#'   (a \code{chains} list, each with a \code{logpost_trace} vector).
#'
#' @return A data frame with one row per (variant, chain, iteration)
#'   triple, and columns \code{variant} (character, using the names of
#'   \code{fits}), \code{chain} (factor, labelled \code{"Chain 1"},
#'   ...), \code{iteration}, and \code{logpost}.
#'
#' @details
#' This is the convergence check that should precede any comparison of
#' ablation variants by performance metric: if one variant's chains are
#' still visibly drifting while another's have stabilized, a
#' performance difference between them may reflect unequal convergence
#' rather than a genuine difference between the models. Plot with
#' \code{variant} faceted and \code{chain} colored (or vice versa) to
#' check both convergence within each variant and agreement across its
#' chains.
#'
#' @examples
#' \dontrun{
#' variants <- c("Sym", "Global", "Flat", "NoPhi", "Full")
#' fits <- lapply(variants, function(v) {
#'   run_multiple_BBMF_ablation(X_noisy, R = 4, n_iter = 100000, burn_in = 20000,
#'                               n_chains = 4, variant = v,
#'                               W_init = W_asso, H_init = H_asso)
#' })
#' names(fits) <- variants
#'
#' logpost_df <- extract_logpost_across_variants(fits)
#' logpost_df$variant <- factor(logpost_df$variant, levels = variants)
#'
#' ggplot2::ggplot(logpost_df, ggplot2::aes(x = iteration, y = logpost, color = chain)) +
#'   ggplot2::geom_line(alpha = 0.7) +
#'   ggplot2::facet_wrap(~ variant, scales = "free_y", ncol = 3) +
#'   ggplot2::theme_minimal()
#' }
#'
#' @seealso \code{\link{extract_afas_across_variants}} for the
#'   companion factor-recovery comparison across the same family.
#'
#' @export
extract_logpost_across_variants <- function(fits) {
  do.call(rbind, lapply(names(fits), function(v) {
    fit <- fits[[v]]
    n_chains <- length(fit$chains)
    do.call(rbind, lapply(seq_len(n_chains), function(c) {
      lp <- fit$chains[[c]]$logpost_trace
      data.frame(variant = v,
                 chain = factor(c, levels = seq_len(n_chains), labels = paste("Chain", seq_len(n_chains))),
                 iteration = seq_along(lp), logpost = lp)
    }))
  }))
}

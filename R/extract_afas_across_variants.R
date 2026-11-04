#' Extract AFAS Traces Across an Entire Ablation Family
#'
#' Builds a long-format data frame of the aggregate factor-alignment
#' similarity (AFAS) at every retained iteration of every chain, across
#' every variant in a fitted ablation family -- the payoff comparison
#' showing how reliably each restricted model recovers the true factor
#' structure, relative to the full model.
#'
#' @param fits A named list of multi-chain fits, one per ablation
#'   variant; see \code{\link{extract_logpost_across_variants}} for the
#'   expected structure.
#' @param H_true The true factor-pattern matrix used to generate the
#'   simulated data, of dimension \eqn{R \times G}.
#'
#' @return A data frame with one row per (variant, chain, iteration)
#'   triple, and columns \code{variant}, \code{chain} (factor,
#'   labelled \code{"Chain 1"}, ...), \code{iteration}, and \code{afas}.
#'
#' @details
#' Calls \code{\link{compute_afas}} once per retained iteration of
#' every chain of every variant, so runtime scales with the total
#' number of retained samples across the whole family; for a
#' full-scale run (e.g.\ 4 chains times 5 variants times tens of
#' thousands of retained iterations) this can take a noticeable amount
#' of time, though each individual call is cheap. Requires the
#' \code{clue} package (via \code{\link{compute_afas}}).
#'
#' As with \code{\link{extract_logpost_across_variants}}, this should
#' generally be examined only after confirming convergence: an AFAS
#' comparison across variants whose chains have not yet stabilized can
#' reflect unequal convergence rather than a genuine difference in
#' factor-recovery ability.
#'
#' @examples
#' \dontrun{
#' afas_df <- extract_afas_across_variants(fits, H_true)
#' afas_df$variant <- factor(afas_df$variant, levels = variants)
#'
#' # Trace view, matching extract_logpost_across_variants()'s layout
#' ggplot2::ggplot(afas_df, ggplot2::aes(x = iteration, y = afas, color = chain)) +
#'   ggplot2::geom_line(alpha = 0.7) +
#'   ggplot2::facet_wrap(~ variant, ncol = 3) +
#'   ggplot2::ylim(0, 1) +
#'   ggplot2::theme_minimal()
#'
#' # Violin + boxplot summary per variant, aggregating over all chains
#' # and iterations at once -- the single clearest "does the ablation
#' # hurt recovery" summary.
#' ggplot2::ggplot(afas_df, ggplot2::aes(x = variant, y = afas, fill = variant)) +
#'   ggplot2::geom_violin(alpha = 0.5, trim = FALSE) +
#'   ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA) +
#'   ggplot2::theme_minimal() + ggplot2::theme(legend.position = "none")
#' }
#'
#' @seealso \code{\link{compute_afas}} for the single-fit version this
#'   loops over, \code{\link{extract_logpost_across_variants}} for the
#'   companion convergence check across the same family.
#'
#' @export
extract_afas_across_variants <- function(fits, H_true) {
  do.call(rbind, lapply(names(fits), function(v) {
    fit <- fits[[v]]
    n_chains <- length(fit$chains)
    do.call(rbind, lapply(seq_len(n_chains), function(c) {
      ch <- fit$chains[[c]]
      n_kept <- dim(ch$samples_H)[3]
      afas_vals <- sapply(seq_len(n_kept), function(idx) compute_afas(ch$samples_H[, , idx], H_true)$afas)
      data.frame(variant = v,
                 chain = factor(c, levels = seq_len(n_chains), labels = paste("Chain", seq_len(n_chains))),
                 iteration = seq_along(afas_vals), afas = afas_vals)
    }))
  }))
}

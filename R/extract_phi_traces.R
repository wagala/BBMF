#' Extract Factor-Level Inclusion Probability (phi) Traces per Chain
#'
#' Builds a long-format data frame of every factor's \eqn{\phi_r}
#' trajectory at every retained iteration of every chain in a fitted
#' multi-chain object, suitable for faceted trace plots (one panel per
#' chain, all \eqn{R} factors overlaid within each panel).
#'
#' @param mc_fit The output of \code{\link{run_multiple_BBMF}} or
#'   \code{\link{run_multiple_BBMF_ablation}}: a list with a
#'   \code{chains} element, each chain holding \code{samples_phi} as
#'   an \eqn{R \times n\_kept} matrix.
#' @param R The number of latent factors (must match the \code{R} used
#'   to fit \code{mc_fit}).
#'
#' @return A data frame with one row per (chain, iteration) pair and
#'   columns \code{factor1}, ..., \code{factorR} (the \eqn{\phi_r}
#'   value for each factor), \code{chain} (factor, labelled
#'   \code{"Chain 1"}, ...), and \code{iteration}.
#'
#' @details
#' Since \code{samples_phi} is already stored as an
#' \eqn{R \times n\_kept} matrix per chain, extraction is a direct
#' transpose (\eqn{n\_kept \times R}) rather than a reconstruction from
#' per-iteration snapshots. For ablation variants with \code{phi} fixed
#' at 1 (\code{"Sym"}, \code{"Global"}, \code{"Flat"}, \code{"NoPhi"};
#' see \code{\link{ablation_flags}}), every column will be identically
#' \code{1} -- confirming this is a genuine sanity check worth running
#' on the real fitted object, not merely on the unit tests that already
#' establish it in isolation.
#'
#' This function expects the \strong{array-based} output of
#' \code{\link{run_multiple_BBMF}}/\code{\link{run_multiple_BBMF_ablation}},
#' not the list-based output of the older list-based sampler.
#'
#' @examples
#' \dontrun{
#' phi_df <- extract_phi_traces(mc_fit, R = R)
#'
#' phi_long <- reshape2::melt(
#'   phi_df, id.vars = c("chain", "iteration"),
#'   measure.vars = paste0("factor", seq_len(R)),
#'   variable.name = "Factor", value.name = "phi"
#' )
#'
#' ggplot2::ggplot(phi_long, ggplot2::aes(x = iteration, y = phi, color = Factor)) +
#'   ggplot2::geom_line(alpha = 0.8) +
#'   ggplot2::facet_wrap(~ chain, ncol = 2) +
#'   ggplot2::ylim(0, 1) +
#'   ggplot2::theme_minimal()
#' }
#'
#' @seealso \code{\link{extract_traces_detailed}} for the companion
#'   likelihood/posterior/prior trace extractor.
#'
#' @export
extract_phi_traces <- function(mc_fit, R) {

  n_chains <- length(mc_fit$chains)

  traces_list <- lapply(seq_len(n_chains), function(c) {
    ch <- mc_fit$chains[[c]]
    phi_mat <- t(ch$samples_phi)   # transpose: n_kept x R (samples_phi is R x n_kept)
    df <- as.data.frame(phi_mat)
    colnames(df) <- paste0("factor", seq_len(R))
    df$chain <- factor(c, levels = seq_len(n_chains), labels = paste("Chain", seq_len(n_chains)))
    df$iteration <- seq_len(nrow(df))
    df
  })

  do.call(rbind, traces_list)
}

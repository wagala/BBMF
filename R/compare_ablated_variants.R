#' Compare Ablated BBMF Variants Pairwise
#'
#' Performs all pairwise paired Wilcoxon signed-rank tests among a specified
#' set of ablated BBMF variants across matched simulation replicates.
#'
#' The function is intended for comparisons among restricted members of the
#' nested BBMF family after excluding BBMF-Full. By default, it compares
#' `"Sym"`, `"Global"`, `"Flat"`, and `"NoPhi"`, giving six pairwise tests.
#' P-values are adjusted jointly for multiple comparisons.
#'
#' Replicates are paired through the `replicate_id` column, so each model
#' variant should have a value of the selected performance metric for the
#' same simulation replicates.
#'
#' @param combined_metrics A data frame containing simulation performance
#'   results. It must include columns `Method`, `replicate_id`, and the
#'   performance metric specified by `metric`.
#' @param metric Character string giving the performance metric to compare.
#'   Default is `"F1"`.
#' @param variants Character vector giving the BBMF variants to compare.
#'   Default is `c("Sym", "Global", "Flat", "NoPhi")`.
#' @param adjust_method Character string specifying the multiple-comparison
#'   adjustment method passed to [stats::p.adjust()]. Default is `"holm"`.
#'
#' @return A data frame with one row for each pairwise comparison and columns:
#' \describe{
#'   \item{Variant1}{First variant in the comparison.}
#'   \item{Variant2}{Second variant in the comparison.}
#'   \item{MeanDiff}{Mean paired difference, `Variant1 - Variant2`.}
#'   \item{V}{Wilcoxon signed-rank test statistic.}
#'   \item{p_raw}{Unadjusted p-value.}
#'   \item{p_adjusted}{P-value after multiple-comparison correction.}
#' }
#'
#' @details
#' Positive values of `MeanDiff` indicate that `Variant1` has a larger
#' average value of the selected metric than `Variant2`. Interpretation of
#' the direction therefore depends on whether larger or smaller values of
#' the chosen metric indicate better performance.
#'
#' @examples
#' \dontrun{
#' compare_ablated_variants(
#'   combined_metrics,
#'   metric = "F1"
#' )
#' }
#'
#' @export
compare_ablated_variants <- function(combined_metrics, metric = "F1",
                                     variants = c("Sym", "Global", "Flat", "NoPhi"),
                                     adjust_method = "holm") {
  
  # Reshape to wide: one row per replicate, one column per variant
  wide <- stats::reshape(
    combined_metrics[
      combined_metrics$Method %in% variants,
      c("Method", metric, "replicate_id")
    ],
    idvar = "replicate_id",
    timevar = "Method",
    direction = "wide"
  )
  
  colnames(wide) <- sub(
    paste0("^", metric, "\\."),
    "",
    colnames(wide)
  )
  
  pairs <- utils::combn(variants, 2, simplify = FALSE)
  
  results <- do.call(rbind, lapply(pairs, function(p) {
    
    test <- stats::wilcox.test(
      wide[[p[1]]],
      wide[[p[2]]],
      paired = TRUE
    )
    
    data.frame(
      Variant1 = p[1],
      Variant2 = p[2],
      MeanDiff = mean(wide[[p[1]]] - wide[[p[2]]]),
      V = unname(test$statistic),
      p_raw = test$p.value
    )
  }))
  
  results$p_adjusted <- stats::p.adjust(
    results$p_raw,
    method = adjust_method
  )
  
  results <- results[order(results$p_raw), ]
  rownames(results) <- NULL
  
  results
}
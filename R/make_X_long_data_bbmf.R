#' Convert BBMF Posterior Reconstruction Summaries to Long Format
#'
#' Converts posterior reconstruction probability and uncertainty matrices
#' into long-format data frames suitable for plotting.
#'
#' @param probability A numeric matrix containing posterior reconstruction
#'   probabilities.
#' @param uncertainty A numeric matrix containing posterior reconstruction
#'   uncertainty scores.
#' @param X_truth A reference matrix used to determine sample and feature
#'   ordering.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{probability}}{
#'     A long-format data frame containing sample, feature, and posterior
#'     reconstruction probability.
#'   }
#'   \item{\code{uncertainty}}{
#'     A long-format data frame containing sample, feature, and posterior
#'     uncertainty.
#'   }
#' }
#'
#' @examples
#' \dontrun{
#' X_long <- make_X_long_data_bbmf(
#'   probability = X_summary$probability,
#'   uncertainty = X_summary$uncertainty,
#'   X_truth = X_true
#' )
#' }
#'
#' @export
#'
make_X_long_data_bbmf <- function(
    probability,
    uncertainty,
    X_truth) {
  
  probability_df <- as.data.frame(
    as.table(probability),
    stringsAsFactors = FALSE
  )
  
  names(probability_df) <- c(
    "Sample",
    "Feature",
    "Probability"
  )
  
  
  uncertainty_df <- as.data.frame(
    as.table(uncertainty),
    stringsAsFactors = FALSE
  )
  
  names(uncertainty_df) <- c(
    "Sample",
    "Feature",
    "Uncertainty"
  )
  
  
  sample_order <- rownames(X_truth)
  feature_order <- colnames(X_truth)
  
  
  probability_df$Sample <- factor(
    probability_df$Sample,
    levels = rev(sample_order)
  )
  
  probability_df$Feature <- factor(
    probability_df$Feature,
    levels = feature_order
  )
  
  
  uncertainty_df$Sample <- factor(
    uncertainty_df$Sample,
    levels = rev(sample_order)
  )
  
  uncertainty_df$Feature <- factor(
    uncertainty_df$Feature,
    levels = feature_order
  )
  
  
  return(
    list(
      probability = probability_df,
      uncertainty = uncertainty_df
    )
  )
}
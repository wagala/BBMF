#' Compute Aligned Mean Similarity Across Posterior Samples
#'
#' Computes the mean factor-wise Jaccard similarity between a reference
#' factor matrix and each posterior sample of a factor matrix after
#' optimal factor alignment using the Hungarian algorithm.
#'
#' @param H_true A binary reference matrix of dimension
#'   \eqn{R \times G}, typically the true factor-feature matrix in a
#'   simulation study.
#' @param ff A three-dimensional binary array of dimension
#'   \eqn{R \times G \times M}, where \eqn{M} is the number of posterior
#'   samples.
#'
#' @return A numeric vector of length \eqn{M}. Each element contains
#'   the mean Jaccard similarity between \code{H_true} and the
#'   corresponding posterior sample after optimal factor alignment.
#'
#' @export
#'
get_mean_sim_per_chain <- function(H_true, ff) {
  
  M <- dim(ff)[3]
  
  mean_sims <- numeric(M)
  
  for (m in seq_len(M)) {
    
    H_m <- ff[, , m]
    
    S <- jaccard_rows(
      H_true,
      H_m
    )
    
    align <- factor_alignment_hungarian(
      S
    )
    
    agg <- aggregate_alignment_similarity(
      S,
      align
    )
    
    mean_sims[m] <- agg$mean_similarity
  }
  
  return(mean_sims)
}
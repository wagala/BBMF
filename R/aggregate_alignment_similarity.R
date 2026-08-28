#' Summarize Similarities After Factor Alignment
#'
#' Extracts the similarities associated with an optimal factor
#' alignment and calculates their mean and total similarity.
#'
#' @param S A square similarity matrix used to obtain the factor
#'   alignment.
#' @param alignment An alignment object returned by
#'   \code{factor_alignment_hungarian()}.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{pair_similarities}}{
#'     Similarities for each matched pair of latent factors.
#'   }
#'   \item{\code{mean_similarity}}{
#'     Mean similarity across the matched factor pairs.
#'   }
#'   \item{\code{sum_similarity}}{
#'     Sum of the similarities across the matched factor pairs.
#'   }
#' }
#'
#' @seealso
#' \code{\link{factor_alignment_hungarian}}
#'
#' @export
#'
aggregate_alignment_similarity <- function(S, alignment) {
  
  assign_vec <- as.integer(
    alignment$assignment
  )
  
  K <- length(assign_vec)
  
  # Similarities corresponding to the matched factor pairs
  pair_sims <- S[
    cbind(
      seq_len(K),
      assign_vec
    )
  ]
  
  return(
    list(
      pair_similarities = pair_sims,
      mean_similarity = mean(pair_sims),
      sum_similarity = sum(pair_sims)
    )
  )
}
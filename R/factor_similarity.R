#' Compare and Align Latent Factor Matrices
#'
#' Computes pairwise Jaccard similarities between the rows of two binary
#' factor matrices, aligns the estimated factors to the reference factors
#' using the Hungarian algorithm, and summarizes the similarities of the
#' matched factor pairs.
#'
#' @param H_ref A binary reference factor matrix of dimension
#'   \eqn{R \times G}. In simulation studies, this may be the true
#'   factor matrix.
#' @param H_est A binary estimated factor matrix of dimension
#'   \eqn{R \times G}.
#' @param ref_names Optional character vector containing names for the
#'   reference factors. If \code{NULL}, names of the form
#'   \code{"F 1"}, \code{"F 2"}, ... are used.
#' @param est_names Optional character vector containing names for the
#'   estimated factors. If \code{NULL}, names of the form
#'   \code{"F 1"}, \code{"F 2"}, ... are used.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{similarity}}{
#'     The original pairwise Jaccard similarity matrix.
#'   }
#'   \item{\code{alignment}}{
#'     The optimal factor alignment returned by
#'     \code{factor_alignment_hungarian()}.
#'   }
#'   \item{\code{permutation}}{
#'     The optimal ordering of the estimated factors.
#'   }
#'   \item{\code{similarity_aligned}}{
#'     The Jaccard similarity matrix after reordering the estimated
#'     factors according to the optimal alignment.
#'   }
#'   \item{\code{pair_similarities}}{
#'     Jaccard similarities for the optimally matched factor pairs.
#'   }
#'   \item{\code{mean_similarity}}{
#'     Mean Jaccard similarity across the matched factor pairs.
#'   }
#'   \item{\code{sum_similarity}}{
#'     Sum of the Jaccard similarities across the matched factor pairs.
#'   }
#' }
#'
#' @details
#' Pairwise Jaccard similarities are first calculated between the rows
#' of \code{H_ref} and \code{H_est}. The Hungarian algorithm is then
#' used to determine the one-to-one assignment of estimated factors
#' that maximizes the total Jaccard similarity.
#'
#' This is useful because latent factor labels are arbitrary. For example,
#' factor 1 in an estimated decomposition may correspond to factor 3 in
#' the reference decomposition.
#'
#' @seealso
#' \code{\link{jaccard_rows}},
#' \code{\link{factor_alignment_hungarian}},
#' \code{\link{aggregate_alignment_similarity}}
#'
#' @export
#'
factor_similarity <- function(
    H_ref,
    H_est,
    ref_names = NULL,
    est_names = NULL) {
  
  H_ref <- as.matrix(H_ref)
  H_est <- as.matrix(H_est)
  
  if (ncol(H_ref) != ncol(H_est)) {
    stop("H_ref and H_est must have the same number of columns.")
  }
  
  if (nrow(H_ref) != nrow(H_est)) {
    stop(
      "factor_similarity() currently requires the same number of factors."
    )
  }
  
  # ---------------------------------------------------------------
  # Pairwise Jaccard similarities
  # ---------------------------------------------------------------
  
  S <- jaccard_rows(
    H_ref,
    H_est
  )
  
  # ---------------------------------------------------------------
  # Hungarian factor alignment
  # ---------------------------------------------------------------
  
  alignment <- factor_alignment_hungarian(
    S
  )
  
  # ---------------------------------------------------------------
  # Aggregate similarity
  # ---------------------------------------------------------------
  
  aggregate <- aggregate_alignment_similarity(
    S,
    alignment
  )
  
  permutation <- as.integer(
    alignment$assignment
  )
  
  # ---------------------------------------------------------------
  # Add factor names
  # ---------------------------------------------------------------
  
  if (is.null(ref_names)) {
    ref_names <- paste0(
      "F ",
      seq_len(nrow(H_ref))
    )
  }
  
  if (is.null(est_names)) {
    est_names <- paste0(
      "F ",
      seq_len(nrow(H_est))
    )
  }
  
  rownames(S) <- ref_names
  colnames(S) <- est_names
  
  # ---------------------------------------------------------------
  # Reorder estimated factors
  # ---------------------------------------------------------------
  
  S_aligned <- S[
    ,
    permutation,
    drop = FALSE
  ]
  
  return(
    list(
      similarity = S,
      alignment = alignment,
      permutation = permutation,
      similarity_aligned = S_aligned,
      pair_similarities = aggregate$pair_similarities,
      mean_similarity = aggregate$mean_similarity,
      sum_similarity = aggregate$sum_similarity
    )
  )
}
#' Align Posterior Samples of H to a Reference Matrix
#'
#' Aligns the rows of each posterior sample of the BBMF factor-feature
#' matrix H to a reference factor matrix using Jaccard similarity and
#' the Hungarian algorithm.
#'
#' @param H_ref A binary reference matrix of dimension
#'   \eqn{R \times G}.
#' @param H_samples A three-dimensional binary array of dimension
#'   \eqn{R \times G \times T}.
#'
#' @return A three-dimensional array with the same dimensions as
#'   \code{H_samples}, with factors aligned to \code{H_ref}.
#'
#' @details
#' For each posterior sample, pairwise Jaccard similarities are computed
#' between the rows of \code{H_ref} and the sampled H matrix. The
#' Hungarian algorithm is then used to determine the optimal factor
#' permutation.
#'
#' @seealso
#' \code{\link{jaccard_rows}},
#' \code{\link{factor_alignment_hungarian}}
#'
#' @export
#'
align_H_samples_bbmf <- function(
    H_ref,
    H_samples) {
  
  H_ref <- as.matrix(H_ref)
  
  if (length(dim(H_samples)) != 3L) {
    stop(
      "H_samples must be a three-dimensional array."
    )
  }
  
  if (nrow(H_ref) != dim(H_samples)[1]) {
    stop(
      "H_ref and H_samples must have the same number of factors."
    )
  }
  
  if (ncol(H_ref) != dim(H_samples)[2]) {
    stop(
      "H_ref and H_samples must have the same number of features."
    )
  }
  
  R <- dim(H_samples)[1]
  G <- dim(H_samples)[2]
  T <- dim(H_samples)[3]
  
  H_aligned <- array(
    0L,
    dim = c(
      R,
      G,
      T
    )
  )
  
  for (t in seq_len(T)) {
    
    H_t <- H_samples[
      ,
      ,
      t
    ]
    
    S <- jaccard_rows(
      H_ref,
      H_t
    )
    
    alignment <- factor_alignment_hungarian(
      S
    )
    
    permutation <- as.integer(
      alignment$assignment
    )
    
    H_aligned[
      ,
      ,
      t
    ] <- H_t[
      permutation,
      ,
      drop = FALSE
    ]
  }
  
  dimnames(H_aligned) <- dimnames(H_samples)
  
  return(H_aligned)
}
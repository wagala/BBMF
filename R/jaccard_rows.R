#' Compute Jaccard Similarities Between Matrix Rows
#'
#' Computes all pairwise Jaccard similarities between the rows of two
#' binary matrices having the same number of columns.
#'
#' @param X1 A binary matrix with features represented by columns.
#' @param X2 A binary matrix with the same number of columns as
#'   \code{X1}.
#'
#' @return A numeric matrix of dimension
#'   \code{nrow(X1)} by \code{nrow(X2)}. Entry \eqn{(i,j)} contains
#'   the Jaccard similarity between row \eqn{i} of \code{X1} and row
#'   \eqn{j} of \code{X2}.
#'
#' @details
#' For binary vectors \eqn{x} and \eqn{y}, the Jaccard similarity is
#'
#' \deqn{
#' J(x,y) =
#' \frac{|x \cap y|}{|x \cup y|}.
#' }
#'
#' If both vectors contain only zeros, their union is empty and the
#' similarity is defined here to be zero.
#'
#' This function is particularly useful for comparing rows of two
#' factor-feature matrices \code{H}.
#'
#' @export
#'
jaccard_rows <- function(X1, X2) {
  
  if (ncol(X1) != ncol(X2)) {
    stop(
      "X1 and X2 must have the same number of columns for row-wise comparison"
    )
  }
  
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  
  S <- matrix(
    0,
    nrow = n1,
    ncol = n2
  )
  
  for (i in seq_len(n1)) {
    
    for (j in seq_len(n2)) {
      
      x <- X1[i, ]
      y <- X2[j, ]
      
      inter <- sum(pmin(x, y))
      union <- sum(pmax(x, y))
      
      S[i, j] <- ifelse(
        union == 0,
        0,
        inter / union
      )
    }
  }
  
  return(S)
}
#' Compute Jaccard Similarities Between Matrix Columns
#'
#' Computes all pairwise Jaccard similarities between the columns of
#' two binary matrices having the same number of rows.
#'
#' @param X1 A binary matrix with samples represented by rows.
#' @param X2 A binary matrix with the same number of rows as
#'   \code{X1}.
#'
#' @return A numeric matrix of dimension
#'   \code{ncol(X1)} by \code{ncol(X2)}. Entry \eqn{(i,j)} contains
#'   the Jaccard similarity between column \eqn{i} of \code{X1} and
#'   column \eqn{j} of \code{X2}.
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
#' This function is particularly useful for comparing columns of two
#' sample-factor matrices \code{W}.
#'
#' @export
#'
jaccard_cols <- function(X1, X2) {
  
  if (nrow(X1) != nrow(X2)) {
    stop(
      "X1 and X2 must have the same number of rows for column-wise comparison"
    )
  }
  
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  
  p1 <- ncol(X1)
  p2 <- ncol(X2)
  
  S <- matrix(
    0,
    nrow = p1,
    ncol = p2
  )
  
  for (i in seq_len(p1)) {
    
    for (j in seq_len(p2)) {
      
      x <- X1[, i]
      y <- X2[, j]
      
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
#' Jaccard Similarity Between Two Binary Vectors
#'
#' Computes the Jaccard similarity coefficient (intersection over union)
#' between two binary vectors, with a defined value when both vectors
#' are entirely zero.
#'
#' @param x,y Binary (0/1 or logical) vectors of the same length.
#'
#' @return A numeric scalar in \eqn{[0,1]}.
#'
#' @details
#' \deqn{J(x,y) = \frac{|x \cap y|}{|x \cup y|} = \frac{\sum(x \wedge y)}{\sum(x \vee y)}.}
#' The ordinary definition is \eqn{0/0} (undefined) when both \code{x}
#' and \code{y} are all-zero -- e.g.\ two factors that touch no loci at
#' all. This is treated as a similarity of \code{1} (two "empty" factors
#' are considered to agree perfectly), rather than propagating a
#' \code{NaN} into any downstream alignment computation
#' (\code{\link{compute_afas}}).
#'
#' @seealso \code{\link{compute_afas}}, which uses this as its
#' pairwise similarity measure before Hungarian alignment.
#'
#' @export
jaccard <- function(x, y) {
  if (sum(x) == 0 && sum(y) == 0) return(1)
  sum(x & y) / sum(x | y)
}

#' Compute the Log-Sum-Exp of a Numeric Vector
#'
#' Computes the logarithm of the sum of exponentials of a numeric
#' vector using a numerically stable formulation.
#'
#' @param vec A numeric vector, typically containing values on the
#'   log scale.
#'
#' @return A numeric scalar equal to
#'   \eqn{\log\left(\sum_i \exp(vec_i)\right)}.
#'
#' @details
#' The function first sorts the values in decreasing order and uses
#' the largest value as a reference. The remaining terms are evaluated
#' relative to this value to improve numerical stability and reduce
#' the risk of overflow.
#'
sumLog <- function(vec) {
  ord <- sort(vec, decreasing = TRUE)
  s <- ord[1]
  s <- s + sum(log1p(exp(ord[-1] - s)))
  return(s)
}
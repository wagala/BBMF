#' Sample from a Truncated Beta Distribution
#'
#' Draws exact samples from a \eqn{\mathrm{Beta}(\mathrm{shape1},
#' \mathrm{shape2})} distribution restricted to the interval
#' \code{(lower, upper)}, via the inverse-CDF method.
#'
#' @param n Number of draws.
#' @param shape1,shape2 Shape parameters of the (untruncated) Beta
#'   distribution.
#' @param lower,upper Lower and upper truncation bounds, within
#'   \eqn{[0,1]}.
#'
#' @return A numeric vector of length \code{n}, with every value
#'   strictly within \code{[lower, upper]}.
#'
#' @details
#' If \eqn{U \sim \mathrm{Unif}(F(\mathrm{lower}), F(\mathrm{upper}))}
#' and \eqn{F} is the \eqn{\mathrm{Beta}(\mathrm{shape1},
#' \mathrm{shape2})} CDF, then \eqn{F^{-1}(U)} has exactly the
#' truncated distribution. This is exact, not approximate, and is used
#' to enforce the identifiability constraint \eqn{p_{11} > p_{10}}
#' (see \code{\link{update_noise}}).
#'
#' @export
rtrunc_beta <- function(n, shape1, shape2, lower, upper) {

  lower <- max(lower, 0)
  upper <- min(upper, 1)

  u <- stats::runif(n, stats::pbeta(lower, shape1, shape2), stats::pbeta(upper, shape1, shape2))
  x <- stats::qbeta(u, shape1, shape2)

  return(x)
}

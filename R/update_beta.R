#' Update the Column-Specific Bernoulli Probabilities Beta
#'
#' Samples the column-specific Bernoulli probabilities \code{beta}
#' from their full conditional Beta distributions given the current
#' binary factor matrix \code{H}.
#'
#' @param H A binary matrix of dimension \eqn{R \times G} containing
#'   the current factor-specific feature patterns.
#' @param b1 First shape parameter of the Beta prior for \code{beta}.
#' @param b2 Second shape parameter of the Beta prior for \code{beta}.
#'
#' @return A numeric vector of length \eqn{G} containing the updated
#'   column-specific probabilities \code{beta}.
#'
#' @details
#' For each column \eqn{g}, the number of active entries is given by
#' \code{sum(H[, g])}, while the number of inactive entries is given by
#' \code{sum(1 - H[, g])}. Each \eqn{\beta_g} is then sampled from its
#' Beta full conditional distribution.
#'
update_beta<-function(H,b1,b2){
  G<-dim(H)[2]
  beta <- rep(NA, G)
  #update beta
  for(g in 1:G){
    beta[g] <- stats::rbeta(1,shape1=sum(H[,g])+b1, shape2=sum(1-H[,g])+b2)
  }
  return(beta)
}
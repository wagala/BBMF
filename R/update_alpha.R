#' Update the Row-Specific Bernoulli Probabilities Alpha
#'
#' Samples the row-specific Bernoulli probabilities \code{alpha} from
#' their full conditional Beta distributions given the current binary
#' factor-activation matrix \code{W}.
#'
#' @param W A binary matrix of dimension \eqn{K \times R} containing
#'   the current factor activations.
#' @param a1 First shape parameter of the Beta prior for \code{alpha}.
#' @param a2 Second shape parameter of the Beta prior for \code{alpha}.
#'
#' @return A numeric vector of length \eqn{K} containing the updated
#'   row-specific probabilities \code{alpha}.
#'
#' @details
#' For each row \eqn{k}, the number of active factors is given by
#' \code{sum(W[k, ])}, while the number of inactive factors is given by
#' \code{sum(1 - W[k, ])}. Each \eqn{\alpha_k} is then sampled from its
#' Beta full conditional distribution.
#'
update_alpha<-function(W,a1,a2){
  K<-dim(W)[1]
  alpha <- rep(NA, K)
  #update alpha
  for (k in 1:K) {
    alpha[k] <- stats::rbeta(1, shape1=sum(W[k,])+a1, 
                      shape2= sum(1-W[k,])+a2)
  }
  return(alpha)
}
#' Sample an Auxiliary Bernoulli Pair Given Their AND
#'
#' Given a binary matrix whose entries are the logical AND of two
#' hidden Bernoulli variables, draws the corresponding pair of hidden
#' variables from their exact conditional distribution. This is the
#' core computational device used to restore conjugacy for the
#' factor-inclusion probabilities \code{phi} (see
#' \code{\link{update_phi}}).
#'
#' @param Wmat A binary matrix of observed AND-results.
#' @param Pmat A matrix of the same dimension as \code{Wmat}, giving
#'   \eqn{P(A=1)} for each entry.
#' @param Qmat A matrix of the same dimension as \code{Wmat}, giving
#'   \eqn{P(B=1)} for each entry.
#'
#' @return A list with elements \code{A} and \code{B}, each a matrix
#'   of the same dimension as \code{Wmat}.
#'
#' @details
#' If \eqn{A \sim \mathrm{Bernoulli}(p)} and
#' \eqn{B \sim \mathrm{Bernoulli}(q)} independently, and
#' \eqn{W = A \wedge B}, then \eqn{W \sim \mathrm{Bernoulli}(pq)}.
#' Conditional on \eqn{W=1}, \eqn{(A,B)=(1,1)} with probability 1.
#' Conditional on \eqn{W=0}, \eqn{(A,B)} takes the values
#' \eqn{(0,0)}, \eqn{(0,1)}, \eqn{(1,0)} with probability
#' proportional to \eqn{(1-p)(1-q)}, \eqn{(1-p)q}, \eqn{p(1-q)}
#' respectively. This function draws exactly from that conditional
#' distribution, entrywise across the whole matrix at once.
#'
#' @export
sample_aux_pair <- function(Wmat, Pmat, Qmat) {

  nr <- nrow(Wmat)
  nc <- ncol(Wmat)
  A <- matrix(0L, nr, nc)
  B <- matrix(0L, nr, nc)

  idx1 <- which(Wmat == 1)
  A[idx1] <- 1L
  B[idx1] <- 1L

  idx0 <- which(Wmat == 0)
  if (length(idx0) > 0) {
    p <- Pmat[idx0]
    q <- Qmat[idx0]
    denom  <- 1 - p * q
    prob00 <- (1 - p) * (1 - q) / denom
    prob01 <- (1 - p) * q / denom

    u <- stats::runif(length(idx0))
    is01 <- u > prob00 & u <= (prob00 + prob01)
    is10 <- u > (prob00 + prob01)

    A[idx0[is01]] <- 0L; B[idx0[is01]] <- 1L
    A[idx0[is10]] <- 1L; B[idx0[is10]] <- 0L
  }

  return(list(A = A, B = B))
}

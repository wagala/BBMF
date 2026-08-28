#' Reconstruct Binary Matrices Across One BBMF Chain
#'
#' Reconstructs the binary matrix \eqn{X} for every posterior sample
#' of \eqn{W} and \eqn{H} from a single BBMF MCMC chain.
#'
#' @param W_samples A three-dimensional binary array of dimension
#'   \eqn{K \times R \times T}, containing posterior samples of
#'   the sample-factor matrix \eqn{W}.
#' @param H_samples A three-dimensional binary array of dimension
#'   \eqn{R \times G \times T}, containing posterior samples of
#'   the factor-feature matrix \eqn{H}.
#'
#' @return A three-dimensional binary array of dimension
#'   \eqn{K \times G \times T}. Slice \code{t} contains the Boolean
#'   reconstruction obtained from the corresponding posterior samples
#'   of \eqn{W} and \eqn{H}.
#'
#' @details
#' For each posterior iteration \eqn{t}, the reconstruction is
#'
#' \deqn{
#' \widetilde{X}^{(t)}
#' =
#' W^{(t)} \circ H^{(t)},
#' }
#'
#' where \eqn{\circ} denotes the Boolean matrix product implemented
#' by \code{bool_prod()}.
#'
#' @seealso \code{\link{bool_prod}}
#'
#' @examples
#' \dontrun{
#' X_samples <- reconstruct_X_one_chain_bbmf(
#'   W_samples = Ws[[1]],
#'   H_samples = Hs[[1]]
#' )
#' }
#'
#' @export
#'
reconstruct_X_one_chain_bbmf <- function(
    W_samples,
    H_samples) {
  
  if (length(dim(W_samples)) != 3L) {
    stop(
      "W_samples must be a three-dimensional array."
    )
  }
  
  if (length(dim(H_samples)) != 3L) {
    stop(
      "H_samples must be a three-dimensional array."
    )
  }
  
  K <- dim(W_samples)[1]
  R <- dim(W_samples)[2]
  T <- dim(W_samples)[3]
  
  R_H <- dim(H_samples)[1]
  G <- dim(H_samples)[2]
  T_H <- dim(H_samples)[3]
  
  if (R != R_H) {
    stop(
      "The factor dimensions of W_samples and H_samples differ."
    )
  }
  
  if (T != T_H) {
    stop(
      "W_samples and H_samples contain different numbers of iterations."
    )
  }
  
  X_reconstructed <- array(
    0L,
    dim = c(
      K,
      G,
      T
    )
  )
  
  for (t in seq_len(T)) {
    
    W_t <- W_samples[
      ,
      ,
      t
    ]
    
    H_t <- H_samples[
      ,
      ,
      t
    ]
    
    X_reconstructed[
      ,
      ,
      t
    ] <- bool_prod(
      W = W_t,
      H = H_t
    )
  }
  
  return(X_reconstructed)
}
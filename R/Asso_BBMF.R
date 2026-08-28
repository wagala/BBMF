#' Initialize BBMF Using the Asso Algorithm
#'
#' Uses the Asso Boolean matrix factorization algorithm to obtain
#' initial estimates of the binary factor matrices \code{W} and
#' \code{H}. The reconstructed binary matrix is then computed using
#' the Boolean matrix product implemented by \code{bool_prod()}.
#'
#' @param X A binary matrix of dimension \eqn{K \times G}, where rows
#'   correspond to samples and columns correspond to features.
#' @param R A positive integer specifying the number of latent factors.
#' @param seed An optional integer used to set the random seed before
#'   fitting the Asso model. The default is \code{NULL}, in which case
#'   the current random-number-generator state is used.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{W_asso}}{
#'     A binary matrix of dimension \eqn{K \times R} containing the
#'     sample-specific factor activations estimated by Asso.
#'   }
#'   \item{\code{H_asso}}{
#'     A binary matrix of dimension \eqn{R \times G} containing the
#'     factor-specific feature patterns estimated by Asso.
#'   }
#'   \item{\code{X_asso}}{
#'     The reconstructed binary matrix obtained from the Boolean
#'     product of \code{W_asso} and \code{H_asso}.
#'   }
#' }
#'
#' @details
#' The input matrix \code{X} is first converted to a logical binary
#' matrix and passed to \code{rBMF::Asso_approximate()}.
#'
#' The Asso algorithm is fitted using a threshold of 0.5, an
#' over-coverage penalty of 1, and a coverage bonus of 1.
#'
#' After estimating the factor matrices, the fitted binary matrix is
#' reconstructed as
#'
#' \deqn{
#' X = W \circ H,
#' }
#'
#' where \eqn{\circ} denotes the Boolean matrix product. Logical AND
#' combines entries within a latent factor, while logical OR combines
#' contributions across latent factors.
#'
#' The package \code{rBMF} must be installed to use this function.
#'
#' @seealso
#' \code{\link{bool_prod}}
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#'
#' X <- matrix(
#'   rbinom(100, size = 1, prob = 0.3),
#'   nrow = 10,
#'   ncol = 10
#' )
#'
#' fit <- Asso_BBMF(
#'   X = X,
#'   R = 2,
#'   seed = 123
#' )
#'
#' fit$W_asso
#' fit$H_asso
#' fit$X_asso
#' }
#'
#' @export
#'
Asso_BBMF <- function(X, R, seed = NULL) {
  
  if (!requireNamespace("rBMF", quietly = TRUE)) {
    stop("Package 'rBMF' is required to use Asso_BBMF().")
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  Xb <- X == 1
  
  Res_Asso <- rBMF::Asso_approximate(
    Xb,
    R,
    list(
      threshold = 0.5,
      penalty_overcovered = 1,
      bonus_covered = 1,
      verbose = 0
    )
  )
  
  W_asso <- as.matrix(Res_Asso$O) * 1
  H_asso <- as.matrix(Res_Asso$B) * 1
  
  X_asso <- bool_prod(W_asso, H_asso)
  
  return(
    list(
      W_asso = W_asso,
      H_asso = H_asso,
      X_asso = X_asso
    )
  )
}
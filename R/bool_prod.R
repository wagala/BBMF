#' Compute the Boolean Matrix Product
#'
#' Computes the Boolean product of two binary matrices \code{W} and
#' \code{H}. Contributions from the latent factors are combined using
#' the logical OR operation.
#'
#' @param W A binary matrix of dimension \eqn{K \times R}.
#' @param H A binary matrix of dimension \eqn{R \times G}.
#'
#' @return An integer matrix of dimension \eqn{K \times G} containing
#'   the Boolean product of \code{W} and \code{H}.
#'
#' @details
#' For each latent factor \eqn{r}, the function combines column
#' \eqn{r} of \code{W} with row \eqn{r} of \code{H}. The resulting
#' contributions are combined across factors using logical OR.
#'
#' The returned matrix has integer storage mode, with entries equal
#' to 0 or 1.
#'
#' @export
#'
bool_prod <- function(W, H){
  K <- nrow(W); G <- ncol(H)
  R <- ncol(W)
  X <- matrix(0L, K, G)
  for(r in seq_len(R)){
    X <- X | (W[, r, drop = FALSE] %*% H[r, , drop = FALSE])
  }
  storage.mode(X) <- "integer"
  X
}
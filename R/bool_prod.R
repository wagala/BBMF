#' Boolean Matrix Product
#'
#' Computes the Boolean (disjunctive) matrix product of two binary
#' matrices \code{W} and \code{H}, i.e. the reconstructed matrix
#' \eqn{Z = W \circ H}, where entry \eqn{Z_{kg} = 1} if
#' and only if at least one factor \eqn{r} satisfies
#' \eqn{W_{kr} = H_{rg} = 1}.
#'
#' @param W A binary matrix of dimension \eqn{K \times R} containing
#'   the current factor activations.
#' @param H A binary matrix of dimension \eqn{R \times G} containing
#'   the current factor-specific feature patterns.
#'
#' @return A binary (0/1) matrix of dimension \eqn{K \times G}
#'   containing the reconstructed entries \eqn{Z_{kg}}.
#'
#' @details
#' The Boolean product is computed efficiently via ordinary
#' (real-valued) matrix multiplication: since \eqn{W_{kr}} and
#' \eqn{H_{rg}} are each 0 or 1, the ordinary matrix product
#' \eqn{(WH)_{kg}} counts the number of factors \eqn{r} for which
#' \eqn{W_{kr} = H_{rg} = 1} simultaneously. Thresholding this count
#' at zero recovers the Boolean OR exactly.
#'
#' @export
bool_prod <- function(W, H) {

  raw_count <- W %*% H
  Z <- (raw_count > 0) * 1

  return(Z)
}

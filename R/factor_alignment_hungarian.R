#' Align Latent Factors Using the Hungarian Algorithm
#'
#' Finds the one-to-one matching between two sets of latent factors
#' that maximizes their total similarity.
#'
#' @param S A square numeric similarity matrix. Entry \eqn{S_{ij}}
#'   represents the similarity between factor \eqn{i} in the first
#'   factorization and factor \eqn{j} in the second factorization.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{assignment}}{
#'     The optimal assignment returned by the Hungarian algorithm.
#'   }
#'   \item{\code{permutation}}{
#'     An integer vector giving the factor ordering in the second
#'     factorization corresponding to the factors in the first.
#'   }
#'   \item{\code{matching}}{
#'     A two-column matrix identifying the matched factor pairs.
#'   }
#'   \item{\code{total_similarity}}{
#'     The sum of the similarities across all optimally matched pairs.
#'   }
#' }
#'
#' @details
#' The Hungarian algorithm solves a minimum-cost assignment problem.
#' The similarity matrix is therefore converted to a cost matrix using
#'
#' \deqn{
#' C_{ij} = \max(S) - S_{ij}.
#' }
#'
#' Minimizing the resulting total cost is equivalent to maximizing
#' the total similarity between matched latent factors.
#'
#' This function currently requires the two factorizations to contain
#' the same number of latent factors.
#'
#' @seealso
#' \code{\link{jaccard_rows}},
#' \code{\link{jaccard_cols}},
#' \code{\link{aggregate_alignment_similarity}}
#'
#' @export
#'
factor_alignment_hungarian <- function(S) {
  
  if (nrow(S) != ncol(S)) {
    stop(
      "Currently implemented for square similarity matrices (same number of factors)."
    )
  }
  
  # Convert similarity to cost
  maxS <- max(S)
  cost <- maxS - S
  
  # Hungarian algorithm
  assignment <- clue::solve_LSAP(cost)
  
  assignment_integer <- as.integer(
    assignment
  )
  
  return(
    list(
      assignment = assignment,
      permutation = assignment_integer,
      matching = cbind(
        model1 = seq_len(nrow(S)),
        model2 = assignment_integer
      ),
      total_similarity = sum(
        S[
          cbind(
            seq_len(nrow(S)),
            assignment_integer
          )
        ]
      )
    )
  )
}
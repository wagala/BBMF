#' Align Rows of One Matrix to Another
#'
#' Aligns the rows of matrix \code{B} to the rows of a reference matrix
#' \code{A} by minimizing the total row-wise absolute difference. The
#' optimal one-to-one assignment is obtained using the Hungarian
#' algorithm.
#'
#' @param A A reference matrix with \eqn{m} rows and \eqn{p} columns.
#' @param B A matrix with the same number of columns as \code{A} and
#'   \eqn{n} rows, where typically \eqn{n \ge m}.
#' @param verbose Logical. If \code{TRUE}, prints a summary of the
#'   optimal row alignment. Default is \code{TRUE}.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{A}}{
#'     The reference matrix \code{A}.
#'   }
#'   \item{\code{B_ordered}}{
#'     Matrix \code{B} with its rows reordered so that the rows
#'     optimally matched to \code{A} appear first, followed by any
#'     unmatched rows.
#'   }
#'   \item{\code{cost}}{
#'     The minimum total alignment cost.
#'   }
#'   \item{\code{matched_B_rows}}{
#'     Integer indices of the rows of \code{B} matched to the rows
#'     of \code{A}.
#'   }
#'   \item{\code{full_B_order}}{
#'     Integer vector giving the complete reordered row indices of
#'     \code{B}.
#'   }
#'   \item{\code{matched_B_row_names}}{
#'     Names of the rows of \code{B} matched to \code{A}.
#'   }
#' }
#'
#' @details
#' The function considers all subsets of \code{B} containing the same
#' number of rows as \code{A}. For each subset, a cost matrix is
#' constructed using the absolute difference between pairs of rows:
#'
#' \deqn{
#' C_{ij}
#' =
#' \sum_g
#' |A_{ig} - B_{jg}|.
#' }
#'
#' The Hungarian algorithm is then used to find the one-to-one
#' assignment that minimizes the total cost.
#'
#' If \code{B} contains more rows than \code{A}, the optimally matched
#' rows are placed first in \code{B_ordered}, while the remaining
#' unmatched rows retain their original relative order.
#'
#' For binary matrices, the absolute-difference cost corresponds to
#' the number of disagreements between two rows.
#'
#' This function is useful for aligning factor-feature matrices when
#' two decompositions contain different numbers of latent factors.
#'
#' @seealso
#' \code{\link{align_cols_to_matrix}},
#' \code{\link{factor_alignment_hungarian}},
#' \code{\link{jaccard_rows}}
#'
#' @export
#'
align_rows_to_matrix <- function(A, B, verbose = TRUE) {
  
  A <- as.matrix(A)
  B <- as.matrix(B)
  
  # Check dimensions
  if (ncol(A) != ncol(B)) {
    stop("A and B must have the same number of columns.")
  }
  
  # Add default names if missing
  if (is.null(rownames(A))) {
    rownames(A) <- paste0("R", 1:nrow(A))
  }
  
  if (is.null(rownames(B))) {
    rownames(B) <- paste0("R", 1:nrow(B))
  }
  
  if (is.null(colnames(A)) || is.null(colnames(B))) {
    colnames(A) <- colnames(B) <- paste0("C", 1:ncol(A))
  }
  
  m <- nrow(A)
  n <- nrow(B)
  
  # Generate all subsets of B rows of size m
  row_subsets <- utils::combn(
    1:n,
    m,
    simplify = FALSE
  )
  
  best_cost <- Inf
  best_matching <- NULL
  best_B_rows <- NULL
  
  for (rows in row_subsets) {
    
    B_sub <- B[
      rows,
      ,
      drop = FALSE
    ]
    
    # Cost: absolute difference between A and B_sub rows
    cost_mat <- outer(
      1:m,
      1:m,
      Vectorize(
        function(i, j) {
          sum(abs(A[i, ] - B_sub[j, ]))
        }
      )
    )
    
    # Hungarian assignment
    assignment <- clue::solve_LSAP(
      cost_mat
    )
    
    total_cost <- sum(
      cost_mat[
        cbind(
          1:m,
          assignment
        )
      ]
    )
    
    if (total_cost < best_cost) {
      best_cost <- total_cost
      best_matching <- assignment
      best_B_rows <- rows
    }
  }
  
  # Reorder matched rows according to the optimal assignment
  matched_B_order <- best_B_rows[
    best_matching
  ]
  
  # Rows of B not matched to A
  nonmatched_B <- setdiff(
    seq_len(n),
    matched_B_order
  )
  
  # Complete row ordering
  full_B_order <- c(
    matched_B_order,
    nonmatched_B
  )
  
  B_ordered <- B[
    full_B_order,
    ,
    drop = FALSE
  ]
  
  rownames(B_ordered) <- c(
    rownames(B)[matched_B_order],
    rownames(B)[nonmatched_B]
  )
  
  colnames(B_ordered) <- colnames(B)
  
  # Prepare output
  out <- list(
    A = A,
    B_ordered = B_ordered,
    cost = best_cost,
    matched_B_rows = matched_B_order,
    full_B_order = full_B_order,
    matched_B_row_names = rownames(B)[matched_B_order]
  )
  
  if (verbose) {
    
    cat("=== Row Alignment Summary ===\n")
    
    cat(
      "Assignment cost:",
      out$cost,
      "\n"
    )
    
    cat(
      "Matched B rows (indices):",
      paste(
        out$matched_B_rows,
        collapse = ", "
      ),
      "\n"
    )
    
    cat(
      "Matched B row names:",
      paste(
        out$matched_B_row_names,
        collapse = ", "
      ),
      "\n\n"
    )
  }
  
  return(out)
}
#' Align Columns of One Binary Matrix to Another
#'
#' Aligns the columns of matrix \code{B} to the columns of a reference
#' matrix \code{A} by minimizing the total column-wise absolute
#' difference. The optimal one-to-one assignment is obtained using
#' the Hungarian algorithm.
#'
#' @param A A reference matrix with \eqn{n} rows and \eqn{p} columns.
#' @param B A matrix with the same number of rows as \code{A} and
#'   \eqn{q} columns, where typically \eqn{q \ge p}.
#' @param verbose Logical. If \code{TRUE}, prints a summary of the
#'   optimal column alignment. Default is \code{TRUE}.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{A}}{
#'     The reference matrix \code{A}.
#'   }
#'   \item{\code{B_ordered}}{
#'     Matrix \code{B} with its columns reordered so that the columns
#'     optimally matched to \code{A} appear first, followed by any
#'     unmatched columns.
#'   }
#'   \item{\code{cost}}{
#'     The minimum total alignment cost.
#'   }
#'   \item{\code{matched_B_cols}}{
#'     Integer indices of the columns of \code{B} matched to the
#'     columns of \code{A}.
#'   }
#'   \item{\code{full_B_order}}{
#'     Integer vector giving the complete reordered column indices
#'     of \code{B}.
#'   }
#'   \item{\code{matched_B_col_names}}{
#'     Names of the columns of \code{B} matched to \code{A}.
#'   }
#' }
#'
#' @details
#' The function considers subsets of \code{B} containing the same
#' number of columns as \code{A}. For each subset, a cost matrix is
#' constructed using the absolute difference between pairs of columns:
#'
#' \deqn{
#' C_{ij}
#' =
#' \sum_k
#' |A_{ki} - B_{kj}|.
#' }
#'
#' The Hungarian algorithm is then used to find the one-to-one
#' assignment minimizing the total cost.
#'
#' If \code{B} contains more columns than \code{A}, the optimally
#' matched columns are placed first in \code{B_ordered}, while the
#' remaining unmatched columns retain their original relative order.
#'
#' For binary matrices, the absolute-difference cost corresponds to
#' the number of disagreements between two columns.
#'
#' @seealso
#' \code{\link{align_rows_to_matrix}},
#' \code{\link{factor_alignment_hungarian}},
#' \code{\link{jaccard_cols}}
#'
#' @export
#'
align_cols_to_matrix <- function(A, B, verbose = TRUE) {
  
  A <- as.matrix(A)
  B <- as.matrix(B)
  
  # Check dimensions
  if (nrow(A) != nrow(B)) {
    stop("A and B must have the same number of rows.")
  }
  
  # Add default names if missing
  if (is.null(colnames(A))) {
    colnames(A) <- paste0("C", 1:ncol(A))
  }
  
  if (is.null(colnames(B))) {
    colnames(B) <- paste0("C", 1:ncol(B))
  }
  
  if (is.null(rownames(A)) || is.null(rownames(B))) {
    rownames(A) <- rownames(B) <- paste0("R", 1:nrow(A))
  }
  
  p <- ncol(A)
  q <- ncol(B)
  
  # Generate all subsets of B columns of size p
  col_subsets <- utils::combn(
    1:q,
    p,
    simplify = FALSE
  )
  
  best_cost <- Inf
  best_matching <- NULL
  best_B_cols <- NULL
  
  for (cols in col_subsets) {
    
    B_sub <- B[, cols, drop = FALSE]
    
    # Cost: absolute difference between A and B_sub columns
    cost_mat <- outer(
      1:p,
      1:p,
      Vectorize(
        function(i, j) {
          sum(abs(A[, i] - B_sub[, j]))
        }
      )
    )
    
    assignment <- clue::solve_LSAP(
      cost_mat
    )
    
    total_cost <- sum(
      cost_mat[
        cbind(
          1:p,
          assignment
        )
      ]
    )
    
    if (total_cost < best_cost) {
      best_cost <- total_cost
      best_matching <- assignment
      best_B_cols <- cols
    }
  }
  
  matched_B_order <- best_B_cols[
    best_matching
  ]
  
  nonmatched_B <- setdiff(
    seq_len(q),
    matched_B_order
  )
  
  full_B_order <- c(
    matched_B_order,
    nonmatched_B
  )
  
  B_ordered <- B[
    ,
    full_B_order,
    drop = FALSE
  ]
  
  colnames(B_ordered) <- c(
    colnames(B)[matched_B_order],
    colnames(B)[nonmatched_B]
  )
  
  rownames(B_ordered) <- rownames(B)
  
  # Prepare output list
  out <- list(
    A = A,
    B_ordered = B_ordered,
    cost = best_cost,
    matched_B_cols = matched_B_order,
    full_B_order = full_B_order,
    matched_B_col_names = colnames(B)[matched_B_order]
  )
  
  if (verbose) {
    
    cat("=== Column Alignment Summary ===\n")
    
    cat(
      "Assignment cost:",
      out$cost,
      "\n"
    )
    
    cat(
      "Matched B columns (indices):",
      paste(
        out$matched_B_cols,
        collapse = ", "
      ),
      "\n"
    )
    
    cat(
      "Matched B column names:",
      paste(
        out$matched_B_col_names,
        collapse = ", "
      ),
      "\n\n"
    )
  }
  
  return(out)
}
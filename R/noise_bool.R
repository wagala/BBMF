#' Introduce Noise into a Boolean Matrix
#'
#' Randomly flips a specified proportion of entries in a binary matrix.
#'
#' @param mat A binary matrix containing entries equal to 0 or 1.
#' @param noise_level Numeric value between 0 and 1 specifying the
#'   proportion of matrix entries to flip.
#' @param seed Optional integer seed for reproducibility. If \code{NULL},
#'   the current random-number state is used.
#'
#' @return A binary matrix with the selected entries flipped.
#'
#' @details
#' The number of entries to flip is calculated as
#' \code{ceiling(noise_level * length(mat))}. The selected entries
#' are sampled without replacement. Entries equal to 0 are changed to 1,
#' and entries equal to 1 are changed to 0.
#'
#' @export
#' 
noise_bool <- function(mat, noise_level, seed = NULL) {
  # Check if noise_level is between 0 and 1
  if (noise_level < 0 || noise_level > 1) {
    stop("noise_level must be between 0 and 1")
  }
  
  # Determine the number of elements to flip
  num_elements <- length(mat)
  num_flips <- ceiling(noise_level * num_elements)
  
  # Optionally set seed for reproducibility
  if (!is.null(seed)) set.seed(seed)
  
  # Randomly select elements to flip
  flip_indices <- sample(seq_len(num_elements), num_flips, replace = FALSE)
  
  # Flip the selected elements
  mat[flip_indices] <- 1 - mat[flip_indices]
  
  return(mat)
}
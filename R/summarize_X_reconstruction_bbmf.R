#' Summarize Posterior BBMF Reconstructions
#'
#' Combines posterior Boolean reconstructions across MCMC chains and
#' computes entry-wise posterior reconstruction probabilities and
#' uncertainty scores.
#'
#' @param reconstructed_chains A list of three-dimensional binary arrays.
#'   Each array must have dimension \eqn{K \times G \times T}, where
#'   the third dimension indexes posterior samples from one MCMC chain.
#' @param X_truth A binary reference matrix of dimension
#'   \eqn{K \times G}. The matrix must have row and column names.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{reconstructed_chains}}{
#'     The reconstructed posterior arrays with dimension names added.
#'   }
#'   \item{\code{reconstructed_all}}{
#'     A single three-dimensional array containing posterior
#'     reconstructions from all chains.
#'   }
#'   \item{\code{probability}}{
#'     A \eqn{K \times G} matrix of posterior reconstruction
#'     probabilities.
#'   }
#'   \item{\code{uncertainty}}{
#'     A \eqn{K \times G} matrix of scaled posterior uncertainty scores.
#'   }
#' }
#'
#' @details
#' For matrix entry \eqn{(k,g)}, the posterior reconstruction probability
#' is estimated by
#'
#' \deqn{
#' \widehat{P}_{kg}
#' =
#' \frac{1}{M}
#' \sum_{m=1}^{M}
#' \widetilde{X}_{kg}^{(m)}.
#' }
#'
#' Posterior reconstruction uncertainty is calculated as
#'
#' \deqn{
#' U_{kg}
#' =
#' 4\widehat{P}_{kg}
#' (1-\widehat{P}_{kg}).
#' }
#'
#' The uncertainty score lies between zero and one. A value of zero
#' corresponds to complete posterior agreement, whereas a value of one
#' occurs when the posterior reconstruction probability is 0.5.
#'
#' @examples
#' \dontrun{
#' X_summary <- summarize_X_reconstruction_bbmf(
#'   reconstructed_chains = X_reconstructed,
#'   X_truth = X_true
#' )
#' }
#'
#' @export
#'
summarize_X_reconstruction_bbmf <- function(
    reconstructed_chains,
    X_truth) {
  
  if (!is.matrix(X_truth)) {
    stop(
      "X_truth must be a matrix."
    )
  }
  
  if (!all(X_truth %in% c(0, 1))) {
    stop(
      "X_truth must contain only zeros and ones."
    )
  }
  
  if (is.null(rownames(X_truth))) {
    stop(
      "X_truth must have row names."
    )
  }
  
  if (is.null(colnames(X_truth))) {
    stop(
      "X_truth must have column names."
    )
  }
  
  if (!is.list(reconstructed_chains)) {
    stop(
      "reconstructed_chains must be a list."
    )
  }
  
  if (length(reconstructed_chains) == 0L) {
    stop(
      "reconstructed_chains is empty."
    )
  }
  
  storage.mode(X_truth) <- "integer"
  
  sample_names <- rownames(X_truth)
  feature_names <- colnames(X_truth)
  
  K <- nrow(X_truth)
  G <- ncol(X_truth)
  
  
  ## Check reconstructed arrays
  for (chain_id in seq_along(reconstructed_chains)) {
    
    X_chain <- reconstructed_chains[[chain_id]]
    
    if (length(dim(X_chain)) != 3L) {
      stop(
        "Chain ",
        chain_id,
        " must be a three-dimensional array."
      )
    }
    
    if (dim(X_chain)[1] != K) {
      stop(
        "Chain ",
        chain_id,
        " has a different number of rows from X_truth."
      )
    }
    
    if (dim(X_chain)[2] != G) {
      stop(
        "Chain ",
        chain_id,
        " has a different number of columns from X_truth."
      )
    }
    
    if (!all(X_chain %in% c(0, 1))) {
      stop(
        "Chain ",
        chain_id,
        " contains values other than zero and one."
      )
    }
  }
  
  
  ## Add dimension names
  reconstructed_chains <- lapply(
    reconstructed_chains,
    function(X_array) {
      
      dimnames(X_array) <- list(
        Sample = sample_names,
        Feature = feature_names,
        Iteration = paste0(
          "Iteration_",
          seq_len(dim(X_array)[3])
        )
      )
      
      X_array
    }
  )
  
  
  ## Combine posterior samples across chains
  reconstructed_all <- abind::abind(
    reconstructed_chains,
    along = 3
  )
  
  dimnames(reconstructed_all) <- list(
    Sample = sample_names,
    Feature = feature_names,
    Posterior_sample = paste0(
      "Posterior_sample_",
      seq_len(dim(reconstructed_all)[3])
    )
  )
  
  
  ## Posterior reconstruction probability
  probability <- apply(
    reconstructed_all,
    MARGIN = c(1, 2),
    FUN = mean
  )
  
  rownames(probability) <- sample_names
  colnames(probability) <- feature_names
  
  
  ## Posterior uncertainty
  uncertainty <- 4 *
    probability *
    (1 - probability)
  
  rownames(uncertainty) <- sample_names
  colnames(uncertainty) <- feature_names
  
  
  return(
    list(
      reconstructed_chains = reconstructed_chains,
      reconstructed_all = reconstructed_all,
      probability = probability,
      uncertainty = uncertainty
    )
  )
}
#' Summarize Posterior BBMF Factor Loadings
#'
#' Combines aligned posterior samples of the BBMF factor-feature matrix
#' H across MCMC chains and computes posterior inclusion probabilities
#' and uncertainty scores.
#'
#' @param H_chains A list of aligned three-dimensional binary arrays.
#'   Each array must have dimension \eqn{R \times G \times T}.
#' @param H_ref A reference binary matrix of dimension
#'   \eqn{R \times G}, used to provide factor and feature names.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{H_all}}{
#'     Posterior samples of H combined across chains.
#'   }
#'   \item{\code{probability}}{
#'     Posterior inclusion probability for each factor-feature entry.
#'   }
#'   \item{\code{uncertainty}}{
#'     Scaled posterior uncertainty for each factor-feature entry.
#'   }
#' }
#'
#' @details
#' For factor \eqn{r} and feature \eqn{g}, the posterior inclusion
#' probability is
#'
#' \deqn{
#' \widehat{P}_{rg}
#' =
#' \frac{1}{M}
#' \sum_{m=1}^{M}
#' H_{rg}^{(m)}.
#' }
#'
#' Posterior uncertainty is calculated as
#'
#' \deqn{
#' U_{rg}
#' =
#' 4\widehat{P}_{rg}
#' (1-\widehat{P}_{rg}).
#' }
#'
#' @export
#'
summarize_H_posterior_bbmf <- function(
    H_chains,
    H_ref) {
  
  if (!is.list(H_chains)) {
    stop(
      "H_chains must be a list."
    )
  }
  
  if (length(H_chains) == 0L) {
    stop(
      "H_chains is empty."
    )
  }
  
  H_ref <- as.matrix(H_ref)
  
  R <- nrow(H_ref)
  G <- ncol(H_ref)
  
  for (chain_id in seq_along(H_chains)) {
    
    H_chain <- H_chains[[chain_id]]
    
    if (length(dim(H_chain)) != 3L) {
      stop(
        "Chain ",
        chain_id,
        " must be a three-dimensional array."
      )
    }
    
    if (
      dim(H_chain)[1] != R ||
      dim(H_chain)[2] != G
    ) {
      stop(
        "Dimensions of H in chain ",
        chain_id,
        " do not match H_ref."
      )
    }
  }
  
  H_all <- abind::abind(
    H_chains,
    along = 3
  )
  
  probability <- apply(
    H_all,
    MARGIN = c(1, 2),
    FUN = mean
  )
  
  uncertainty <- 4 *
    probability *
    (1 - probability)
  
  factor_names <- rownames(H_ref)
  feature_names <- colnames(H_ref)
  
  if (!is.null(factor_names)) {
    rownames(probability) <- factor_names
    rownames(uncertainty) <- factor_names
  }
  
  if (!is.null(feature_names)) {
    colnames(probability) <- feature_names
    colnames(uncertainty) <- feature_names
  }
  
  return(
    list(
      H_all = H_all,
      probability = probability,
      uncertainty = uncertainty
    )
  )
}
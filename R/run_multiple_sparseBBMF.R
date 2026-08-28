#' Run Multiple Sparse BBMF Chains in Parallel
#'
#' Runs multiple independent sparse Bayesian Boolean Matrix Factorization
#' MCMC chains in parallel using \code{run_sparseBBMF()}.
#'
#' @param X An observed binary matrix of dimension \eqn{K \times G}.
#' @param W_o Initial binary matrix of dimension \eqn{K \times R}.
#' @param H_o Initial binary matrix of dimension \eqn{R \times G}.
#' @param psi_o Initial binary vector of length \eqn{G} containing
#'   feature-specific sparsity indicators.
#' @param n_iter Total number of Gibbs sampling iterations for each chain.
#' @param burn_in Number of initial iterations discarded as burn-in.
#' @param n_chains Number of independent MCMC chains to run.
#' @param a1 First shape parameter of the Beta prior for \code{alpha}.
#' @param a2 Second shape parameter of the Beta prior for \code{alpha}.
#' @param b1 First shape parameter of the Beta distribution used when
#'   \eqn{\psi_g = 1}.
#' @param b2 Second shape parameter of the Beta distribution used when
#'   \eqn{\psi_g = 1}.
#' @param c1 First shape parameter of the Beta distribution used when
#'   \eqn{\psi_g = 0}.
#' @param c2 Second shape parameter of the Beta distribution used when
#'   \eqn{\psi_g = 0}.
#' @param d1 First shape parameter of the Beta prior for \eqn{\pi}.
#' @param d2 Second shape parameter of the Beta prior for \eqn{\pi}.
#'
#' @return A named list of length \code{n_chains}, with each element
#'   containing the output from \code{run_sparseBBMF()}.
#'
#' @details
#' The chains are run in parallel using \code{pbmcapply::pbmclapply()}.
#' The number of cores is restricted to the smaller of \code{n_chains}
#' and the number of available processor cores reported by
#' \code{parallel::detectCores()}.
#'
#' A chain-specific random seed is set using \code{1001 + chain_id}.
#'
#'@export
#'
run_multiple_sparseBBMF <- function(X,W_o,H_o,psi_o,n_iter,
                                    burn_in,n_chains,a1,a2,b1,b2,c1,c2,d1,d2) {
  # Define a helper function to run a single chain with a specific seed
  run_chain <- function(chain_id) {
    set.seed(1001 + chain_id)  # Unique seed for each chain
    
    results <- run_sparseBBMF(
      X,W_o,H_o,psi_o,n_iter,
      burn_in,a1,a2,b1,b2,c1,c2,d1,d2
    )
    
    return(results)
  }
  
  # Use pbmclapply to run chains in parallel and show progress
  all_chain_results <- pbmcapply::pbmclapply(
    X        = 1:n_chains,
    FUN      = run_chain,
    mc.cores = min(n_chains, parallel::detectCores()),
    mc.style = "ETA"
  )
  
  # Rename list elements with chain labels
  names(all_chain_results) <- paste0("Chain_", 1:n_chains)
  
  return(all_chain_results)
}
#' Run Multiple BBMF Chains Sequentially
#'
#' Runs multiple independent Bayesian Boolean Matrix Factorization
#' MCMC chains sequentially using \code{run_BBMF()}.
#'
#' @param X An observed binary matrix of dimension \eqn{K \times G}.
#' @param W_o Initial binary matrix of dimension \eqn{K \times R}
#'   containing the factor activations.
#' @param H_o Initial binary matrix of dimension \eqn{R \times G}
#'   containing the factor-specific feature patterns.
#' @param n_iter Total number of Gibbs sampling iterations for each chain.
#' @param burn_in Number of initial iterations discarded as burn-in for
#'   each chain.
#' @param n_chains Number of independent MCMC chains to run.
#' @param a1 First shape parameter of the Beta prior for the
#'   row-specific probabilities \code{alpha}.
#' @param a2 Second shape parameter of the Beta prior for the
#'   row-specific probabilities \code{alpha}.
#' @param b1 First shape parameter of the Beta prior for the
#'   column-specific probabilities \code{beta}.
#' @param b2 Second shape parameter of the Beta prior for the
#'   column-specific probabilities \code{beta}.
#'
#' @return A named list of length \code{n_chains}. Each element contains
#'   the output from \code{run_BBMF()} for the corresponding MCMC chain.
#'   The list elements are named \code{Chain_1}, \code{Chain_2}, and so on.
#'
#' @details
#' The chains are run sequentially. Before each chain is started,
#' a chain-specific random seed is set using \code{1001 + chain_id}.
#' This gives each chain a different random-number sequence while
#' preserving reproducibility.
#'
#' Each chain starts from the supplied initial matrices \code{W_o}
#' and \code{H_o}.
#'
run_multiple_BBMF_sequential<- function(X,W_o,H_o,n_iter,burn_in,n_chains,a1,a2,b1,b2) {
  # Initialize an empty list to store results from each chain
  all_chain_results <- vector("list", n_chains)
  
  # Run each chain sequentially
  for (chain_id in 1:n_chains) {
    set.seed(1001 + chain_id)  # Unique seed for each chain
    
    # Run the MCMC 
    results <- run_BBMF(X,W_o,H_o,n_iter,burn_in,a1,a2,b1,b2)
    
    # Label and store the results
    chain_label <- paste0("Chain_", chain_id)
    all_chain_results[[chain_label]] <- results
  }
  
  return(all_chain_results)
}
#' Run Multiple BBMF Chains in Parallel
#'
#' Runs multiple independent Bayesian Boolean Matrix Factorization
#' MCMC chains in parallel using \code{run_BBMF()}.
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
#' Each MCMC chain is executed in parallel using \code{pbmclapply()}.
#' Before each chain is started, a chain-specific random seed is set
#' using \code{1000 + chain_id}, giving each chain a different
#' random-number sequence while preserving reproducibility.
#'
#' The number of processor cores is determined using \code{detectCores()}.
#' An ETA-style progress indicator is displayed while the chains are
#' being executed.
#'
#' Each chain starts from the supplied initial matrices \code{W_o}
#' and \code{H_o}.
#'
run_multiple_BBMF_parallel <- function(X,W_o,H_o,n_iter,burn_in, n_chains,a1,a2,b1,b2) {
  # Define a helper function to run a single chain with a specific seed
  run_chain <- function(chain_id) {
    set.seed(1000 + chain_id)  # Unique seed for each chain
    results <- run_BBMF(X,W_o,H_o,n_iter,burn_in,a1,a2,b1,b2)
    return(results)
  }
  
  # Use pbmclapply to run chains in parallel and show progress
  all_chain_results <- pbmcapply::pbmclapply(1:n_chains, run_chain, 
                                             mc.cores = parallel::detectCores(), mc.style = "ETA")
  
  # Rename list elements with chain labels
  names(all_chain_results) <- paste0("Chain_", 1:n_chains)
  
  return(all_chain_results)
}
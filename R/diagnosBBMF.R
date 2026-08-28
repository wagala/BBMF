#' Compute BBMF Posterior Diagnostics
#'
#' Computes iteration-wise diagnostic quantities from posterior samples
#' produced by the Bayesian Boolean Matrix Factorization algorithm.
#' For each MCMC chain, the function evaluates the log-likelihood,
#' log-prior contributions for \code{W} and \code{H}, and the resulting
#' log-posterior value at each retained iteration.
#'
#' @param X An observed binary matrix of dimension \eqn{K \times G}.
#' @param BBMFresult A list containing posterior samples from multiple
#'   BBMF chains. Each chain is expected to contain \code{samples_W},
#'   \code{samples_H}, \code{samples_alpha}, \code{samples_beta},
#'   \code{samples_p11}, and \code{samples_p10}.
#' @param n_iter Total number of stored iterations to evaluate.
#' @param burn_in Number of initial iterations discarded as burn-in.
#' @param num_chains Number of MCMC chains contained in
#'   \code{BBMFresult}.
#'
#' @return A list of length \code{num_chains}. Each element is a data
#'   frame containing:
#' \describe{
#'   \item{\code{log_likhood}}{
#'     Log-likelihood of the observed matrix \code{X}.
#'   }
#'   \item{\code{log_postr}}{
#'     Log-posterior value obtained from the log-likelihood and
#'     prior contributions for \code{W} and \code{H}.
#'   }
#'   \item{\code{priorW}}{
#'     Log-prior probability of the factor-activation matrix \code{W}.
#'   }
#'   \item{\code{priorH}}{
#'     Log-prior probability of the factor-feature matrix \code{H}.
#'   }
#' }
#'
#' @details
#' For each retained posterior sample, the function extracts the sampled
#' matrices \code{W} and \code{H}, the corresponding prior probabilities
#' \code{alpha} and \code{beta}, and the observation probabilities
#' \code{p11} and \code{p10}.
#'
#' The log-posterior quantity is calculated as
#'
#' \deqn{
#' \log p(W,H \mid X)
#' =
#' \log p(X \mid W,H)
#' +
#' \log p(W)
#' +
#' \log p(H),
#' }
#'
#' up to terms not included in the function.
#'
#' @seealso
#' \code{\link{compute_likelihood}},
#' \code{\link{compute_log_priorW}},
#' \code{\link{compute_log_priorH}}
#'
#' @export
#'
diagnosBBMF <- function(X, BBMFresult, n_iter, burn_in, num_chains) {
  
  df_list <- list()
  
  for (chain in seq_len(num_chains)) {
    
    samples_W <- BBMFresult[[chain]]$samples_W
    samples_H <- BBMFresult[[chain]]$samples_H
    samples_alpha <- BBMFresult[[chain]]$samples_alpha
    samples_beta <- BBMFresult[[chain]]$samples_beta
    samples_p11 <- BBMFresult[[chain]]$samples_p11
    samples_p10 <- BBMFresult[[chain]]$samples_p10
    
    iter_len <- n_iter - burn_in
    
    log_likhood <- numeric(iter_len)
    log_postr <- numeric(iter_len)
    priorW <- numeric(iter_len)
    priorH <- numeric(iter_len)
    
    pb <- utils::txtProgressBar(
      min = 0,
      max = iter_len,
      style = 3
    )
    
    for (i in seq_len(iter_len)) {
      
      W <- samples_W[, , i]
      H <- samples_H[, , i]
      alpha <- samples_alpha[, i]
      beta <- samples_beta[, i]
      
      p11 <- samples_p11[i]
      p10 <- samples_p10[i]
      
      log_likhood[i] <- compute_likelihood(
        W,
        H,
        X,
        p11,
        p10
      )
      
      priorW[i] <- compute_log_priorW(
        W,
        alpha
      )
      
      priorH[i] <- compute_log_priorH(
        H,
        beta
      )
      
      log_postr[i] <-
        log_likhood[i] +
        priorW[i] +
        priorH[i]
      
      utils::setTxtProgressBar(
        pb,
        i
      )
    }
    
    close(pb)
    
    results_df <- data.frame(
      log_likhood = log_likhood,
      log_postr = log_postr,
      priorW = priorW,
      priorH = priorH
    )
    
    df_list[[chain]] <- results_df
    
    message(
      "Chain ",
      chain,
      " done!"
    )
  }
  
  return(df_list)
}
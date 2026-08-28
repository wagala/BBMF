#' Run the Bayesian Boolean Matrix Factorization Gibbs Sampler
#'
#' Runs a Gibbs sampler for the Bayesian Boolean Matrix Factorization
#' model using an observed binary matrix \code{X} and initial values
#' for the binary factor matrices \code{W} and \code{H}.
#'
#' @param X An observed binary matrix of dimension \eqn{K \times G}.
#' @param W_o Initial binary matrix of dimension \eqn{K \times R}
#'   containing the factor activations.
#' @param H_o Initial binary matrix of dimension \eqn{R \times G}
#'   containing the factor-specific feature patterns.
#' @param n_iter Total number of Gibbs sampling iterations.
#' @param burn_in Number of initial iterations discarded as burn-in.
#' @param a1 First shape parameter of the Beta prior for the
#'   row-specific probabilities \code{alpha}.
#' @param a2 Second shape parameter of the Beta prior for the
#'   row-specific probabilities \code{alpha}.
#' @param b1 First shape parameter of the Beta prior for the
#'   column-specific probabilities \code{beta}.
#' @param b2 Second shape parameter of the Beta prior for the
#'   column-specific probabilities \code{beta}.
#'
#' @return A list containing posterior samples collected after burn-in:
#' \describe{
#'   \item{\code{samples_p10}}{A numeric vector containing posterior
#'     samples of \code{p10}.}
#'   \item{\code{samples_p11}}{A numeric vector containing posterior
#'     samples of \code{p11}.}
#'   \item{\code{samples_alpha}}{A matrix containing posterior samples
#'     of the row-specific probabilities \code{alpha}.}
#'   \item{\code{samples_beta}}{A matrix containing posterior samples
#'     of the column-specific probabilities \code{beta}.}
#'   \item{\code{samples_W}}{A three-dimensional array containing
#'     posterior samples of \code{W}.}
#'   \item{\code{samples_H}}{A three-dimensional array containing
#'     posterior samples of \code{H}.}
#' }
#'
#' @details
#' At each Gibbs sampling iteration, the function updates \code{alpha}
#' and \code{beta}, samples the observation probabilities \code{p11}
#' and \code{p10}, and subsequently updates the binary factor matrices
#' \code{W} and \code{H}.
#'
#' Samples are stored only after the burn-in period. A progress bar is
#' displayed while the Gibbs sampler is running.
#'
run_BBMF<-function(X,W_o,H_o,n_iter,burn_in,a1,a2,b1,b2){
  #Gibbs Sampling
  W<-W_o
  H<-H_o
  R<-ncol(W_o)
  K<-nrow(X)
  G<-ncol(X)
  
  #Arrays for storing the samples
  samples_W <-array(0, dim = c(K, R, n_iter - burn_in))
  samples_H <- array(0, dim = c(R, G, n_iter - burn_in))
  samples_alpha <- matrix(0, nrow = K, ncol = n_iter-burn_in)
  samples_beta <- matrix(0, nrow = G, ncol = n_iter - burn_in)
  samples_p11 <- numeric(n_iter - burn_in)
  samples_p10 <- numeric(n_iter - burn_in)
  
  # Initialize the progress bar
  pb_chain <- progress::progress_bar$new(format = "working [:bar] :percent in :elapsed",
                               total = n_iter, clear = FALSE, width = 100)
  
  for (iter in 1:n_iter){
    
    # Update alpha and beta
    alpha<-update_alpha(W,a1,a2)
    beta<-update_beta(H,b1,b2)
    
    # Compute p11 and p10
    pvals<-update_p(W,H,X)
    #pvals<-update_p11_p10(W,H,X)
    p11 <- pvals[1]
    p10 <- pvals[2]
    
    # Update W and H
    W<-update_W(W,H,X,alpha,p11=p11,p10=p10)
    H<-update_H(W,H,X,beta,p11=p11,p10=p10)
    # W<-update_W_Block_rows(W,H,X,alpha,p11=p11,p10=p10)
    # H<-update_H_Block_columns(W,H,X,beta,p11=p11,p10=p10)
    
    # Update the progress bar
    pb_chain$tick()
    
    #---------------------------------------------------------------------------
    #Store samples after burn-in
    if (iter > burn_in) {
      samples_W[,,iter-burn_in]<-W
      samples_H[,,iter-burn_in]<-H
      samples_alpha[,iter-burn_in]<-alpha
      samples_beta[,iter-burn_in] <-beta
      samples_p11[iter - burn_in] <- p11
      samples_p10[iter - burn_in] <- p10
    }
  }
  
  #Terminate and clear the progress bar
  pb_chain$terminate()
  
  #Return a List of Arrays
  return(list(samples_p10=samples_p10,samples_p11=samples_p11,samples_alpha=samples_alpha,
              samples_beta=samples_beta,samples_W=samples_W,samples_H=samples_H))  
}
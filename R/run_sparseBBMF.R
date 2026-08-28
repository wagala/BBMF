#' Run the Sparse Bayesian Boolean Matrix Factorization Gibbs Sampler
#'
#' Runs a Gibbs sampler for the sparse Bayesian Boolean Matrix
#' Factorization model using an observed binary matrix \code{X},
#' initial factor matrices \code{W_o} and \code{H_o}, and initial
#' sparsity indicators \code{psi_o}.
#'
#' @param X An observed binary matrix of dimension \eqn{K \times G}.
#' @param W_o Initial binary matrix of dimension \eqn{K \times R}
#'   containing the factor activations.
#' @param H_o Initial binary matrix of dimension \eqn{R \times G}
#'   containing the factor-specific feature patterns.
#' @param psi_o Initial binary vector of length \eqn{G} containing the
#'   feature-specific sparsity indicators.
#' @param n_iter Total number of Gibbs sampling iterations.
#' @param burn_in Number of initial iterations discarded as burn-in.
#' @param a1 First shape parameter of the Beta prior for the
#'   row-specific probabilities \code{alpha}.
#' @param a2 Second shape parameter of the Beta prior for the
#'   row-specific probabilities \code{alpha}.
#' @param b1 First shape parameter of the Beta distribution for
#'   \eqn{\beta_g} when \eqn{\psi_g = 1}.
#' @param b2 Second shape parameter of the Beta distribution for
#'   \eqn{\beta_g} when \eqn{\psi_g = 1}.
#' @param c1 First shape parameter of the Beta distribution for
#'   \eqn{\beta_g} when \eqn{\psi_g = 0}.
#' @param c2 Second shape parameter of the Beta distribution for
#'   \eqn{\beta_g} when \eqn{\psi_g = 0}.
#' @param d1 First shape parameter of the Beta prior for the global
#'   sparsity probability \eqn{\pi}.
#' @param d2 Second shape parameter of the Beta prior for the global
#'   sparsity probability \eqn{\pi}.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{samples_p10}}{Posterior samples of \code{p10}.}
#'   \item{\code{samples_p11}}{Posterior samples of \code{p11}.}
#'   \item{\code{samples_alpha}}{Posterior samples of the row-specific
#'     probabilities \code{alpha}.}
#'   \item{\code{samples_beta}}{Posterior samples of the feature-specific
#'     probabilities \eqn{\beta_g}.}
#'   \item{\code{samples_pi}}{Posterior samples of the global sparsity
#'     probability \eqn{\pi}.}
#'   \item{\code{samples_psi}}{Posterior samples of the feature-specific
#'     sparsity indicators \eqn{\psi_g}.}
#'   \item{\code{samples_W}}{Posterior samples of the binary factor
#'     matrix \code{W}.}
#'   \item{\code{samples_H}}{Posterior samples of the binary factor
#'     matrix \code{H}.}
#' }
#'
#'@export
#'
run_sparseBBMF<- function(X, W_o, H_o, psi_o, n_iter, burn_in,
                          a1, a2, b1, b2, c1, c2, d1, d2) {
  # Initialize W and H with the initial values provided
  W <- W_o
  H <- H_o
  
  # Dimensions for W and H
  K <- nrow(W_o)
  R <- ncol(W_o)
  G <- ncol(X)
  psi_g<-psi_o
  
  # Pre-allocate space for samples, including space for initial values
  num_samples <- n_iter - burn_in + 1
  samples_W <- array(0, dim = c(K, R, num_samples))
  samples_H <- array(0, dim = c(R, G, num_samples))
  samples_alpha <- matrix(0, nrow = K, ncol = num_samples)
  samples_beta <- matrix(0, nrow = G, ncol = num_samples)
  samples_pi <- numeric(num_samples)
  samples_psi <- matrix(0, nrow = G, ncol = num_samples)
  samples_p11 <- numeric(num_samples)
  samples_p10 <- numeric(num_samples)
  
  # Store initial values
  samples_W[,,1] <- W
  samples_H[,,1] <- H
  samples_alpha[,1] <- update_alpha(W, a1, a2)
  samples_beta[,1] <- update_beta_g(H, psi_o, b1, b2, c1, c2)
  samples_pi[1] <- update_pi(psi_o, d1, d2)
  samples_psi[,1] <- psi_o
  pvals <- update_p(W, H, X)
  samples_p11[1] <- pvals[1]
  samples_p10[1] <- pvals[2]
  
  # Initialize the progress bar
  pb_chain <- progress::progress_bar$new(
    format = "working [:bar] :percent in :elapsed",
    total = n_iter,
    clear = FALSE,
    width = 100
  )
  
  for (iter in 1:n_iter) {
    # Perform updates for W, H, and other parameters
    alpha <- update_alpha(W, a1, a2)
    pi_s <- update_pi(psi_g, d1, d2)
    beta_g <- update_beta_g(H, psi_g, b1, b2, c1, c2)
    psi_g <- update_psi(beta_g, b1, b2, c1, c2, pi_s)
    pvals <- update_p(W, H, X)
    p11 <- pvals[1]
    p10 <- pvals[2]
    
    # Update W and H
    # W<-update_W(W,H,X,alpha,p11=p11,p10=p10)
    # H<-update_H(W,H,X,beta= beta_g,p11=p11,p10=p10)
    W<-update_W_Block_rows(W,H,X,alpha,p11=p11,p10=p10)
    H<-update_H_Block_columns(W,H,X,beta=beta_g,p11=p11,p10=p10)
    
    # Update the progress bar
    pb_chain$tick()
    
    # Store the results
    if (iter > burn_in) {
      idx <- iter - burn_in + 1
      samples_W[,,idx] <- W
      samples_H[,,idx] <- H
      samples_alpha[,idx] <- alpha
      samples_beta[,idx] <- beta_g
      samples_pi[idx] <- pi_s
      samples_psi[,idx] <- psi_g
      samples_p11[idx] <- p11
      samples_p10[idx] <- p10
    }
  }
  
  # Terminate the progress bar
  pb_chain$terminate()
  
  # Return a List of Arrays
  return(list(
    samples_p10 = samples_p10,
    samples_p11 = samples_p11,
    samples_alpha = samples_alpha,
    samples_beta = samples_beta,
    samples_pi = samples_pi,
    samples_psi = samples_psi,
    samples_W = samples_W,
    samples_H = samples_H
  ))
}
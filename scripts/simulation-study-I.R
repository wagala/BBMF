## =============================================================================
## Simulation Study I
## Adolphus Wagala
## August 26, 2026
## =============================================================================


## =============================================================================
## Load BBMF functions
## =============================================================================

r_files <- list.files(
  path = "R",
  pattern = "\\.R$",
  full.names = TRUE
)

invisible(lapply(r_files, source))


## =============================================================================
## Simulation setup
## =============================================================================

set.seed(101)

# Number of samples
K <- 70

# Number of latent factors
R <- 4

# Chromosome arms:
# 1p, 1q, ..., 22p, 22q
arms <- unlist(
  lapply(
    1:22,
    function(i) paste0(i, c("p", "q"))
  )
)

# Number of genomic features
G <- length(arms)

# Sample labels
sample_ids <- paste0(
  "X",
  seq_len(K)
)


## =============================================================================
## Simulate H
## =============================================================================

# High-density probabilities for the four factors
theta_high <- c(
  0.80,
  0.45,
  0.65,
  0.05
)

# Low-density probabilities for the four factors
theta_low <- c(
  0.05,
  0.05,
  0.25,
  0.05
)

# Genomic regions associated with each factor
regions <- list(
  1:22,     # Factor 1
  20:40,    # Factor 2
  25:44,    # Factor 3
  1:44      # Factor 4
)

# Initialize H
# Dimensions: R x G
H_true <- matrix(
  0,
  nrow = R,
  ncol = G
)

# Generate factor-feature patterns
for (r in seq_len(R)) {
  
  for (g in seq_len(G)) {
    
    if (g %in% regions[[r]]) {
      
      H_true[r, g] <- stats::rbinom(
        1,
        size = 1,
        prob = theta_high[r]
      )
      
    } else {
      
      H_true[r, g] <- stats::rbinom(
        1,
        size = 1,
        prob = theta_low[r]
      )
    }
  }
}


## =============================================================================
## Simulate factor prevalences
## =============================================================================

# Beta prior hyperparameters for factor prevalence
a_W <- 2
b_W <- 6

# Factor-specific prevalences
# Denoted by nu_r in the manuscript
nu <- stats::rbeta(
  R,
  shape1 = a_W,
  shape2 = b_W
)

# Display factor prevalences
nu


## =============================================================================
## Simulate W
## =============================================================================

# Initialize W
# Dimensions: K x R
W_true <- matrix(
  0,
  nrow = K,
  ncol = R
)

# Generate sample-factor activations
for (r in seq_len(R)) {
  
  W_true[, r] <- stats::rbinom(
    K,
    size = 1,
    prob = nu[r]
  )
}


## =============================================================================
## Generate the noise-free Boolean matrix
## =============================================================================

X_true <- bool_prod(
  W_true,
  H_true
)


## =============================================================================
## Add matrix labels
## =============================================================================

rownames(X_true) <- sample_ids
colnames(X_true) <- arms

rownames(W_true) <- sample_ids
colnames(W_true) <- paste0(
  "Factor_",
  seq_len(R)
)

rownames(H_true) <- paste0(
  "Factor_",
  seq_len(R)
)

colnames(H_true) <- arms


## =============================================================================
## Add noise
## =============================================================================

# Randomly flip 20% of entries
X_noisy <- noise_bool(
  X_true,
  noise_level = 0.20,
  seed = 123
)

rownames(X_noisy) <- sample_ids
colnames(X_noisy) <- arms


## =============================================================================
## Prepare simulated matrices for visualization
## =============================================================================

W_true_long <- reshape2::melt(
  W_true
)

W_true_long$Var2 <- as.factor(
  W_true_long$Var2
)


H_true_long <- reshape2::melt(
  H_true
)

H_true_long$Var1 <- as.factor(
  H_true_long$Var1
)


X_true_long <- reshape2::melt(
  X_true
)

X_true_long$Var1 <- as.factor(
  X_true_long$Var1
)


X_noisy_long <- reshape2::melt(
  X_noisy
)

X_noisy_long$Var1 <- as.factor(
  X_noisy_long$Var1
)


## =============================================================================
## Visualize the true W matrix
## =============================================================================

pW <- matrixPlot(
  long_data = W_true_long,
  x_value = Var2,
  y_value = Var1,
  value = value,
  title = "True W"
)

pW


## =============================================================================
## Visualize the true H matrix
## =============================================================================

pH <- matrixPlot(
  long_data = H_true_long,
  x_value = Var2,
  y_value = Var1,
  value = value,
  title = "True H"
)

pH <- pH +
  ggplot2::theme(
    text = ggplot2::element_text(
      size = 7
    ),
    axis.text.x = ggplot2::element_text(
      size = 10,
      angle = 90,
      vjust = 0.5
    ),
    axis.text.y = ggplot2::element_text(
      size = 10
    )
  )

pH


## =============================================================================
## Visualize the true X matrix
## =============================================================================

pX <- matrixPlot(
  long_data = X_true_long,
  x_value = Var2,
  y_value = Var1,
  value = value,
  title = "True X"
)

pX <- pX +
  ggplot2::theme(
    text = ggplot2::element_text(
      size = 7
    ),
    axis.text.x = ggplot2::element_text(
      size = 8,
      angle = 90,
      vjust = 0.5
    ),
    axis.text.y = ggplot2::element_text(
      size = 8
    )
  )

pX


## =============================================================================
## Visualize the noisy X matrix
## =============================================================================

pX_noisy <- matrixPlot(
  long_data = X_noisy_long,
  x_value = Var2,
  y_value = Var1,
  value = value,
  title = "Observed X with 20% Noise"
)

pX_noisy <- pX_noisy +
  ggplot2::theme(
    text = ggplot2::element_text(
      size = 7
    ),
    axis.text.x = ggplot2::element_text(
      size = 8,
      angle = 90,
      vjust = 0.5
    ),
    axis.text.y = ggplot2::element_text(
      size = 8
    )
  )

pX_noisy


## =============================================================================
## Obtain initial values using Asso
## =============================================================================

asso_fit <- Asso_BBMF(
  X = X_noisy,
  R = R,
  seed = 101
)

# Extract Asso estimates
W_asso <- asso_fit$W_asso
H_asso <- asso_fit$H_asso
X_asso <- asso_fit$X_asso


## =============================================================================
## Add labels to Asso matrices
## =============================================================================

rownames(W_asso) <- sample_ids

colnames(W_asso) <- paste0(
  "Factor_",
  seq_len(R)
)

rownames(H_asso) <- paste0(
  "Factor_",
  seq_len(R)
)

colnames(H_asso) <- arms

rownames(X_asso) <- sample_ids
colnames(X_asso) <- arms


## =============================================================================
## Prepare Asso matrices for visualization
## =============================================================================

W_asso_long <- reshape2::melt(
  W_asso
)

W_asso_long$Var2 <- as.factor(
  W_asso_long$Var2
)


H_asso_long <- reshape2::melt(
  H_asso
)

H_asso_long$Var1 <- as.factor(
  H_asso_long$Var1
)


X_asso_long <- reshape2::melt(
  X_asso
)


## =============================================================================
## Visualize Asso W
## =============================================================================

pW_asso <- matrixPlot(
  long_data = W_asso_long,
  x_value = Var2,
  y_value = Var1,
  value = value,
  title = "Asso Estimate of W"
)

pW_asso


## =============================================================================
## Visualize Asso H
## =============================================================================

pH_asso <- matrixPlot(
  long_data = H_asso_long,
  x_value = Var2,
  y_value = Var1,
  value = value,
  title = "Asso Estimate of H"
)

pH_asso


## =============================================================================
## Visualize Asso reconstruction
## =============================================================================

pX_asso <- matrixPlot(
  long_data = X_asso_long,
  x_value = Var2,
  y_value = Var1,
  value = value,
  title = "Asso Reconstruction"
)

pX_asso


## =============================================================================
## Sparse BBMF settings
## =============================================================================

XX1 <- X_noisy

n_iter <- 1000
burn_in <- 0
n_chains <- 4

G <- ncol(XX1)


## =============================================================================
## Initialize sparsity indicators
## =============================================================================

set.seed(2026)

psi_o <- stats::rbinom(
  G,
  size = 1,
  prob = 0.5
)

psi_o


## =============================================================================
## Fit one sparse BBMF chain
## =============================================================================

# This is useful for checking that the sampler runs correctly
# before starting multiple parallel chains.

# test_fit <- run_sparseBBMF(
#   X = XX1,
#   W_o = W_asso,
#   H_o = H_asso,
#   psi_o = psi_o,
#   n_iter = n_iter,
#   burn_in = burn_in,
#   a1 = 1,
#   a2 = 1,
#   b1 = 1,
#   b2 = 1,
#   c1 = 1,
#   c2 = 1,
#   d1 = 1,
#   d2 = 1
# )


## =============================================================================
## Run multiple sparse BBMF chains
## =============================================================================

fit_BBMF <- run_multiple_sparseBBMF(
  X = XX1,
  W_o = W_asso,
  H_o = H_asso,
  psi_o = psi_o,
  n_iter = n_iter,
  burn_in = burn_in,
  n_chains = n_chains,
  a1 = 1,
  a2 = 1,
  b1 = 1,
  b2 = 1,
  c1 = 1,
  c2 = 1,
  d1 = 1,
  d2 = 1
)


## =============================================================================
## Inspect fitted object
## =============================================================================

names(fit_BBMF)

names(fit_BBMF$Chain_1)


## =============================================================================
## Example: dimensions of posterior samples from Chain 1
## =============================================================================

dim(
  fit_BBMF$Chain_1$samples_W
)

dim(
  fit_BBMF$Chain_1$samples_H
)

dim(
  fit_BBMF$Chain_1$samples_alpha
)

dim(
  fit_BBMF$Chain_1$samples_beta
)

length(
  fit_BBMF$Chain_1$samples_pi
)

dim(
  fit_BBMF$Chain_1$samples_psi
)

length(
  fit_BBMF$Chain_1$samples_p11
)

length(
  fit_BBMF$Chain_1$samples_p10
)
## =============================================================================
## BBMF diagnostics
## =============================================================================

# run_sparseBBMF() stores the initial state plus n_iter iterations.
# With n_iter = 1000 and burn_in = 0, each chain has 1001 stored samples.
nn <- n_iter + 1

diagnosticBBMF1 <- diagnosBBMF(
  X = XX1,
  BBMFresult = fit_BBMF,
  n_iter = nn,
  burn_in = 0,
  num_chains = n_chains
)

# Group each diagnostic quantity across the four chains
grouped_data1 <- group_diagnost_List(
  diagnosticBBMF1
)


## =============================================================================
## Helper function for diagnostic trace plots
## =============================================================================

plot_diagnostic <- function(data,
                            prefix,
                            plot_title,
                            y_label,
                            n_chains = 4) {
  
  # Remove the initial stored state.
  # Remaining rows correspond to MCMC iterations 1, ..., n_iter.
  data <- data[-1, , drop = FALSE]
  
  # Add iteration number
  data$iteration <- seq_len(nrow(data))
  
  # Convert to long format
  data_long <- tidyr::pivot_longer(
    data,
    cols = tidyselect::starts_with(prefix),
    names_to = "diagnostic",
    values_to = "Value"
  )
  
  # Chain labels
  data_long$chain <- factor(
    data_long$diagnostic,
    levels = paste0(prefix, "_", seq_len(n_chains)),
    labels = paste0("chain_", seq_len(n_chains))
  )
  
  # Display "chain 1" instead of "chain_1"
  chain_lab <- function(x) {
    gsub("^chain_", "chain ", x)
  }
  
  ggplot2::ggplot(
    data_long,
    ggplot2::aes(
      x = iteration,
      y = Value,
      color = chain,
      group = chain
    )
  ) +
    ggplot2::geom_line(
      linewidth = 0.25
    ) +
    ggplot2::facet_wrap(
      ~ chain,
      scales = "fixed",
      nrow = 2,
      ncol = 2,
      labeller = ggplot2::labeller(
        chain = chain_lab
      )
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = 12
      ),
      legend.title = ggplot2::element_blank(),
      legend.position = "none",
      text = ggplot2::element_text(
        size = 12
      ),
      axis.text.y = ggplot2::element_text(
        size = 10
      ),
      axis.text.x = ggplot2::element_text(
        size = 10,
        angle = 90,
        vjust = 0.5,
        hjust = 0.5
      ),
      strip.background = ggplot2::element_rect(
        fill = "white",
        colour = "grey70"
      ),
      strip.text = ggplot2::element_text(
        size = 12,
        colour = "black"
      )
    ) +
    ggplot2::labs(
      title = plot_title,
      x = "Iteration",
      y = y_label
    )
}


## =============================================================================
## Log-likelihood trace
## =============================================================================

plot1 <- plot_diagnostic(
  data = grouped_data1$log_likhood,
  prefix = "log_likhood",
  plot_title = "Log Likelihood",
  y_label = "Log likelihood",
  n_chains = n_chains
)

plot1


## =============================================================================
## Log-posterior trace
## =============================================================================

plot2 <- plot_diagnostic(
  data = grouped_data1$log_postr,
  prefix = "log_postr",
  plot_title = "Log Posterior",
  y_label = "Log posterior",
  n_chains = n_chains
)

plot2


## =============================================================================
## Log-prior for W
## =============================================================================

plot3 <- plot_diagnostic(
  data = grouped_data1$priorW,
  prefix = "priorW",
  plot_title = "Log Prior for W",
  y_label = "Log prior for W",
  n_chains = n_chains
)

plot3


## =============================================================================
## Log-prior for H
## =============================================================================

plot4 <- plot_diagnostic(
  data = grouped_data1$priorH,
  prefix = "priorH",
  plot_title = "Log Prior for H",
  y_label = "Log prior for H",
  n_chains = n_chains
)

plot4


## =============================================================================
## Combine diagnostic plots
## =============================================================================

diagnMM <- gridExtra::grid.arrange(
  plot1,
  plot2,
  plot3,
  plot4,
  nrow = 2,
  top = ""
)
## =============================================================================
## Select the posterior sample with the largest log-posterior value
## =============================================================================

# Extract log-posterior values for all chains
logPost <- grouped_data1$log_postr

# Remove the initial stored state
logPost_mcmc <- logPost[-1, , drop = FALSE]

# Find the overall maximum log-posterior value across all chains
max_logPost <- max(
  as.matrix(logPost_mcmc),
  na.rm = TRUE
)

max_logPost


## =============================================================================
## Locate the chain and iteration containing the maximum
## =============================================================================

best_position <- which(
  as.matrix(logPost_mcmc) == max_logPost,
  arr.ind = TRUE
)[1, ]

# Row = MCMC iteration
best_iteration <- best_position[1]

# Column = chain
best_chain <- best_position[2]

best_chain
best_iteration


## =============================================================================
## Convert diagnostic iteration to stored-sample index
## =============================================================================

# samples_W and samples_H contain the initial state in position 1.
# Therefore MCMC iteration j corresponds to stored sample j + 1.
best_sample_index <- best_iteration + 1

best_sample_index


## =============================================================================
## Extract the best posterior sample
## =============================================================================

best_fit <- fit_BBMF[[best_chain]]

W_max <- best_fit$samples_W[
  ,
  ,
  best_sample_index
]

H_max <- best_fit$samples_H[
  ,
  ,
  best_sample_index
]


## =============================================================================
## Boolean reconstruction
## =============================================================================

X_tilde_BBMF <- bool_prod(
  W_max,
  H_max
)


## =============================================================================
## Summary of selected posterior sample
## =============================================================================

cat(
  "Best chain:", best_chain, "\n",
  "MCMC iteration:", best_iteration, "\n",
  "Stored sample index:", best_sample_index, "\n",
  "Maximum log-posterior:", max_logPost, "\n"
)
## =============================================================================
## Add matrix names
## =============================================================================

factor_names <- paste0(
  "Factor_",
  seq_len(ncol(W_max))
)

# W: samples x factors
rownames(W_max) <- rownames(X_true)
colnames(W_max) <- factor_names

# H: factors x genomic features
rownames(H_max) <- factor_names
colnames(H_max) <- colnames(X_true)

# Boolean reconstruction
X_tilde_BBMF <- bool_prod(
  W_max,
  H_max
)

# X_tilde: samples x genomic features
rownames(X_tilde_BBMF) <- rownames(X_true)
colnames(X_tilde_BBMF) <- colnames(X_true)

## =============================================================================
## Factor similarity: Truth vs Asso
## =============================================================================

sim_asso <- factor_similarity(
  H_ref = H_true,
  H_est = H_asso
)

sim_asso$alignment$matching

sim_asso$pair_similarities

sim_asso$mean_similarity

sim_asso$sum_similarity


Asso_FactorsH <- plot_factor_similarity(
  similarity_result = sim_asso,
  x_label = "Asso Estimated Factors",
  y_label = "True Latent Factors",
  title = "True H vs Asso Estimated H"
)

Asso_FactorsH
## =============================================================================
## Factor similarity: Truth vs sparse BBMF
## =============================================================================

sim_BBMF <- factor_similarity(
  H_ref = H_true,
  H_est = H_max
)

sim_BBMF$alignment$matching

sim_BBMF$pair_similarities

sim_BBMF$mean_similarity

sim_BBMF$sum_similarity

BBMF_FactorsH <- plot_factor_similarity(
  similarity_result = sim_BBMF,
  x_label = "BBMF Estimated Factors",
  y_label = "True Latent Factors",
  title = "True H vs BBMF Estimated H"
)

BBMF_FactorsH

gridExtra::grid.arrange(
  Asso_FactorsH,
  BBMF_FactorsH,
  ncol = 2
)

## =============================================================================
## Align estimated factors to the true factors
## =============================================================================

## -----------------------------------------------------------------------------
## Asso alignment
## -----------------------------------------------------------------------------

sim_asso <- factor_similarity(
  H_ref = H_true,
  H_est = H_asso
)

perm_asso <- sim_asso$permutation

# Align rows of H
H_asso_aligned <- H_asso[
  perm_asso,
  ,
  drop = FALSE
]

# Apply the same permutation to columns of W
W_asso_aligned <- W_asso[
  ,
  perm_asso,
  drop = FALSE
]


## -----------------------------------------------------------------------------
## Sparse BBMF alignment
## -----------------------------------------------------------------------------

sim_BBMF <- factor_similarity(
  H_ref = H_true,
  H_est = H_max
)

perm_BBMF <- sim_BBMF$permutation

# Align rows of H
H_max_aligned <- H_max[
  perm_BBMF,
  ,
  drop = FALSE
]

# Apply the same permutation to columns of W
W_max_aligned <- W_max[
  ,
  perm_BBMF,
  drop = FALSE
]


## =============================================================================
## Add consistent factor, sample, and feature names
## =============================================================================

factor_names <- paste0(
  "Factor_",
  seq_len(R)
)


## True matrices
rownames(W_true) <- rownames(X_true)
colnames(W_true) <- factor_names

rownames(H_true) <- factor_names
colnames(H_true) <- colnames(X_true)


## Asso matrices
rownames(W_asso_aligned) <- rownames(X_true)
colnames(W_asso_aligned) <- factor_names

rownames(H_asso_aligned) <- factor_names
colnames(H_asso_aligned) <- colnames(X_true)


## BBMF matrices
rownames(W_max_aligned) <- rownames(X_true)
colnames(W_max_aligned) <- factor_names

rownames(H_max_aligned) <- factor_names
colnames(H_max_aligned) <- colnames(X_true)


## =============================================================================
## Reconstruct X after alignment
## =============================================================================

X_asso_aligned <- bool_prod(
  W_asso_aligned,
  H_asso_aligned
)

rownames(X_asso_aligned) <- rownames(X_true)
colnames(X_asso_aligned) <- colnames(X_true)


X_tilde_BBMF_aligned <- bool_prod(
  W_max_aligned,
  H_max_aligned
)

rownames(X_tilde_BBMF_aligned) <- rownames(X_true)
colnames(X_tilde_BBMF_aligned) <- colnames(X_true)


## =============================================================================
## Check that alignment does not change the reconstruction
## =============================================================================

identical(
  X_asso,
  X_asso_aligned
)

identical(
  X_tilde_BBMF,
  X_tilde_BBMF_aligned
)


## =============================================================================
## Convert aligned matrices to long format
## =============================================================================

## True matrices
W_true_df <- reshape2::melt(W_true)
H_true_df <- reshape2::melt(H_true)

## Asso matrices
W_asso_aligned_df <- reshape2::melt(
  W_asso_aligned
)

H_asso_aligned_df <- reshape2::melt(
  H_asso_aligned
)

X_asso_aligned_df <- reshape2::melt(
  X_asso_aligned
)

## BBMF matrices
W_max_aligned_df <- reshape2::melt(
  W_max_aligned
)

H_max_aligned_df <- reshape2::melt(
  H_max_aligned
)

X_tilde_BBMF_aligned_df <- reshape2::melt(
  X_tilde_BBMF_aligned
)


## -----------------------------------------------------------------------------
## Preserve factor ordering
## -----------------------------------------------------------------------------

W_true_df$Var2 <- factor(
  W_true_df$Var2,
  levels = factor_names
)

H_true_df$Var1 <- factor(
  H_true_df$Var1,
  levels = factor_names
)


W_asso_aligned_df$Var2 <- factor(
  W_asso_aligned_df$Var2,
  levels = factor_names
)

H_asso_aligned_df$Var1 <- factor(
  H_asso_aligned_df$Var1,
  levels = factor_names
)


W_max_aligned_df$Var2 <- factor(
  W_max_aligned_df$Var2,
  levels = factor_names
)

H_max_aligned_df$Var1 <- factor(
  H_max_aligned_df$Var1,
  levels = factor_names
)


## =============================================================================
## Plot W matrices
## =============================================================================

W_true_graph <- matrixPlot(
  long_data = W_true_df,
  x_value = Var2,
  y_value = Var1,
  value = value,
  title = "True W"
)


W_asso_graph <- matrixPlot(
  long_data = W_asso_aligned_df,
  x_value = Var2,
  y_value = Var1,
  value = value,
  title = "Asso Estimate of W"
)


W_BBMF_graph <- matrixPlot(
  long_data = W_max_aligned_df,
  x_value = Var2,
  y_value = Var1,
  value = value,
  title = "Sparse BBMF Estimate of W"
)


## =============================================================================
## Plot H matrices
## =============================================================================

H_true_graph <- matrixPlot(
  long_data = H_true_df,
  x_value = Var2,
  y_value = Var1,
  value = value,
  title = "True H"
) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      size = 10,
      angle = 90,
      vjust = 0.5
    ),
    axis.text.y = ggplot2::element_text(
      size = 10
    )
  )


H_asso_graph <- matrixPlot(
  long_data = H_asso_aligned_df,
  x_value = Var2,
  y_value = Var1,
  value = value,
  title = "Asso Estimate of H"
) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      size = 10,
      angle = 90,
      vjust = 0.5
    ),
    axis.text.y = ggplot2::element_text(
      size = 10
    )
  )


H_BBMF_graph <- matrixPlot(
  long_data = H_max_aligned_df,
  x_value = Var2,
  y_value = Var1,
  value = value,
  title = "Sparse BBMF Estimate of H"
) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      size = 10,
      angle = 90,
      vjust = 0.5
    ),
    axis.text.y = ggplot2::element_text(
      size = 10
    )
  )


## =============================================================================
## Plot X reconstructions
## =============================================================================

X_asso_graph <- matrixPlot(
  long_data = X_asso_aligned_df,
  x_value = Var2,
  y_value = Var1,
  value = value,
  title = "Asso Reconstruction"
) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      size = 8,
      angle = 90,
      vjust = 0.5
    ),
    axis.text.y = ggplot2::element_text(
      size = 8
    )
  )


X_BBMF_graph <- matrixPlot(
  long_data = X_tilde_BBMF_aligned_df,
  x_value = Var2,
  y_value = Var1,
  value = value,
  title = "Sparse BBMF Reconstruction"
) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      size = 8,
      angle = 90,
      vjust = 0.5
    ),
    axis.text.y = ggplot2::element_text(
      size = 8
    )
  )


## =============================================================================
## Compare X matrices
## =============================================================================

gridExtra::grid.arrange(
  pX,
  pX_noisy,
  X_asso_graph,
  X_BBMF_graph,
  ncol = 2
)


## =============================================================================
## Compare aligned W matrices
## =============================================================================

gridExtra::grid.arrange(
  W_true_graph,
  W_asso_graph,
  W_BBMF_graph,
  ncol = 3
)


## =============================================================================
## Compare aligned H matrices
## =============================================================================

gridExtra::grid.arrange(
  H_true_graph,
  H_asso_graph,
  H_BBMF_graph,
  ncol = 1
)

## =============================================================================
## Aggregate Factor Alignment Similarity (AFAS)
## Simulation Study I
## =============================================================================

# AFAS measures recovery of the latent factor structure in H.
# For each posterior sample:
#   1. Compute pairwise Jaccard similarities between rows of H_true
#      and rows of the sampled H.
#   2. Align factors using the Hungarian algorithm.
#   3. Average the Jaccard similarities of the optimally matched factors.
#
# AFAS ranges from 0 to 1, with larger values indicating better
# recovery of the true latent factor structure.


## =============================================================================
## Compute AFAS for every posterior sample in every chain
## =============================================================================

afas_by_chain <- lapply(
  fit_BBMF,
  function(chain_fit) {
    
    get_mean_sim_per_chain(
      H_true = H_true,
      ff = chain_fit$samples_H
    )
  }
)


## -----------------------------------------------------------------------------
## Combine chains into one matrix
## -----------------------------------------------------------------------------

afas_matrix <- do.call(
  cbind,
  afas_by_chain
)

# Use informative chain names
colnames(afas_matrix) <- paste0(
  "Chain_",
  seq_len(ncol(afas_matrix))
)

dim(afas_matrix)


## =============================================================================
## Remove the initial stored state
## =============================================================================

# run_sparseBBMF() stores the initial H in position 1.
# Remove it so that rows correspond to MCMC iterations 1, ..., n_iter.

afas_matrix <- afas_matrix[
  -1,
  ,
  drop = FALSE
]

# Confirm dimensions
dim(afas_matrix)


## =============================================================================
## Prepare AFAS data for visualization
## =============================================================================

afas_df <- data.frame(
  iteration = seq_len(nrow(afas_matrix)),
  afas_matrix,
  check.names = FALSE
)

afas_long <- tidyr::pivot_longer(
  afas_df,
  cols = tidyselect::starts_with("Chain_"),
  names_to = "chain",
  values_to = "AFAS"
)

# Display chain labels as "Chain 1", ..., "Chain 4"
afas_long$chain <- factor(
  afas_long$chain,
  levels = paste0(
    "Chain_",
    seq_len(ncol(afas_matrix))
  ),
  labels = paste(
    "Chain",
    seq_len(ncol(afas_matrix))
  )
)


## =============================================================================
## AFAS summary by chain
## =============================================================================

afas_summary <- do.call(
  rbind,
  lapply(
    seq_len(ncol(afas_matrix)),
    function(j) {
      
      x <- afas_matrix[, j]
      
      data.frame(
        Chain = paste("Chain", j),
        Mean = mean(x, na.rm = TRUE),
        SD = stats::sd(x, na.rm = TRUE),
        Median = stats::median(x, na.rm = TRUE),
        Q025 = stats::quantile(
          x,
          probs = 0.025,
          na.rm = TRUE,
          names = FALSE
        ),
        Q975 = stats::quantile(
          x,
          probs = 0.975,
          na.rm = TRUE,
          names = FALSE
        )
      )
    }
  )
)

afas_summary


## =============================================================================
## AFAS trace plots
## =============================================================================

AFAS_trace <- ggplot2::ggplot(
  afas_long,
  ggplot2::aes(
    x = iteration,
    y = AFAS
  )
) +
  ggplot2::geom_line(
    linewidth = 0.3,
    alpha = 0.8
  ) +
  ggplot2::facet_wrap(
    ~ chain,
    nrow = 2,
    ncol = 2,
    scales = "fixed"
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1)
  ) +
  ggplot2::labs(
    x = "Iteration",
    y = "AFAS score",
    title = "Aggregate Factor Alignment Similarity"
  ) +
  ggplot2::theme_bw(
    base_size = 13
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      hjust = 0.5,
      size = 14
    ),
    strip.background = ggplot2::element_rect(
      fill = "white",
      colour = "grey70"
    ),
    strip.text = ggplot2::element_text(
      size = 12
    )
  )

AFAS_trace


## =============================================================================
## Posterior distribution of AFAS
## =============================================================================

AFAS_density <- ggplot2::ggplot(
  afas_long,
  ggplot2::aes(
    x = AFAS,
    fill = chain,
    color = chain
  )
) +
  ggplot2::geom_density(
    alpha = 0.25,
    linewidth = 0.6
  ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 1)
  ) +
  ggplot2::labs(
    x = "AFAS score",
    y = "Density",
    title = "Posterior Distribution of AFAS"
  ) +
  ggplot2::theme_bw(
    base_size = 14
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      hjust = 0.5
    ),
    legend.title = ggplot2::element_blank()
  )

AFAS_density


## =============================================================================
## AFAS violin and boxplots
## =============================================================================

AFAS_violin <- ggplot2::ggplot(
  afas_long,
  ggplot2::aes(
    x = chain,
    y = AFAS,
    fill = chain
  )
) +
  ggplot2::geom_violin(
    trim = FALSE,
    alpha = 0.4
  ) +
  ggplot2::geom_boxplot(
    width = 0.15,
    outlier.size = 0.5
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1)
  ) +
  ggplot2::labs(
    x = "Chain",
    y = "AFAS score",
    title = "AFAS Across MCMC Chains"
  ) +
  ggplot2::theme_bw(
    base_size = 14
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      hjust = 0.5
    ),
    legend.position = "none"
  )

AFAS_violin


## =============================================================================
## Combined AFAS diagnostics
## =============================================================================

gridExtra::grid.arrange(
  AFAS_trace,
  AFAS_density,
  AFAS_violin,
  nrow = 2
)


## =============================================================================
## AFAS of the selected maximum log-posterior H
## =============================================================================

AFAS_H_max <- factor_similarity(
  H_ref = H_true,
  H_est = H_max
)$mean_similarity

AFAS_H_max

cat(
  "AFAS for the selected maximum log-posterior H:",
  round(AFAS_H_max, 3),
  "\n"
)

## =============================================================================
## Reconstruction performance across BMF methods
## Simulation Study I
## =============================================================================

X_hat_list <- list(
  Asso = X_asso,
  Sparse_BBMF = X_tilde_BBMF
)


## =============================================================================
## Compute reconstruction metrics
## =============================================================================

Metrics_X_ExptI <- do.call(
  rbind,
  lapply(
    names(X_hat_list),
    function(method) {
      
      metrics <- bmf_metrics(
        X = X_true,
        X_hat = X_hat_list[[method]]
      )
      
      data.frame(
        Method = method,
        t(metrics),
        row.names = NULL,
        check.names = FALSE
      )
    }
  )
)


## =============================================================================
## Round results for presentation
## =============================================================================

Metrics_X_ExptI_round <- Metrics_X_ExptI

numeric_cols <- vapply(
  Metrics_X_ExptI_round,
  is.numeric,
  logical(1)
)

Metrics_X_ExptI_round[
  numeric_cols
] <- round(
  Metrics_X_ExptI_round[
    numeric_cols
  ],
  digits = 3
)


## =============================================================================
## Display reconstruction performance
## =============================================================================

Metrics_X_ExptI_round

## =============================================================================
## Posterior reconstruction and uncertainty
## Simulation Study I
## =============================================================================

## -----------------------------------------------------------------------------
## Extract posterior W and H samples
## -----------------------------------------------------------------------------

Ws_ExptI <- extract_samples_bbmf(
  model_object = fit_BBMF,
  sample_name = "samples_W"
)

Hs_ExptI <- extract_samples_bbmf(
  model_object = fit_BBMF,
  sample_name = "samples_H"
)


## Inspect dimensions
lapply(
  Ws_ExptI,
  dim
)

lapply(
  Hs_ExptI,
  dim
)


## =============================================================================
## Reconstruct X for every posterior sample
## =============================================================================

X_reconstructed_ExptI <- lapply(
  seq_along(Ws_ExptI),
  function(chain_id) {
    
    reconstruct_X_one_chain_bbmf(
      W_samples = Ws_ExptI[[chain_id]],
      H_samples = Hs_ExptI[[chain_id]]
    )
  }
)

names(X_reconstructed_ExptI) <- names(fit_BBMF)


## Inspect dimensions
lapply(
  X_reconstructed_ExptI,
  dim
)


## =============================================================================
## Verify one posterior reconstruction
## =============================================================================

X_manual_check_ExptI <- bool_prod(
  W = Ws_ExptI[[1]][
    ,
    ,
    1
  ],
  H = Hs_ExptI[[1]][
    ,
    ,
    1
  ]
)

identical(
  X_reconstructed_ExptI[[1]][
    ,
    ,
    1
  ],
  X_manual_check_ExptI
)


## =============================================================================
## Posterior reconstruction summaries
## =============================================================================

X_summary_ExptI <- summarize_X_reconstruction_bbmf(
  reconstructed_chains = X_reconstructed_ExptI,
  X_truth = X_true
)


X_all_reconstructed_ExptI <-
  X_summary_ExptI$reconstructed_all

X_probability_ExptI <-
  X_summary_ExptI$probability

X_uncertainty_ExptI <-
  X_summary_ExptI$uncertainty


## Inspect posterior summaries
range(
  X_probability_ExptI
)

range(
  X_uncertainty_ExptI
)


## =============================================================================
## Prepare posterior summaries for plotting
## =============================================================================

X_long_ExptI <- make_X_long_data_bbmf(
  probability = X_probability_ExptI,
  uncertainty = X_uncertainty_ExptI,
  X_truth = X_true
)


probability_X_ExptI_df <-
  X_long_ExptI$probability

uncertainty_X_ExptI_df <-
  X_long_ExptI$uncertainty


## =============================================================================
## Posterior reconstruction probability
## =============================================================================

plot_X_probability_ExptI <- plot_X_posterior_bbmf(
  long_data = probability_X_ExptI_df,
  type = "probability",
  title = "Posterior Reconstruction Probabilities"
)

plot_X_probability_ExptI


## =============================================================================
## Posterior reconstruction uncertainty
## =============================================================================

plot_X_uncertainty_ExptI <- plot_X_posterior_bbmf(
  long_data = uncertainty_X_ExptI_df,
  type = "uncertainty",
  title = "Posterior Reconstruction Uncertainty"
)

plot_X_uncertainty_ExptI


## =============================================================================
## Display both posterior summaries
## =============================================================================

patchwork::wrap_plots(
  plot_X_probability_ExptI,
  plot_X_uncertainty_ExptI,
  ncol = 2
)

## =============================================================================
## Posterior inclusion probabilities and uncertainty for H
## =============================================================================

## Align posterior H samples within each chain
Hs_aligned_ExptI <- lapply(
  Hs_ExptI,
  function(H_chain) {
    
    align_H_samples_bbmf(
      H_ref = H_true,
      H_samples = H_chain
    )
  }
)

names(Hs_aligned_ExptI) <- names(Hs_ExptI)


## =============================================================================
## Posterior summaries
## =============================================================================

H_summary_ExptI <- summarize_H_posterior_bbmf(
  H_chains = Hs_aligned_ExptI,
  H_ref = H_true
)

H_probability_ExptI <-
  H_summary_ExptI$probability

H_uncertainty_ExptI <-
  H_summary_ExptI$uncertainty


## Inspect
round(
  H_probability_ExptI,
  2
)

round(
  H_uncertainty_ExptI,
  3
)


## =============================================================================
## Posterior inclusion probability plot
## =============================================================================

plot_H_probability_ExptI <- plot_H_posterior_bbmf(
  matrix_summary = H_probability_ExptI,
  type = "probability",
  title = "Posterior Inclusion Probabilities for H"
)

plot_H_probability_ExptI


## =============================================================================
## Posterior uncertainty plot
## =============================================================================

plot_H_uncertainty_ExptI <- plot_H_posterior_bbmf(
  matrix_summary = H_uncertainty_ExptI,
  type = "uncertainty",
  title = "Posterior Uncertainty for H"
)

plot_H_uncertainty_ExptI


## =============================================================================
## Display together
## =============================================================================

patchwork::wrap_plots(
  plot_H_probability_ExptI,
  plot_H_uncertainty_ExptI,
  ncol = 1
)
## =============================================================================
## End of Simulation Study I
## =============================================================================
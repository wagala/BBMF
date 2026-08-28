# BBMF

**Bayesian Boolean Matrix Factorization in R**

`BBMF` is an R package for Bayesian Boolean matrix factorization of binary data. The package implements the methodology presented in:

> Adolphus Wagala, Mehmet Samur, and Giovanni Parmigiani.  
> **A Bayesian Boolean Matrix Factorization with Application to Copy Number Analysis in Cancer.**  
> arXiv:2606.17491.

The methodology was motivated by the analysis of chromosomal copy number alterations in cancer, but the package can be applied to binary matrix data more generally.

## Overview

For a binary matrix \(X\) of dimension \(K \times G\), BBMF represents

\[
X \approx W \circ H,
\]

where

- \(W\) is a \(K \times R\) binary sample-factor matrix,
- \(H\) is an \(R \times G\) binary factor-feature matrix,
- \(\circ\) denotes Boolean matrix multiplication.

The package provides tools for:

- sparse Bayesian Boolean matrix factorization;
- posterior inference using Markov chain Monte Carlo;
- Boolean reconstruction;
- latent-factor alignment;
- Jaccard-based factor comparison;
- Aggregate Factor Alignment Similarity (AFAS);
- reconstruction diagnostics;
- posterior inclusion probabilities;
- posterior uncertainty quantification;
- visualization of binary latent structures.

## Installation

The development version can be installed from GitHub using:

```r
install.packages("remotes")

remotes::install_github("wagala/BBMF")
```

Then load the package with:

```r
library(BBMF)
```

## Basic Boolean matrix multiplication

The function `bool_prod()` computes a Boolean matrix product.

```r
W <- matrix(
  c(
    1, 0,
    0, 1,
    1, 1
  ),
  nrow = 3,
  byrow = TRUE
)

H <- matrix(
  c(
    1, 0, 1, 0,
    0, 1, 1, 0
  ),
  nrow = 2,
  byrow = TRUE
)

X <- bool_prod(W, H)

X
```

## Sparse BBMF

The main fitting functions are:

```r
run_sparseBBMF()
run_multiple_sparseBBMF()
```

`run_sparseBBMF()` fits a single MCMC chain, while `run_multiple_sparseBBMF()` can be used to fit multiple chains.

A typical workflow is:

```text
Observed binary matrix X
          |
          v
Initial Boolean factorization
          |
          v
Sparse BBMF
          |
          +----> MCMC diagnostics
          |
          +----> Posterior samples of W and H
          |          |
          |          +----> factor alignment
          |          +----> factor similarity / AFAS
          |          +----> posterior uncertainty
          |
          +----> Boolean reconstruction
                     |
                     +----> reconstruction metrics
                     +----> posterior probabilities
```

## Factor alignment and similarity

Because the ordering of latent factors is arbitrary, BBMF provides functions for aligning and comparing factor matrices:

```r
jaccard_rows()
jaccard_cols()
factor_alignment_hungarian()
factor_similarity()
plot_factor_similarity()
```

The package uses the Hungarian algorithm to obtain an optimal correspondence between latent factors.

## Reconstruction assessment

Boolean reconstructions can be evaluated using:

```r
bmf_metrics()
```

The reported measures include:

- specificity;
- F1 score;
- Matthews correlation coefficient;
- reconstruction error rate.

## Posterior uncertainty

BBMF provides functions for summarizing posterior uncertainty in both the reconstructed matrix \(X\) and the latent factor-feature matrix \(H\).

For \(X\):

```r
extract_samples_bbmf()
reconstruct_X_one_chain_bbmf()
summarize_X_reconstruction_bbmf()
plot_X_posterior_bbmf()
```

For \(H\):

```r
align_H_samples_bbmf()
summarize_H_posterior_bbmf()
plot_H_posterior_bbmf()
```

Posterior samples of \(H\) are aligned before computing posterior inclusion probabilities to account for factor-label switching.

## Vignette

A worked example demonstrating the main package functionality is provided in the package vignette:

**Getting Started with BBMF**

After installing the package, available vignettes can be listed using:

```r
vignette(package = "BBMF")
```

## Software dependencies

`BBMF` uses several existing R packages for supporting functionality.

In particular, [`rBMF`](https://CRAN.R-project.org/package=rBMF) is used for the Asso Boolean matrix factorization used for initialization and comparison.

The package also uses R packages including `ggplot2`, `clue`, `abind`, `reshape2`, `pbmcapply`, and others for computation, factor alignment, data processing, and visualization.

## Citation

If you use `BBMF`, please cite:

> Wagala, A., Samur, M., and Parmigiani, G. (2026).  
> **A Bayesian Boolean Matrix Factorization with Application to Copy Number Analysis in Cancer.**  
> arXiv:2606.17491.

## Authors

### Software

**Adolphus Wagala**

### Methodology

- Adolphus Wagala
- Mehmet Samur
- Giovanni Parmigiani

## License

This package is released under the MIT License.
# BBMF

**A General Nested Family of Bayesian Boolean Matrix Factorizations in R**

`BBMF` is an R package for Bayesian Boolean matrix factorization of binary data.

The package implements the methodology developed in:

> **A General Nested Family of Bayesian Boolean Matrix Factorizations with Application to Copy Number Analysis in Cancer**  
> Adolphus Wagala, Mehmet Samur, and Giovanni Parmigiani  
> 2026

The methodology was motivated by chromosome-arm copy-number alterations in cancer, but the framework is applicable to binary matrix data more generally.

## Overview

For an observed binary matrix

$$
X \in \{0,1\}^{K \times G},
$$

BBMF represents latent structure through

$$
Z = W \circ H,
$$

where

- $W \in \{0,1\}^{K \times R}$ is the sample-factor matrix,
- $H \in \{0,1\}^{R \times G}$ is the factor-feature matrix,
- $R$ is the factorization rank,
- $\circ$ denotes Boolean matrix multiplication.

For each entry,

$$
Z_{kg}
=
\bigvee_{r=1}^{R}
\left(
W_{kr} \wedge H_{rg}
\right).
$$

Thus, logical AND combines the sample and feature memberships within a factor, while logical OR combines contributions across factors.

## The nested BBMF family

`BBMF` provides five nested Bayesian Boolean matrix factorization models.

| Model | Description |
|---|---|
| **BBMF-Full** | Full hierarchical model with factor-specific sparsity, locus-level mixture selection, sample- and locus-specific heterogeneity, and asymmetric observation noise |
| **BBMF-NoPhi** | Removes factor-specific sparsity |
| **BBMF-Flat** | Additionally removes the locus-level mixture selection |
| **BBMF-Global** | Additionally replaces sample- and locus-specific activation probabilities with shared global parameters |
| **BBMF-Sym** | Additionally imposes symmetric observation noise |

The models form the sequence

$$
\text{BBMF-Full}
\rightarrow
\text{BBMF-NoPhi}
\rightarrow
\text{BBMF-Flat}
\rightarrow
\text{BBMF-Global}
\rightarrow
\text{BBMF-Sym}.
$$

Each successive model imposes one additional restriction while retaining the same Boolean factorization framework.

## Installation

Install the development version from GitHub with:

```r
install.packages("remotes")

remotes::install_github(
  "wagala/BBMF",
  build_vignettes = TRUE,
  dependencies = TRUE


# BBMF

**Bayesian Boolean Matrix Factorization in R**

`BBMF` is an R package for Bayesian Boolean matrix factorization of binary data.

The package implements the methods presented in:

> **A Bayesian Boolean Matrix Factorization with Application to Copy Number Analysis in Cancer**  
> Adolphus Wagala, Mehmet Samur, and Giovanni Parmigiani  
> arXiv:2606.17491

The methodology was motivated by the analysis of chromosomal copy number alterations in cancer, but the package can be applied to binary matrix data more generally.

## Overview

For a binary matrix $X$ of dimension $K \times G$, BBMF represents

$$
X \approx W \circ H,
$$

where

- $W$ is a $K \times R$ binary sample-factor matrix,
- $H$ is an $R \times G$ binary factor-feature matrix,
- $\circ$ denotes Boolean matrix multiplication.

In Boolean matrix multiplication, multiplication is replaced by logical AND and addition is replaced by logical OR.

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
- visualization of latent Boolean structure.

## Installation

The development version of `BBMF` can be installed from GitHub using:

```r
install.packages("remotes")

remotes::install_github(
  "wagala/BBMF",
  build_vignettes = TRUE,
  dependencies = TRUE
)
```

Then load the package:

```r
library(BBMF)
```

## Boolean matrix multiplication

The function `bool_prod()` computes the Boolean product of two binary matrices.

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

Here,

$$
X = W \circ H.
$$

## Sparse Bayesian Boolean Matrix Factorization

The main sparse BBMF fitting functions are:

```r
run_sparseBBMF()
run_multiple_sparseBBMF()
```

`run_sparseBBMF()` fits a single MCMC chain, while `run_multiple_sparseBBMF()` can be used to fit multiple chains.

A typical analysis follows the workflow

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
                     +----> posterior uncertainty
```

## Initialization with Asso

`BBMF` provides the function

```r
Asso_BBMF()
```

which uses the Asso Boolean matrix factorization implemented in the `rBMF` package.

The resulting factorization can be used to initialize BBMF or as a comparison method.

## MCMC diagnostics

Posterior sampling diagnostics can be obtained using:

```r
diagnosBBMF()
group_diagnost_List()
```

These functions summarize quantities such as the log-likelihood and log-posterior across MCMC samples and chains.

## Factor alignment

Latent factor labels are arbitrary. Therefore, factors should be aligned before comparing different factorizations or posterior samples.

BBMF provides:

```r
jaccard_rows()
jaccard_cols()
factor_alignment_hungarian()
factor_similarity()
plot_factor_similarity()
```

The package uses Jaccard similarity together with the Hungarian assignment algorithm to identify an optimal correspondence between latent factors.

## Aggregate Factor Alignment Similarity

The Aggregate Factor Alignment Similarity (AFAS) summarizes factor recovery after optimal alignment.

For an estimated factor matrix $H^{(m)}$, AFAS is defined as

```math
\mathrm{AFAS}^{(m)}
=
\frac{1}{R}
\sum_{r=1}^{R}
J\left(
H^{\mathrm{ref}}_{r\cdot},
H^{(m)}_{\sigma_m^*(r)\cdot}
\right).
```

Here, $J(\cdot,\cdot)$ denotes Jaccard similarity and $\sigma_m^*$ denotes the optimal factor assignment obtained using the Hungarian algorithm.

AFAS values closer to 1 indicate greater similarity between the estimated and reference latent factor structures.

Posterior AFAS values can be calculated using:

```r
get_mean_sim_per_chain()
```

## Reconstruction assessment

Boolean reconstruction performance can be evaluated using:

```r
bmf_metrics()
```

The function reports:

- specificity;
- F1 score;
- Matthews correlation coefficient;
- reconstruction error rate.

## Posterior reconstruction probabilities

Posterior samples of $W$ and $H$ can be used to reconstruct $X$ across the posterior distribution.

Relevant functions include:

```r
extract_samples_bbmf()
reconstruct_X_one_chain_bbmf()
summarize_X_reconstruction_bbmf()
make_X_long_data_bbmf()
plot_X_posterior_bbmf()
```

For entry $(k,g)$, the posterior reconstruction probability is

```math
\widehat{P}_{kg}
=
\frac{1}{M}
\sum_{m=1}^{M}
\widetilde{X}_{kg}^{(m)}.
```

Posterior uncertainty is summarized using

```math
U_{kg}
=
4\widehat{P}_{kg}
\left(
1-\widehat{P}_{kg}
\right).
```

Values close to 0 indicate strong posterior agreement, while values close to 1 indicate greater posterior uncertainty.

## Posterior inclusion probabilities for H

Because latent factor labels may switch across posterior samples, posterior samples of $H$ should be aligned before averaging.

BBMF provides:

```r
align_H_samples_bbmf()
summarize_H_posterior_bbmf()
plot_H_posterior_bbmf()
```

For factor $r$ and feature $g$, the posterior inclusion probability is

```math
\widehat{P}^{H}_{rg}
=
\frac{1}{M}
\sum_{m=1}^{M}
H_{rg}^{(m)}.
```

These summaries allow users to identify factor-feature associations that are consistently supported across posterior samples.

## Visualization

The package includes functions for visualizing:

- binary matrices;
- latent factor structures;
- aligned factor similarity;
- posterior reconstruction probabilities;
- posterior reconstruction uncertainty;
- posterior inclusion probabilities for $H$;
- posterior uncertainty in $H$.

For example:

```r
matrixPlot()
plot_factor_similarity()
plot_X_posterior_bbmf()
plot_H_posterior_bbmf()
```

## Main functions

| Function | Purpose |
|---|---|
| `bool_prod()` | Boolean matrix multiplication |
| `noise_bool()` | Add random noise to binary matrices |
| `Asso_BBMF()` | Obtain an Asso Boolean factorization |
| `run_sparseBBMF()` | Fit one sparse BBMF chain |
| `run_multiple_sparseBBMF()` | Fit multiple sparse BBMF chains |
| `diagnosBBMF()` | Compute MCMC diagnostic quantities |
| `group_diagnost_List()` | Organize diagnostics across chains |
| `matrixPlot()` | Visualize binary matrices |
| `jaccard_rows()` | Compute row-wise Jaccard similarities |
| `jaccard_cols()` | Compute column-wise Jaccard similarities |
| `factor_alignment_hungarian()` | Align latent factors |
| `factor_similarity()` | Compare factor matrices |
| `plot_factor_similarity()` | Visualize factor similarity |
| `get_mean_sim_per_chain()` | Compute AFAS across posterior samples |
| `bmf_metrics()` | Evaluate reconstruction performance |
| `extract_samples_bbmf()` | Extract posterior arrays |
| `reconstruct_X_one_chain_bbmf()` | Reconstruct $X$ across posterior samples |
| `summarize_X_reconstruction_bbmf()` | Summarize posterior reconstruction |
| `plot_X_posterior_bbmf()` | Plot posterior summaries for $X$ |
| `align_H_samples_bbmf()` | Align posterior samples of $H$ |
| `summarize_H_posterior_bbmf()` | Summarize posterior inclusion in $H$ |
| `plot_H_posterior_bbmf()` | Plot posterior summaries for $H$ |

## Vignette

A worked example demonstrating the main package functionality is included with the package.

After installation, available vignettes can be listed using:

```r
vignette(package = "BBMF")
```

The introductory vignette demonstrates:

- simulation of binary data;
- sparse BBMF fitting;
- MCMC diagnostics;
- latent-factor alignment;
- AFAS;
- reconstruction assessment;
- posterior reconstruction probabilities;
- posterior inclusion probabilities;
- uncertainty visualization.

## Software dependencies

`BBMF` builds on several existing R packages.

In particular, [`rBMF`](https://CRAN.R-project.org/package=rBMF) provides the Asso Boolean matrix factorization used for initialization and comparison.

Other packages support computation, factor alignment, data processing, parallel computation, and visualization, including:

- `abind`
- `clue`
- `ggplot2`
- `pbmcapply`
- `progress`
- `reshape2`

Additional packages such as `patchwork` are used for arranging visualizations in examples and vignettes.

## Citation

If you use `BBMF` in your research, please cite:

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

`BBMF` is released under the MIT License.
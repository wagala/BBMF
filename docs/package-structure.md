# BBMF Package Structure

The BBMF codebase is organized into four main components: core functions,
the standard BBMF model, the sparse BBMF model, and utilities for
simulation, visualization, initialization, and diagnostics.

```text
BBMF
│
├── Core functions
│   ├── bool_prod()
│   ├── compute_likelihood()
│   ├── compute_log_priorW()
│   ├── compute_log_priorH()
│   ├── compute_log_postW()
│   ├── compute_log_posteriorH()
│   └── sumLog()
│
├── Standard BBMF
│   ├── update_W()
│   ├── update_H()
│   ├── update_beta()
│   ├── run_BBMF()
│   ├── run_multiple_BBMF_sequential()
│   └── run_multiple_BBMF_parallel()
│
├── Sparse BBMF
│   ├── update_pi()
│   ├── update_psi()
│   ├── update_beta_g()
│   ├── update_W_Block_rows()
│   ├── update_H_Block_columns()
│   ├── run_sparseBBMF()
│   └── run_multiple_sparseBBMF()
│
└── Utilities and diagnostics
    ├── Asso_BBMF()
    ├── noise_bool()
    ├── matrixPlot()
    ├── diagnosBBMF()
    └── group_diagnost_List()
```

## Model hierarchy

The sparse BBMF implementation is the primary model in the package.

The standard BBMF implementation is retained mainly for baseline
comparisons, methodological evaluation, and reproducibility.

For most applications, the recommended workflow is:

```text
Observed binary matrix X
        │
        ▼
Asso_BBMF()
        │
        ├── W initial
        └── H initial
        │
        ▼
run_sparseBBMF()
        │
        ▼
Posterior samples
        │
        ├── W
        ├── H
        ├── alpha
        ├── beta
        ├── pi
        ├── psi
        ├── p11
        └── p10
        │
        ▼
diagnosBBMF()
        │
        ▼
Posterior diagnostics
```

## Supporting directories

```text
R/
    Reusable package functions.

scripts/
    Complete simulation and analysis scripts.

vignettes/
    User-facing Quarto tutorials and reproducible examples.

docs/
    Developer notes, package architecture, and implementation details.
```
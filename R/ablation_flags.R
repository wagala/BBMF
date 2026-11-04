#' Ablation Variant Flags
#'
#' Maps a variant name from the nested ablation family (Section 2.4)
#' to the four restriction switches that define it.
#'
#' @param variant One of \code{"Full"}, \code{"NoPhi"}, \code{"Flat"},
#'   \code{"Global"}, \code{"Sym"}.
#'
#' @return A list with logical elements \code{fix_phi},
#'   \code{disable_spike_slab}, \code{global_priors},
#'   \code{symmetric_noise}.
#'
#' @details
#' Each variant restricts the row above it in the nested family by one
#' further modeling choice:
#' \describe{
#'   \item{\code{Full}}{The proposed model; no restrictions.}
#'   \item{\code{NoPhi}}{\eqn{\phi_r \equiv 1} for all \eqn{r}
#'     (factor-level sparsity disabled).}
#'   \item{\code{Flat}}{Additionally \eqn{\psi_g \equiv 1} for all
#'     \eqn{g} (spike-and-slab disabled; \eqn{\beta_g} drawn only from
#'     the slab component).}
#'   \item{\code{Global}}{Additionally a single shared
#'     \eqn{\alpha_0,\beta_0} in place of per-sample/per-feature
#'     rates.}
#'   \item{\code{Sym}}{Additionally \eqn{p_{10} = 1-p_{11}} (single
#'     noise parameter); coincides with OrMachine's model up to
#'     differences in noise-parameter inference and sampler mechanics.}
#' }
#'
#' @export
ablation_flags <- function(variant = c("Full", "NoPhi", "Flat", "Global", "Sym")) {
  variant <- match.arg(variant)
  switch(
    variant,
    Full   = list(fix_phi = FALSE, disable_spike_slab = FALSE, global_priors = FALSE, symmetric_noise = FALSE),
    NoPhi  = list(fix_phi = TRUE,  disable_spike_slab = FALSE, global_priors = FALSE, symmetric_noise = FALSE),
    Flat   = list(fix_phi = TRUE,  disable_spike_slab = TRUE,  global_priors = FALSE, symmetric_noise = FALSE),
    Global = list(fix_phi = TRUE,  disable_spike_slab = TRUE,  global_priors = TRUE,  symmetric_noise = FALSE),
    Sym    = list(fix_phi = TRUE,  disable_spike_slab = TRUE,  global_priors = TRUE,  symmetric_noise = TRUE)
  )
}

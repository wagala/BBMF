#' Extract Posterior Samples from BBMF Chains
#'
#' Extracts a specified three-dimensional posterior sample array from
#' each BBMF MCMC chain. The initial stored state can optionally be
#' removed.
#'
#' @param model_object A list containing BBMF MCMC chain results.
#'   Each chain must contain the posterior sample array specified by
#'   \code{sample_name}.
#' @param sample_name Character string giving the name of the posterior
#'   sample array to extract, for example \code{"samples_W"} or
#'   \code{"samples_H"}.
#' @param drop_initial Logical. If \code{TRUE}, the first stored state
#'   is removed. Default is \code{TRUE}.
#'
#' @return A list of three-dimensional arrays, one for each MCMC chain.
#'
#' @details
#' The sparse BBMF sampler stores the initial state in the first position
#' of its posterior arrays. Setting \code{drop_initial = TRUE} removes
#' this state before posterior summaries are calculated.
#'
#' @examples
#' \dontrun{
#' Ws <- extract_samples_bbmf(
#'   model_object = fit_BBMF,
#'   sample_name = "samples_W"
#' )
#'
#' Hs <- extract_samples_bbmf(
#'   model_object = fit_BBMF,
#'   sample_name = "samples_H"
#' )
#' }
#'
#' @export
#'
extract_samples_bbmf <- function(
    model_object,
    sample_name,
    drop_initial = TRUE) {
  
  samples <- lapply(
    model_object,
    function(chain) {
      
      if (!sample_name %in% names(chain)) {
        stop(
          "The object ",
          sample_name,
          " was not found in one or more chains."
        )
      }
      
      sample_array <- chain[[sample_name]]
      
      if (length(dim(sample_array)) != 3L) {
        stop(
          sample_name,
          " must be a three-dimensional array."
        )
      }
      
      if (drop_initial) {
        
        if (dim(sample_array)[3] <= 1L) {
          stop(
            sample_name,
            " does not contain posterior samples after the initial state."
          )
        }
        
        sample_array <- sample_array[
          ,
          ,
          -1,
          drop = FALSE
        ]
      }
      
      sample_array
    }
  )
  
  names(samples) <- names(model_object)
  
  return(samples)
}
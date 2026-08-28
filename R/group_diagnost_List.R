#' Group BBMF Diagnostic Results Across MCMC Chains
#'
#' Rearranges a list of diagnostic data frames so that posterior
#' diagnostic quantities with the same name are grouped together
#' across MCMC chains.
#'
#' @param df_list A list of data frames, typically returned by
#'   \code{diagnosBBMF()}. Each data frame should have identical
#'   column names and the same number of rows.
#'
#' @return A named list of data frames. Each element corresponds to
#'   one diagnostic quantity in \code{df_list}. The columns of each
#'   returned data frame contain that diagnostic quantity for the
#'   different MCMC chains.
#'
#' @details
#' Suppose each data frame in \code{df_list} contains the columns
#' \code{log_likhood}, \code{log_postr}, \code{priorW},
#' \code{priorH}, and \code{BIC}. The function reorganizes the
#' results so that all chains for a particular diagnostic are stored
#' together.
#'
#' For example, the \code{log_likhood} element of the returned list
#' will contain columns named
#' \code{log_likhood_1}, \code{log_likhood_2}, and so on, where the
#' suffix identifies the MCMC chain.
#'
#' @return A list with one element for each column in the original
#'   diagnostic data frames. Each element is a data frame with one
#'   column per MCMC chain.
#'
#' @seealso
#' \code{\link{diagnosBBMF}}
#'
#' @export
#'
group_diagnost_List <- function(df_list) {
  
  # Get the column names from the first data frame
  column_names <- colnames(df_list[[1]])
  
  # Create a list to store the grouped data frames
  grouped_list <- lapply(
    seq_along(column_names),
    function(i) {
      
      group_df <- data.frame(
        do.call(
          cbind,
          lapply(
            df_list,
            function(df) df[[i]]
          )
        )
      )
      
      # Set column names indicating the diagnostic and chain
      colnames(group_df) <- paste(
        column_names[i],
        seq_along(df_list),
        sep = "_"
      )
      
      return(group_df)
    }
  )
  
  # Name each list element according to the diagnostic quantity
  names(grouped_list) <- column_names
  
  return(grouped_list)
}
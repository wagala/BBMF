#' Plot an Aligned Factor Similarity Heatmap
#'
#' Creates a heatmap of Jaccard similarities between reference and
#' estimated latent factors after optimal factor alignment.
#'
#' @param similarity_result A result returned by
#'   \code{factor_similarity()}.
#' @param x_label Character string specifying the x-axis label.
#' @param y_label Character string specifying the y-axis label.
#' @param title Optional plot title.
#' @param threshold Numeric value determining the switch between black
#'   and white text inside the heatmap cells. Default is \code{0.5}.
#' @param digits Number of decimal places displayed for similarity
#'   values. Default is \code{2}.
#'
#' @return A \code{ggplot2} object.
#'
#' @details
#' Similarity values greater than \code{threshold} are displayed with
#' white text, while smaller values are displayed with black text to
#' improve readability.
#'
#' @seealso
#' \code{\link{factor_similarity}}
#'
#' @export
#'
plot_factor_similarity <- function(
    similarity_result,
    x_label = "Estimated Factors",
    y_label = "Reference Factors",
    title = NULL,
    threshold = 0.5,
    digits = 2) {
  
  S <- similarity_result$similarity_aligned
  
  df <- reshape2::melt(
    S,
    varnames = c(
      "Reference_Factors",
      "Estimated_Factors"
    ),
    value.name = "Similarity"
  )
  
  # Text color for readability
  df$txt_col <- ifelse(
    df$Similarity > threshold,
    "white",
    "black"
  )
  
  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = Estimated_Factors,
      y = Reference_Factors,
      fill = Similarity
    )
  ) +
    ggplot2::geom_tile(
      color = "grey80"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = sprintf(
          paste0("%.", digits, "f"),
          Similarity
        ),
        color = txt_col
      ),
      size = 4
    ) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = "darkblue",
      limits = c(0, 1)
    ) +
    ggplot2::labs(
      x = x_label,
      y = y_label,
      title = title
    ) +
    ggplot2::theme_minimal(
      base_size = 14
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = 14
      ),
      axis.text.x = ggplot2::element_text(
        size = 14
      ),
      axis.text.y = ggplot2::element_text(
        size = 14
      ),
      legend.position = "none"
    ) +
    ggplot2::coord_fixed(
      ratio = 1
    )
}
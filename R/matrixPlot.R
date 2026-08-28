#' Plot a Binary Matrix as a Heatmap
#'
#' Creates a tile-based visualization of a binary matrix represented in
#' long format. Matrix entries are displayed using two fill colors.
#'
#' @param long_data A data frame containing the matrix data in long format.
#' @param x_value Unquoted column name in \code{long_data} to be displayed
#'   on the x-axis.
#' @param y_value Unquoted column name in \code{long_data} to be displayed
#'   on the y-axis.
#' @param value Unquoted column name containing the matrix entry values.
#' @param title Character string specifying the plot title. Defaults to
#'   an empty string.
#'
#' @return A \code{ggplot} object representing the matrix as a heatmap.
#'
#' @details
#' The function displays the data using tiles, reverses the ordering of
#' the y-axis, and uses a fixed aspect ratio so that each matrix entry is
#' represented by a square tile.
#'
#' Matrix values are treated as factors and displayed using cyan and
#' black fill colors.
#'
#'@export
#'
matrixPlot <- function(long_data, x_value, y_value, value, title = "") {
  
  ggplot2::ggplot(
    data = long_data,
    ggplot2::aes(
      x = {{ x_value }},
      y = forcats::fct_rev({{ y_value }}),
      fill = as.factor({{ value }})
    )
  ) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::labs(x = " ", y = " ", title = title) +
    ggplot2::scale_fill_manual(values = c("cyan", "black")) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 7),
      axis.text.x = ggplot2::element_text(
        size = 7,
        angle = 90,
        vjust = 0.5
      ),
      axis.text.y = ggplot2::element_text(size = 7),
      axis.text = ggplot2::element_text(size = 1),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 12),
      legend.text = ggplot2::element_text(size = 4),
      legend.title = ggplot2::element_blank(),
      legend.key.size = grid::unit(0.25, "cm"),
      legend.spacing.x = grid::unit(0.2, "cm"),
      legend.spacing.y = grid::unit(0.2, "cm")
    ) +
    ggplot2::coord_fixed(ratio = 1)
}
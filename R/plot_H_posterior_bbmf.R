#' Plot BBMF Posterior Summaries for H
#'
#' Creates a heatmap of posterior inclusion probabilities or posterior
#' uncertainty for the BBMF factor-feature matrix H.
#'
#' @param matrix_summary A numeric matrix containing posterior inclusion
#'   probabilities or uncertainty scores.
#' @param type Character string. Either \code{"probability"} or
#'   \code{"uncertainty"}.
#' @param title Optional plot title.
#' @param x_label Character string specifying the x-axis label.
#' @param y_label Character string specifying the y-axis label.
#' @param base_size Base font size. Default is \code{12}.
#'
#' @return A \code{ggplot2} object.
#'
#' @export
#'
plot_H_posterior_bbmf <- function(
    matrix_summary,
    type = c(
      "probability",
      "uncertainty"
    ),
    title = NULL,
    x_label = "Chromosome arm",
    y_label = "Latent factor",
    base_size = 12) {
  
  type <- match.arg(type)
  
  plot_df <- as.data.frame(
    as.table(matrix_summary),
    stringsAsFactors = FALSE
  )
  
  names(plot_df) <- c(
    "Factor",
    "Feature",
    "Value"
  )
  
  plot_df$Factor <- factor(
    plot_df$Factor,
    levels = rev(
      rownames(matrix_summary)
    )
  )
  
  plot_df$Feature <- factor(
    plot_df$Feature,
    levels = colnames(matrix_summary)
  )
  
  
  if (type == "probability") {
    
    if (is.null(title)) {
      title <- "Posterior Inclusion Probabilities for H"
    }
    
    fill_scale <- ggplot2::scale_fill_viridis_c(
      option = "C",
      limits = c(0, 1),
      breaks = seq(
        0,
        1,
        by = 0.25
      ),
      name = ""
    )
  }
  
  
  if (type == "uncertainty") {
    
    if (is.null(title)) {
      title <- "Posterior Uncertainty for H"
    }
    
    fill_scale <- ggplot2::scale_fill_viridis_c(
      option = "B",
      limits = c(0, 1),
      breaks = seq(
        0,
        1,
        by = 0.25
      ),
      name = ""
    )
  }
  
  
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = Feature,
      y = Factor,
      fill = Value
    )
  ) +
    ggplot2::geom_tile(
      color = "white",
      linewidth = 0.25
    ) +
    fill_scale +
    ggplot2::labs(
      title = title,
      x = x_label,
      y = y_label
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal(
      base_size = base_size
    ) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 9
      ),
      axis.text.y = ggplot2::element_text(
        size = 10
      ),
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = 12
      ),
      legend.position = "bottom"
    )
  
  return(p)
}
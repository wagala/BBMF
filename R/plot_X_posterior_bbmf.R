#' Plot BBMF Posterior Reconstruction Summaries
#'
#' Creates a heatmap of either posterior reconstruction probabilities
#' or posterior reconstruction uncertainty for a binary matrix
#' factorization.
#'
#' @param long_data A long-format data frame containing the columns
#'   \code{Sample}, \code{Feature}, and either \code{Probability} or
#'   \code{Uncertainty}. Such data can be produced using
#'   \code{make_X_long_data_bbmf()}.
#' @param type Character string specifying the quantity to plot.
#'   Must be either \code{"probability"} or \code{"uncertainty"}.
#' @param title Optional character string specifying the plot title.
#'   If \code{NULL}, an appropriate default title is used.
#' @param x_label Character string specifying the x-axis label.
#'   Default is \code{"Chromosome arm"}.
#' @param y_label Character string specifying the y-axis label.
#'   Default is \code{NULL}.
#' @param base_size Base font size used by the plot theme.
#'   Default is \code{12}.
#'
#' @return A \code{ggplot2} object.
#'
#' @details
#' For \code{type = "probability"}, the heatmap displays the posterior
#' probability that each reconstructed matrix entry equals one.
#'
#' For \code{type = "uncertainty"}, the heatmap displays the scaled
#' posterior reconstruction uncertainty
#'
#' \deqn{
#' U_{kg} =
#' 4\widehat{P}_{kg}(1-\widehat{P}_{kg}).
#' }
#'
#' Both quantities range from zero to one.
#'
#' @seealso
#' \code{\link{make_X_long_data_bbmf}},
#' \code{\link{summarize_X_reconstruction_bbmf}}
#'
#' @examples
#' \dontrun{
#' plot_X_posterior_bbmf(
#'   long_data = X_long$probability,
#'   type = "probability"
#' )
#'
#' plot_X_posterior_bbmf(
#'   long_data = X_long$uncertainty,
#'   type = "uncertainty"
#' )
#' }
#'
#' @export
#'
plot_X_posterior_bbmf <- function(
    long_data,
    type = c("probability", "uncertainty"),
    title = NULL,
    x_label = "Chromosome arm",
    y_label = NULL,
    base_size = 12) {
  
  type <- match.arg(type)
  
  
  ## ---------------------------------------------------------------------------
  ## Check required columns
  ## ---------------------------------------------------------------------------
  
  required_columns <- c(
    "Sample",
    "Feature"
  )
  
  if (!all(required_columns %in% names(long_data))) {
    stop(
      "long_data must contain the columns 'Sample' and 'Feature'."
    )
  }
  
  
  ## ---------------------------------------------------------------------------
  ## Posterior reconstruction probability
  ## ---------------------------------------------------------------------------
  
  if (type == "probability") {
    
    if (!"Probability" %in% names(long_data)) {
      stop(
        "For type = 'probability', long_data must contain ",
        "a 'Probability' column."
      )
    }
    
    if (is.null(title)) {
      title <- "Posterior Reconstruction Probabilities"
    }
    
    p <- ggplot2::ggplot(
      long_data,
      ggplot2::aes(
        x = Feature,
        y = Sample,
        fill = Probability
      )
    ) +
      ggplot2::geom_tile(
        color = "white",
        linewidth = 0.25
      ) +
      ggplot2::scale_fill_viridis_c(
        option = "C",
        limits = c(0, 1),
        breaks = seq(
          0,
          1,
          by = 0.25
        ),
        name = "",
        guide = ggplot2::guide_colorbar(
          title.position = "top",
          title.hjust = 0.5,
          barwidth = grid::unit(
            5,
            "cm"
          ),
          barheight = grid::unit(
            0.3,
            "cm"
          )
        )
      )
  }
  
  
  ## ---------------------------------------------------------------------------
  ## Posterior reconstruction uncertainty
  ## ---------------------------------------------------------------------------
  
  if (type == "uncertainty") {
    
    if (!"Uncertainty" %in% names(long_data)) {
      stop(
        "For type = 'uncertainty', long_data must contain ",
        "an 'Uncertainty' column."
      )
    }
    
    if (is.null(title)) {
      title <- "Posterior Reconstruction Uncertainty"
    }
    
    p <- ggplot2::ggplot(
      long_data,
      ggplot2::aes(
        x = Feature,
        y = Sample,
        fill = Uncertainty
      )
    ) +
      ggplot2::geom_tile(
        color = "white",
        linewidth = 0.25
      ) +
      ggplot2::scale_fill_viridis_c(
        option = "B",
        limits = c(0, 1),
        breaks = seq(
          0,
          1,
          by = 0.25
        ),
        name = "",
        guide = ggplot2::guide_colorbar(
          title.position = "top",
          title.hjust = 0.5,
          barwidth = grid::unit(
            5,
            "cm"
          ),
          barheight = grid::unit(
            0.3,
            "cm"
          )
        )
      )
  }
  
  
  ## ---------------------------------------------------------------------------
  ## Common plotting options
  ## ---------------------------------------------------------------------------
  
  p <- p +
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
        size = 8
      ),
      axis.text.y = ggplot2::element_text(
        size = 8
      ),
      plot.title = ggplot2::element_text(
        size = 12,
        hjust = 0.5
      ),
      legend.position = "bottom"
    )
  
  return(p)
}
#' make_statistic_dotplot: Statistics-derived dot plot
#'
#' @description A DotPlot visualization showing key statistics of a two-group
#' comparison
#'
#' @param data A dataframe containing the statistics of interest
#' @param x_axis A character string for the column name of the x-axis
#' @param y_axis A character string for the column name of the y-axis
#' @param size A character string for the column name defining size of the
#' bubbles/dots
#' @param color  A character string for the column name defining the color
#' @param palette A character string defining a color palette from RColorBrewer
#' or a list of colors defining low, middle (optional), and high values
#' @export
#'
#' @return A ggplot2 figure
#'
make_statistic_dotplot <- function(
  data,
  x_axis,
  y_axis,
  size,
  color,
  palette = "RdYlBu"
) {
  #identify range limits for color
  maxAbs <- NULL
  maxAbs <- max(abs(data[[color]]))
  # Create the plot using dynamic column names
  if (
    !is.list(palette) && (palette %in% rownames(RColorBrewer::brewer.pal.info))
  ) {
    dotplot <- ggplot2::ggplot(
      data,
      ggplot2::aes(
        x = .data[[x_axis]], # Dynamic x-axis column
        y = .data[[y_axis]], # Dynamic y-axis column
        size = .data[[size]], # Dynamic size column
        color = .data[[color]] # Dynamic color column
      )
    ) +
      ggplot2::geom_point() + # Add points (customize as needed)
      ggplot2::scale_color_distiller(
        palette = palette,
        limits = c(-maxAbs, maxAbs)
      ) + # Apply color palette
      ggplot2::theme_minimal() # Basic theme (customize as needed)
  } else {
    colPalette <- unlist(palette) #Turns palette into character string
    lowCol <- midCol <- highCol <- NULL
    if (length(colPalette) > 3) {
      stop("Too many colors submitted")
    } else if (length(colPalette) == 3) {
      lowCol <- colPalette[1]
      midCol <- colPalette[2]
      highCol <- colPalette[3]
    } else if (length(colPalette) == 2) {
      lowCol <- colPalette[1]
      highCol <- colPalette[2]
    } else {
      lowCol <- "white"
      highCol <- colPalette
    }
    if (is.null(midCol)) {
      dotplot <- ggplot2::ggplot(
        data,
        ggplot2::aes(
          x = .data[[x_axis]], # Dynamic x-axis column
          y = .data[[y_axis]], # Dynamic y-axis column
          size = .data[[size]], # Dynamic size column
          color = .data[[color]] # Dynamic color column
        )
      ) +
        ggplot2::geom_point() + # Add points (customize as needed)
        ggplot2::scale_color_gradient(
          low = lowCol,
          high = highCol,
          limits = c(-maxAbs, maxAbs)
        ) + # Apply color palette
        ggplot2::theme_minimal() # Basic theme (customize as needed)
    } else {
      dotplot <- ggplot2::ggplot(
        data,
        ggplot2::aes(
          x = .data[[x_axis]], # Dynamic x-axis column
          y = .data[[y_axis]], # Dynamic y-axis column
          size = .data[[size]], # Dynamic size column
          color = .data[[color]] # Dynamic color column
        )
      ) +
        ggplot2::geom_point() + # Add points (customize as needed)
        ggplot2::scale_color_gradient2(
          low = lowCol,
          mid = midCol,
          high = highCol,
          limits = c(-maxAbs, maxAbs)
        ) + # Apply color palette
        ggplot2::theme_minimal() # Basic theme (customize as needed)
    }
  }
  return(dotplot)
}

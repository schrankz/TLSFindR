#' Plot TLS niche overlay on spatial image
#'
#' @description
#' Create a spatial overlay that highlights TLS candidate niches vs other niches.
#'
#' @param obj Seurat object with `Niche_Type` metadata.
#' @param structure_label Character. Label used to identify the TLS-rich niche. Default "TLS niche".
#' @param tissue_overlay Logical. If TRUE, use SpatialDimPlot (requires histology).
#' @param structure_color Character. Color of structure_label on spatial image. Default "#D81B60".
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' niche_image <- make_niche_plot(obj, 
#' structure_label = "TLS niche", 
#' tissue_overlay = TRUE, 
#' structure_color = "blue")
#' }
#'
#' @export
#' @importFrom Seurat ImageDimPlot SpatialDimPlot Images
#' @importFrom ggplot2 theme_minimal theme labs
make_niche_plot <- function(
  obj,
  structure_label = "TLS niche",
  tissue_overlay = FALSE,
  structure_color = "#D81B60"
) {

  # Group labels
  obj$Labels <- ifelse(
    obj$Niche_Type == structure_label,
    structure_label,
    "Other niches"
  )

  lab_levels <- sort(unique(as.character(obj$Labels)))
  colors <- setNames(rep("gray80", length(lab_levels)), lab_levels)
  colors[structure_label] <- structure_color

  # Decide which plotting function to use
  has_images <- length(Seurat::Images(obj)) > 0

  if (tissue_overlay && has_images) {
    p <- Seurat::SpatialDimPlot(
      obj,
      group.by = "Labels",
      cols = colors,
      pt.size.factor = 1
    )
  } else if (tissue_overlay && !has_images) {
    warning("tissue_overlay = TRUE but no spatial image found. Falling back to ImageDimPlot().")
    p <- Seurat::ImageDimPlot(
      obj,
      group.by = "Labels",
      size = 1,
      cols = colors
    )
  } else {
    p <- Seurat::ImageDimPlot(
      obj,
      group.by = "Labels",
      size = 1,
      cols = colors
    )
  }

  # Common theming
  p +
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      legend.text = element_text(color = "black"),
      legend.title = element_text(color = "black")
    ) +
    labs(color = "Niche_Type")
}
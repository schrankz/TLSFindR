#' Plot heatmap of cell type proportions across niches
#'
#' @description
#' Create a heatmap showing cell type proportions across niches.
#'
#' @param seurat_obj A Seurat object with metadata columns specified by `niche_slot` and `label_slot`.
#' @param niche_slot Character. Metadata column for niche labels (default 'Niche_Type').
#' @param label_slot Character. Metadata column for cell type labels (default 'SingleR_label').
#' @param normalize_by Character. One of "none", "row", or "column" to normalize counts accordingly.
#' @param wrap_labels Logical. Whether to wrap row and column labels onto multiple lines.
#' @param wrap_width Integer. Character width used for wrapping labels.
#' @param wrap_by Character. How to wrap labels: "space", "underscore", or "custom".
#' @param wrap_pattern Character. Custom regex pattern used when wrap_by = "custom".
#' @param row_fontsize Numeric. Font size for row labels.
#' @param column_fontsize Numeric. Font size for column labels.
#' @param max_label_chars Integer or NULL. Truncate labels longer than this length.
#' @param row_label_fun Function or NULL. Custom function to transform row labels.
#' @param column_label_fun Function or NULL. Custom function to transform column labels.
#'
#' @return A ComplexHeatmap object (invisible return).
#'
#' @examples
#' \dontrun{
#' plot <- plot_niche_heatmaps(obj, 
#' niche_slot = "Niche_Type", 
#' label_slot = "SingleR_label", 
#' normalize_by = "row")
#' }
#'
#' @export
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom grid gpar
plot_niche_heatmaps <- function(
  seurat_obj,
  niche_slot = "Niche_Type",
  label_slot = "SingleR_label",
  normalize_by = c("none", "row", "column"),
  row_fontsize = 12,
  column_fontsize = 12,
  max_label_chars = NULL,
  wrap_labels = FALSE,
  wrap_width = 20,
  wrap_by = c("space", "underscore", "custom"),
  wrap_pattern = NULL,
  row_label_fun = NULL,
  column_label_fun = NULL
) {
  normalize_by <- match.arg(normalize_by)
  meta <- seurat_obj@meta.data

  # Build contingency table
  comp <- table(meta[[niche_slot]], meta[[label_slot]])

  # Normalization
  if (normalize_by == "row") {
    comp <- prop.table(comp, margin = 1) * 100
  } else if (normalize_by == "column") {
    comp <- prop.table(comp, margin = 2) * 100
  }

  # --- Label processing ----
  wrap_text <- function(x, width, by, pattern = NULL) {
  if (by == "space") {
    sapply(x, function(lbl) {
      paste(strwrap(lbl, width = width), collapse = "\n")
    })
  } else if (by == "underscore") {
    sapply(x, function(lbl) {
      paste(strsplit(lbl, "_")[[1]], collapse = "\n")
    })
  } else if (by == "custom" && !is.null(pattern)) {
    sapply(x, function(lbl) {
      paste(strsplit(lbl, pattern)[[1]], collapse = "\n")
    })
  } else {
    x
  }
}
  row_labels <- rownames(comp)
  col_labels <- colnames(comp)
  
# Truncation
  if (!is.null(max_label_chars)) {
    truncate <- function(x) {
      ifelse(
        nchar(x) > max_label_chars,
        paste0(substr(x, 1, max_label_chars - 1), "..."),
        x
      )
    }
    row_labels <- truncate(row_labels)
    col_labels <- truncate(col_labels)
  }
  
  # Wrapping
  if (wrap_labels) {
    wrap_by <- match.arg(wrap_by)
    row_labels <- wrap_text(row_labels, wrap_width, wrap_by, wrap_pattern)
    col_labels <- wrap_text(col_labels, wrap_width, wrap_by, wrap_pattern)
  }

  if (!is.null(row_label_fun)) {
    row_labels <- row_label_fun(row_labels)
  }

  if (!is.null(column_label_fun)) {
    col_labels <- column_label_fun(col_labels)
  }

  # Create heatmap
  ht <- ComplexHeatmap::Heatmap(
    comp,
    name = switch(
      normalize_by,
      "none" = "Count",
      "row" = "% Cells in Niche",
      "column" = "% Cells by Niche"
    ),
    row_labels = row_labels,
    column_labels = col_labels,
    row_names_gp = grid::gpar(fontsize = row_fontsize),
    column_names_gp = grid::gpar(fontsize = column_fontsize)
  )

  return(ht)
}
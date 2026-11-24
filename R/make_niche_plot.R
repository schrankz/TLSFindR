#' Plot TLS niche overlay on spatial image
#'
#' @description
#' Create a spatial overlay that highlights TLS candidate niches vs other niches.
#' If TLS were counted with `summarize_niches(count_TLS = TRUE)` and added with
#' `add_niche_metadata()`, the plot can instead show only the counted TLS, show
#' fragments excluded by `min_TLS_cells` in a separate color, and label each TLS
#' with its number.
#'
#' @param obj Seurat object with `Niche_Type` metadata. For `plot_mode = "counted"`
#'   or `label_TLS = TRUE`, also `TLS_number` and `TLS_status` metadata (added by
#'   `add_niche_metadata()` from `summarize_niches(count_TLS = TRUE)` output).
#' @param structure_label Character. Label used to identify the TLS-rich niche. Default "TLS niche".
#' @param tissue_overlay Logical. If TRUE, use SpatialDimPlot (requires histology).
#' @param structure_color Character. Color of the TLS niche (or counted TLS) on the spatial
#'   image. Default "#D81B60".
#' @param plot_mode Character. `"niche"` (default) colors every cell in the TLS niche.
#'   `"counted"` colors only cells in counted TLS, i.e. structures with at least
#'   `min_TLS_cells` cells.
#' @param show_excluded Logical. Only used when `plot_mode = "counted"`. If TRUE, TLS-niche
#'   fragments that were too small to count are shown in `excluded_color`; if FALSE
#'   (default) they are shown as "Other niches".
#' @param excluded_color Character. Color for excluded fragments. Default "#FFC107".
#' @param other_color Character. Color for all other cells. Default "gray80".
#' @param label_TLS Logical. If TRUE, write each counted TLS's number at its center.
#'   Numbers match the `TLS_number` metadata column and `TLS_summary` table
#'   (1 = largest TLS in that image). Default FALSE.
#' @param label_size Numeric. Text size of TLS number labels. Default 3.
#' @param label_color Character. Text color of TLS number labels. Default "black".
#' @param label_fill Character or NA. Background color of the label box, which keeps
#'   numbers readable on top of colored cells. Use NA for plain text. Default "white".
#' @param image Character or NULL. Name of the image/FOV to plot. If NULL (default),
#'   Seurat's default image is used (all images for `SpatialDimPlot`).
#' @param pt_size Numeric. Point size (`size` in ImageDimPlot, `pt.size.factor` in
#'   SpatialDimPlot). Default 1.
#'
#' @return A ggplot object (or a patchwork object if several images are plotted).
#'
#' @examples
#' \dontrun{
#' # Original behavior: highlight the whole TLS niche
#' niche_image <- make_niche_plot(obj,
#' structure_label = "TLS niche",
#' tissue_overlay = TRUE,
#' structure_color = "blue")
#'
#' # After counting: show only counted TLS, excluded fragments in yellow, numbered
#' res <- summarize_niches(obj, cell_label = "SingleR_label", niche_label = "niches",
#'                         Bcell_label = "B cells", count_TLS = TRUE, min_TLS_cells = 20)
#' obj <- add_niche_metadata(obj, res)
#' make_niche_plot(obj, plot_mode = "counted", show_excluded = TRUE, label_TLS = TRUE)
#' }
#'
#' @export
#' @importFrom Seurat ImageDimPlot SpatialDimPlot Images
#' @importFrom ggplot2 theme_minimal theme labs element_line element_rect element_text
#'   geom_text geom_label aes ggplot_build unit
make_niche_plot <- function(
  obj,
  structure_label = "TLS niche",
  tissue_overlay = FALSE,
  structure_color = "#D81B60",
  plot_mode = c("niche", "counted"),
  show_excluded = FALSE,
  excluded_color = "#FFC107",
  other_color = "gray80",
  label_TLS = FALSE,
  label_size = 3,
  label_color = "black",
  label_fill = "white",
  image = NULL,
  pt_size = 1
) {

  plot_mode <- match.arg(plot_mode)
  meta <- obj@meta.data

  if (!"Niche_Type" %in% colnames(meta)) {
    stop("No `Niche_Type` metadata found. Run add_niche_metadata() first.")
  }

  needs_counts <- plot_mode == "counted" || label_TLS
  if (needs_counts && !"TLS_number" %in% colnames(meta)) {
    stop("plot_mode = \"counted\" and label_TLS = TRUE need TLS counts. Run ",
         "summarize_niches(..., count_TLS = TRUE), then add_niche_metadata() with its output.")
  }
  if (plot_mode == "counted" && show_excluded && !"TLS_status" %in% colnames(meta)) {
    warning("No `TLS_status` metadata found, so excluded fragments cannot be shown. ",
            "Re-run summarize_niches(count_TLS = TRUE) and add_niche_metadata() with ",
            "this version of TLSFindR.", call. = FALSE)
    show_excluded <- FALSE
  }
  if (plot_mode == "niche" && show_excluded) {
    message("show_excluded is only used with plot_mode = \"counted\"; ignoring it.")
    show_excluded <- FALSE
  }

  # ---- Group labels and colors ----
  other_lab <- "Other niches"
  if (plot_mode == "niche") {
    tls_lab <- structure_label
    is_tls <- !is.na(meta$Niche_Type) & meta$Niche_Type == structure_label
    labs_vec <- ifelse(is_tls, tls_lab, other_lab)
    colors <- stats::setNames(c(structure_color, other_color), c(tls_lab, other_lab))
  } else {
    tls_lab <- "Counted TLS"
    excl_lab <- "Excluded fragments"
    is_tls <- !is.na(meta$TLS_number)
    labs_vec <- ifelse(is_tls, tls_lab, other_lab)
    if (show_excluded) {
      is_excl <- !is.na(meta$TLS_status) & meta$TLS_status == "Excluded fragment"
      labs_vec[is_excl & !is_tls] <- excl_lab
    }
    colors <- stats::setNames(c(structure_color, excluded_color, other_color),
                              c(tls_lab, excl_lab, other_lab))
  }
  colors <- colors[names(colors) %in% labs_vec]
  obj$Labels <- factor(labs_vec, levels = names(colors))

  # ---- Plot ----
  has_images <- length(Seurat::Images(obj)) > 0
  if (tissue_overlay && !has_images) {
    warning("tissue_overlay = TRUE but no spatial image found. Falling back to ImageDimPlot().")
    tissue_overlay <- FALSE
  }

  if (tissue_overlay) {
    args <- list(obj, group.by = "Labels", cols = colors, pt.size.factor = pt_size)
    if (!is.null(image)) args$images <- image
    p <- do.call(Seurat::SpatialDimPlot, args)
  } else {
    args <- list(obj, group.by = "Labels", size = pt_size, cols = colors)
    if (!is.null(image)) args$fov <- image
    p <- do.call(Seurat::ImageDimPlot, args)
  }

  # ---- Common theming ----
  style <- list(
    ggplot2::theme_minimal(),
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "gray90"),
      panel.grid.minor = ggplot2::element_line(color = "gray95"),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.key = ggplot2::element_rect(fill = "white", color = NA),
      legend.text = ggplot2::element_text(color = "black"),
      legend.title = ggplot2::element_text(color = "black")
    ),
    ggplot2::labs(color = "Niche_Type", fill = "Niche_Type")
  )

  tls_num <- stats::setNames(meta$TLS_number, rownames(meta))
  finish <- function(sp) {
    sp <- sp + style
    if (label_TLS) sp <- .add_tls_number_labels(sp, tls_num, label_size, label_color, label_fill)
    sp
  }

  if (inherits(p, "patchwork")) {
    for (i in seq_along(p)) p[[i]] <- finish(p[[i]])
    p
  } else {
    finish(p)
  }
}

#' Add TLS numbers at the center of each counted TLS in one spatial ggplot
#'
#' Positions are taken from the plot itself (via `ggplot_build()`), so labels line
#' up regardless of axis flipping or image scaling done by Seurat.
#'
#' @keywords internal
#' @noRd
.add_tls_number_labels <- function(p, tls_num, label_size, label_color, label_fill = "white") {
  built <- ggplot2::ggplot_build(p)

  pos <- NULL
  for (i in seq_along(p$layers)) {
    ld <- p$layers[[i]]$data
    if (!is.data.frame(ld) || nrow(ld) == 0) ld <- p$data
    if (!is.data.frame(ld) || nrow(ld) == 0) next
    cells <- if ("cell" %in% colnames(ld)) as.character(ld$cell) else rownames(ld)
    # ImageDimPlot prefixes cell names with the boundary name (e.g. "centroids_cell1")
    if (!any(cells %in% names(tls_num)) && "boundary" %in% colnames(ld)) {
      cells <- substring(cells, nchar(as.character(ld$boundary)) + 2)
    }
    bd <- built$data[[i]]
    if (!all(c("x", "y") %in% colnames(bd)) || nrow(bd) != length(cells)) next
    if (!any(cells %in% names(tls_num))) next
    pos <- data.frame(cell = cells, x = bd$x, y = bd$y, stringsAsFactors = FALSE)
    break
  }

  if (is.null(pos)) {
    warning("Could not find point positions in the plot, so TLS numbers were not added.",
            call. = FALSE)
    return(p)
  }

  pos$TLS_number <- tls_num[pos$cell]
  pos <- pos[!is.na(pos$TLS_number), , drop = FALSE]
  if (nrow(pos) == 0) return(p)

  centers <- stats::aggregate(cbind(x, y) ~ TLS_number, data = pos, FUN = stats::median)

  # Snap each label to the TLS cell nearest the median, so it sits on the TLS even if curved
  centers[, c("x", "y")] <- t(vapply(seq_len(nrow(centers)), function(i) {
    sub <- pos[pos$TLS_number == centers$TLS_number[i], ]
    j <- which.min((sub$x - centers$x[i])^2 + (sub$y - centers$y[i])^2)
    c(sub$x[j], sub$y[j])
  }, numeric(2)))

  if (is.na(label_fill)) {
    p + ggplot2::geom_text(
      data = centers,
      mapping = ggplot2::aes(x = x, y = y, label = TLS_number),
      inherit.aes = FALSE, size = label_size, color = label_color, fontface = "bold"
    )
  } else {
    p + ggplot2::geom_label(
      data = centers,
      mapping = ggplot2::aes(x = x, y = y, label = TLS_number),
      inherit.aes = FALSE, size = label_size, color = label_color, fill = label_fill,
      fontface = "bold", label.padding = ggplot2::unit(0.12, "lines")
    )
  }
}

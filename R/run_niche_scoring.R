#' Score niches using signatures via GSVA/ssGSEA
#'
#' @description
#' Compute enrichment scores (ssGSEA/GSVA) for provided gene signatures on niche-aggregated expression profiles and generate either a bar or dot plot depicting scores per niche.
#'
#' @param seurat_obj A Seurat object.
#' @param signatures A named list of gene vectors (signatures).
#' @param assay Character. Which assay to use (e.g., "Xenium" or "SCT").
#' @param slot Character. Slot to pull from (default "data").
#' @param make_plot Logical. Whether to return a ggplot2 object (bar or dot plot) depicting scores across niches.
#' @param plot_title Character. Title for plot.
#' @param niche_label Character. Metadata column with niche assignments. Default "niches".
#' @param labeled Logical. Add text labels to points in plot.
#' @param plot_type Character. Type of plot to generate. 
#'   Either "bar" or "dot". Default is "bar".
#'
#' @param palette Optional named character vector of colors to use for niches.
#'   If NULL (default), ggplot2 default discrete colors are used.
#'   Named vectors are recommended to ensure correct color mapping.
#'
#' @return A list with `ssgsea_results` (data.frame) and `plot` (ggplot or NULL).
#'
#' @examples
#' \dontrun{
#' out <- run_niche_scoring(obj, human_tls_signatures, assay = "SCT", slot = "data")
#' out$plot
#' }
#'
#' @export
#' @importFrom GSVA gsva
#' @importFrom Matrix rowMeans
#' @importFrom Seurat GetAssayData
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_point geom_col scale_size 
#'   scale_fill_manual scale_color_manual
#'   facet_wrap theme_minimal ylab xlab ggtitle geom_text
#' @importFrom tidytext reorder_within scale_x_reordered
#' @importFrom dplyr %>%
run_niche_scoring <- function(
  seurat_obj,
  signatures,
  assay = c("Xenium", "SCT"),
  slot = "data",
  make_plot = TRUE,
  plot_type = c("bar", "dot"),
  plot_title = "TLS Signatures Per Niche",
  niche_label = "niches",# default
  labeled = FALSE,
  palette = NULL
) {
  
  # Check if niche column exists
  if (!(niche_label %in% colnames(seurat_obj@meta.data))) {
    stop(paste("No column named", niche_label, "found in Seurat metadata."))
  }
  
  plot_type <- match.arg(plot_type)
  
  expr <- GetAssayData(seurat_obj, assay = assay, slot = slot)
  # Extract clean vector of niche identifiers
  niches_vec <- seurat_obj@meta.data[[niche_label]] %>% as.character()
  niches <- unique(niches_vec)
  
  niche_expr_list <- lapply(niches, function(n) {
    cells <- colnames(seurat_obj)[niches_vec == n]
    if (length(cells) > 2) {
      Matrix::rowMeans(expr[, cells, drop = FALSE])
    } else {
      NULL
    }
  })
  names(niche_expr_list) <- niches
  niche_expr_list <- niche_expr_list[!sapply(niche_expr_list, is.null)]

  expr_mat <- do.call(cbind, niche_expr_list)

  ssgsea_param <- GSVA::ssgseaParam(expr_mat, geneSets = signatures)
  ssgsea_scores <- GSVA::gsva(ssgsea_param)

  ssgsea_df <- as.data.frame(t(ssgsea_scores))
  ssgsea_df$niches <- rownames(ssgsea_df)

  ssgsea_long <- reshape2::melt(
    ssgsea_df,
    id.vars = "niches",
    variable.name = "Signature",
    value.name = "NES"
  )

  plot <- NULL
  if (make_plot) {

    if (plot_type == "bar") {

        plot <- ggplot(
         ssgsea_long,
         aes(
            y = tidytext::reorder_within(niches, NES, Signature),
            x = NES,
            fill = niches
    )
    ) +
        geom_col(show.legend = FALSE) +
         tidytext::scale_y_reordered() +
         facet_wrap(~Signature, scales = "free_y") +
            theme_minimal() +
            theme(
             axis.text.y = element_text(size = 9),
             axis.title.x = element_text(size = 11),
             axis.title.y = element_text(size = 11),
             strip.text = element_text(size = 12, face = "bold")
            ) +
            ylab("Niche") +
            xlab("Normalized Enrichment Score") +
            ggtitle(plot_title)
    
        if (!is.null(palette) && length(palette) > 0) {
  
            niche_levels <- unique(ssgsea_long$niches)
  
            if (!is.null(names(palette))) {
             # Named palette → match to niches
             palette_use <- palette[niche_levels]
          } else {
            # Unnamed → assign in order
              palette_use <- palette[seq_along(niche_levels)]
             names(palette_use) <- niche_levels
             }
        if (!is.null(palette) && length(palette) < length(unique(ssgsea_long$niches))) {
        warning("Provided palette has fewer colors than niches.")
        }
  
         plot <- plot + scale_fill_manual(values = palette_use)
        }
    }


    if (plot_type == "dot") {
        plot <- ggplot(ssgsea_long, aes(x = NES, y = reorder(Signature, NES))) +
         geom_point(aes(color = niches), size = 3) +
        theme_minimal() +
        ylab("Signature") +
        xlab("Normalized Enrichment Score") +
        ggtitle(plot_title)
    
        if (!is.null(palette) && length(palette) > 0) {
  
        niche_levels <- unique(ssgsea_long$niches)
  
        if (!is.null(names(palette))) {
            # Named palette → match to niches
            palette_use <- palette[niche_levels]
        } else {
         # Unnamed → assign in order
            palette_use <- palette[seq_along(niche_levels)]
            names(palette_use) <- niche_levels
            }
        if (!is.null(palette) && length(palette) < length(unique(ssgsea_long$niches))) {
        warning("Provided palette has fewer colors than niches.")
        }
  
        plot <- plot + scale_fill_manual(values = palette_use)
        }
    }
}
  return(list(ssgsea_results = ssgsea_long, plot = plot))
}
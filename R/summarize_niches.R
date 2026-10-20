#' Summarize niche composition and identify candidate TLS niches
#'
#' @description
#' Summarize cell type counts and proportions per niche and classify niches based on dominant and subdominant cell type content. Candidate TLS niches are identified
#' based on a specified proportion of B cells or whether the primary cell type is B cells.
#'
#' @param obj Seurat object with `obj[[niche_label]]` and `obj[[cell_label]]` accessible as metadata or assays.
#' @param cell_label Character. Metadata column name with cell labels (e.g., "SingleR_label").
#' @param niche_label Character. Name of niche assignment column (e.g., "niches").
#' @param Bcell_label Character. Label to use for B cells in the cell_label column. Default "B cells".
#' @param B_cutoff Numeric. Minimum proportion of B cells to call a niche a "TLS niche". Default 0.1.
#' @param structure_label Label to call the B cell-rich niche. Default "TLS niche"
#' @param max_subdominant Integer (0–5). Maximum number of subdominant cell types to report per niche.
#' @param count_TLS Logical. If TRUE, count physically separate TLS structures within the
#'   TLS niche(s). Cells in TLS niches are linked when they lie within a linking distance of
#'   each other, and each separate group of linked cells is one structure. Default FALSE.
#' @param min_TLS_cells Integer. Minimum number of cells (or spots, for spot-based data) a
#'   separate structure must contain to be counted as a TLS. Smaller fragments are excluded
#'   and reported in `n_fragments_excluded`. Only used when `count_TLS = TRUE`. Default 20.
#' @param TLS_max_dist Numeric or NULL. Maximum distance, in the object's coordinate units,
#'   between two cells for them to be part of the same structure. If NULL (default), it is set
#'   per image to `TLS_dist_multiplier` times the median nearest-neighbor distance between all
#'   cells/spots in that image.
#' @param TLS_dist_multiplier Numeric. Multiplier used to set the linking distance when
#'   `TLS_max_dist` is NULL. Larger values bridge bigger gaps and merge nearby structures;
#'   smaller values split structures at small gaps. Default 2.
#' @param image Character vector or NULL. Image/FOV name(s) to count TLS in. If NULL
#'   (default), all images in the object are used, and each is counted separately.
#'
#' @return If `count_TLS = FALSE`, a data.frame summarizing niche-level information including
#'   candidate TLS identification, dominant, and subdominant cell types and proportions.
#'
#'   If `count_TLS = TRUE`, a list with:
#'   \itemize{
#'     \item `niche_summary`: the niche-level data.frame described above.
#'     \item `TLS_counts`: one row per image with the number of TLS (`n_TLS`), the number of
#'       fragments below `min_TLS_cells` (`n_fragments_excluded`), and the linking distance used.
#'     \item `TLS_summary`: one row per counted TLS with its size, B cell count and proportion,
#'       centroid coordinates, and the niche(s) it came from.
#'     \item `TLS_cells`: the cell-to-TLS assignments (`cell`, `TLS_ID`).
#'   }
#'   This list can be passed directly to `add_niche_metadata()`, which will also add a
#'   `TLS_ID` metadata column.
#'
#' @examples
#' \dontrun{
#' summary <- summarize_niches(obj, 
#' cell_label = "SingleR_label", 
#' niche_label = "niches", 
#' Bcell_label = "B cells")
#'
#' # Also count physically separate TLS, keeping structures with at least 20 cells
#' res <- summarize_niches(obj,
#' cell_label = "SingleR_label",
#' niche_label = "niches",
#' Bcell_label = "B cells",
#' count_TLS = TRUE,
#' min_TLS_cells = 20)
#' res$TLS_counts
#' res$TLS_summary
#' }
#'
#' @export
#' @importFrom dplyr group_by summarise mutate arrange ungroup
#' @importFrom tidyr pivot_wider
#' @importFrom glue glue
summarize_niches <- function(
  obj,
  cell_label,
  niche_label,
  Bcell_label,
  B_cutoff = 0.1,
  structure_label = "TLS niche",
  max_subdominant = 1,
  count_TLS = FALSE,
  min_TLS_cells = 20,
  TLS_max_dist = NULL,
  TLS_dist_multiplier = 2,
  image = NULL
) {

  stopifnot(max_subdominant >= 0, max_subdominant <= 5)
  stopifnot(is.logical(count_TLS), length(count_TLS) == 1)
  if (count_TLS) {
    stopifnot(is.numeric(min_TLS_cells), min_TLS_cells >= 1)
    stopifnot(is.null(TLS_max_dist) || (is.numeric(TLS_max_dist) && TLS_max_dist > 0))
    stopifnot(is.numeric(TLS_dist_multiplier), TLS_dist_multiplier > 0)
  }

  niche_info <- data.frame(
    Niche = obj[[niche_label]][,1],
    CellType = obj[[cell_label]][,1]
  )

  # Compute proportions
  niche_props <- niche_info %>%
    dplyr::group_by(Niche, CellType) %>%
    dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(Niche) %>%
    dplyr::mutate(Proportion = Count / sum(Count)) %>%
    dplyr::arrange(Niche, dplyr::desc(Proportion))

  niche_labels <- niche_props %>%
    dplyr::group_by(Niche) %>%
    dplyr::mutate(
      Bcell_Prop = sum(Proportion[CellType == Bcell_label]),
      Rank = dplyr::row_number()
    ) %>%
    dplyr::summarise(
      Dominant = CellType[1],
      DominantProp = Proportion[1],
      Subdominant = list(CellType[Rank > 1 & Rank <= (1 + max_subdominant)]),
      SubdominantProp = list(Proportion[Rank > 1 & Rank <= (1 + max_subdominant)]),
      Bcell_Prop = dplyr::first(Bcell_Prop),
      .groups = "drop"
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      Niche_Type = dplyr::case_when(
        Dominant == Bcell_label | Bcell_Prop >= B_cutoff ~ structure_label,
        length(Subdominant[[1]]) > 0 ~
          glue::glue(
            "{Dominant}-dominant; {paste(glue::glue('{Subdominant}-subdominant'), collapse = ', ')} niche"
          ),
        TRUE ~ glue::glue("{Dominant}-dominant niche")
      )
    ) %>%
    dplyr::ungroup()

  # Join labels back to proportions for wide output
  niche_summary <- niche_props %>%
    tidyr::pivot_wider(
      names_from = CellType,
      values_from = c(Proportion, Count),
      names_glue = "{CellType}_{.value}",
      values_fill = 0
    ) %>%
    dplyr::left_join(niche_labels, by = "Niche") %>%
    dplyr::distinct(Niche, .keep_all = TRUE)

  if (!count_TLS) {
    return(niche_summary)
  }

  # ---- Count physically separate TLS structures ----
  cells <- rownames(obj[[niche_label]])
  niches_vec <- stats::setNames(as.character(obj[[niche_label]][, 1]), cells)
  cell_types <- stats::setNames(as.character(obj[[cell_label]][, 1]), cells)

  tls_niches <- as.character(niche_labels$Niche[niche_labels$Niche_Type == structure_label])
  if (length(tls_niches) == 0) {
    message("No niches were labeled '", structure_label, "', so no TLS were counted.")
  }
  tls_cells <- cells[niches_vec %in% tls_niches]

  tls <- .count_tls_structures(
    obj = obj,
    tls_cells = tls_cells,
    cell_types = cell_types,
    niches = niches_vec,
    Bcell_label = Bcell_label,
    min_cells = min_TLS_cells,
    max_dist = TLS_max_dist,
    dist_multiplier = TLS_dist_multiplier,
    images = image
  )

  message(sprintf(
    "Counted %d TLS (>= %d cells each) across %d image(s); %d smaller fragment(s) excluded.",
    sum(tls$TLS_counts$n_TLS), min_TLS_cells, nrow(tls$TLS_counts),
    sum(tls$TLS_counts$n_fragments_excluded)
  ))

  list(
    niche_summary = niche_summary,
    TLS_counts = tls$TLS_counts,
    TLS_summary = tls$TLS_summary,
    TLS_cells = tls$TLS_cells
  )
}
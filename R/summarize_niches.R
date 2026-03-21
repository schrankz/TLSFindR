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
#'
#' @return A data.frame summarizing niche-level information including candidate TLS identification, dominant, and subdominant cell types and proportions.
#'
#' @examples
#' \dontrun{
#' summary <- summarize_niches(obj, 
#' cell_label = "SingleR_label", 
#' niche_label = "niches", 
#' Bcell_label = "B cells")
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
  max_subdominant = 1
) {

  stopifnot(max_subdominant >= 0, max_subdominant <= 5)

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
      Bcell_Prop = first(Bcell_Prop),
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

  return(niche_summary)
}
#' Add niche summary information to Seurat metadata
#'
#' @description
#' Merge niche-level summary information back into cell-level metadata and add to the Seurat object.
#'
#' @param obj A Seurat object.
#' @param niche_summary A data.frame returned by `summarize_niches()` with (at minimum) Niche and Niche_Type columns.
#'
#' @return Seurat object with added metadata columns.
#'
#' @examples
#' \dontrun{
#' obj <- add_niche_metadata(obj, niche_summary)
#' }
#'
#' @export
#' @importFrom dplyr select
#' @importFrom Seurat AddMetaData
add_niche_metadata <- function(obj, niche_summary){
  cell_niches <- obj[["niches"]]  %>%
  as.data.frame()
  colnames(cell_niches) <- "niches"
  trimmed_niche_summary <- dplyr::select(niche_summary, Niche_Type, Niche)
  niches_meta <- left_join(cell_niches, trimmed_niche_summary, by = c("niches" = "Niche"))

  obj <- AddMetaData(obj, niches_meta)
  return(obj)
}
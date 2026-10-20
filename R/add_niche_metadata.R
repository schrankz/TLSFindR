#' Add niche summary information to Seurat metadata
#'
#' @description
#' Merge niche-level summary information back into cell-level metadata and add to the Seurat object.
#'
#' @param obj A Seurat object.
#' @param niche_summary The output of `summarize_niches()`: either the niche-level data.frame
#'   (with at minimum Niche and Niche_Type columns), or the list returned when
#'   `count_TLS = TRUE`. In the latter case a `TLS_ID` column is also added, giving the
#'   structure each cell belongs to (NA for cells not in a counted TLS).
#'
#' @return Seurat object with added metadata columns.
#'
#' @examples
#' \dontrun{
#' obj <- add_niche_metadata(obj, niche_summary)
#' }
#'
#' @export
#' @importFrom dplyr select left_join
#' @importFrom Seurat AddMetaData
add_niche_metadata <- function(obj, niche_summary){
  # Accept the list returned by summarize_niches(count_TLS = TRUE)
  tls_cells <- NULL
  if (is.list(niche_summary) && !is.data.frame(niche_summary)) {
    tls_cells <- niche_summary$TLS_cells
    niche_summary <- niche_summary$niche_summary
  }

  cell_niches <- obj[["niches"]]  %>%
  as.data.frame()
  colnames(cell_niches) <- "niches"
  trimmed_niche_summary <- dplyr::select(niche_summary, Niche_Type, Niche)
  niches_meta <- dplyr::left_join(cell_niches, trimmed_niche_summary, by = c("niches" = "Niche"))

  obj <- AddMetaData(obj, niches_meta)

  if (!is.null(tls_cells)) {
    tls_id <- stats::setNames(rep(NA_character_, ncol(obj)), colnames(obj))
    tls_id[tls_cells$cell] <- tls_cells$TLS_ID
    obj <- AddMetaData(obj, metadata = tls_id, col.name = "TLS_ID")
  }

  return(obj)
}
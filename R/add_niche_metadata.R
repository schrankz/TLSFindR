#' Add niche summary information to Seurat metadata
#'
#' @description
#' Merge niche-level summary information back into cell-level metadata and add to the Seurat object.
#'
#' @param obj A Seurat object.
#' @param niche_summary The output of `summarize_niches()`: either the niche-level data.frame
#'   (with at minimum Niche and Niche_Type columns), or the list returned when
#'   `count_TLS = TRUE`. In the latter case three more columns are added:
#'   `TLS_ID` (e.g. "fov_TLS_1"; NA for cells not in a counted TLS), `TLS_number`
#'   (the TLS number within its image, 1 = largest), and `TLS_status`
#'   ("Counted TLS", "Excluded fragment", or NA for cells outside TLS niches).
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
  excluded_cells <- NULL
  if (is.list(niche_summary) && !is.data.frame(niche_summary)) {
    tls_cells <- niche_summary$TLS_cells
    excluded_cells <- niche_summary$excluded_cells
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

    tls_num <- stats::setNames(rep(NA_integer_, ncol(obj)), colnames(obj))
    if (!is.null(tls_cells$TLS_number)) {
      tls_num[tls_cells$cell] <- as.integer(tls_cells$TLS_number)
    } else {
      # Output from an older version: recover the number from the ID
      tls_num[tls_cells$cell] <- as.integer(sub(".*_TLS_", "", tls_cells$TLS_ID))
    }
    obj <- AddMetaData(obj, metadata = tls_num, col.name = "TLS_number")

    tls_status <- stats::setNames(rep(NA_character_, ncol(obj)), colnames(obj))
    if (!is.null(excluded_cells) && nrow(excluded_cells) > 0) {
      tls_status[excluded_cells$cell] <- "Excluded fragment"
    }
    tls_status[tls_cells$cell] <- "Counted TLS"
    obj <- AddMetaData(obj, metadata = tls_status, col.name = "TLS_status")
  }

  return(obj)
}
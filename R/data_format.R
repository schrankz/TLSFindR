#' Convert SpatialExperiment to Seurat
#'
#' @description
#' Convert a `Seurat` or `SpatialExperiment` object to a Seurat object.
#' If the input is already a `Seurat` object, it is returned unchanged.
#'
#' @param obj An R object. Expected to be a `Seurat` or `SpatialExperiment`.
#'
#' @return A `Seurat` object.
#'
#' @examples
#' \dontrun{
#' seurat_obj <- data_format(my_spatial_experiment)
#' }
#'
#' @export
#' @importFrom Seurat as.Seurat
data_format <- function(obj){

  check_pkgs <- function(pkgs){
    missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if(length(missing)){
      stop(sprintf("Missing required packages: %s.",
                   paste(missing, collapse=", ")))
    }
    invisible(TRUE)
  }

  check_pkgs(c("Seurat", "SpatialExperiment"))

  # If already Seurat, return unchanged
  if (inherits(obj, "Seurat")) {
    return(obj)
  }

  # Only allow SpatialExperiment conversion
  if (inherits(obj, "SpatialExperiment")) {
    seurat_obj <- tryCatch(
      Seurat::as.Seurat(obj, assay = NULL, data = "counts"),
      error = function(e) Seurat::as.Seurat(obj)
    )
    return(seurat_obj)
  }

  stop("Input must be a Seurat or SpatialExperiment object.")
}
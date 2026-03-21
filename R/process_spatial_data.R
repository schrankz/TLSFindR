#' Normalize and scale spatial transcriptomic data
#'
#' @description
#' Normalize and scale a Seurat object using either SCTransform or log-normalize + FindVariableFeatures + ScaleData.
#'
#' @param seurat_obj A Seurat object.
#' @param assay Character. Name of the assay to normalize (e.g., "Xenium" or "SCT").
#' @param method Character. One of "sctransform" or "lognorm". Default "sctransform".
#'
#' @return The normalized and scaled Seurat object.
#'
#' @examples
#' \dontrun{
#' obj <- process_spatial_data(seurat_obj, assay = "Xenium", method = "sctransform")
#' }
#'
#' @export
#' @importFrom Seurat SCTransform NormalizeData FindVariableFeatures ScaleData
process_spatial_data <- function(seurat_obj, assay, method = c("sctransform","lognorm")){
  method <- match.arg(method)
  if(method == "sctransform"){
    seurat_obj <- Seurat::SCTransform(seurat_obj, assay = assay, verbose = FALSE)
    print("Using SCTransform")
  } else {
    seurat_obj <- Seurat::NormalizeData(seurat_obj, assay = assay, verbose = FALSE)
    seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, verbose = FALSE)
    seurat_obj <- Seurat::ScaleData(seurat_obj, verbose = FALSE)
    print("Using Log Normalization")
  }
  return(seurat_obj)
}
#' Build niche assay using supplied groupings
#'
#' @description
#' Convenience wrapper to build niche/neighbor-based assay using available BuildNicheAssay function (expected to exist in the environment).
#'
#' @param seurat_obj A Seurat object.
#' @param optimal_k Numeric. Number of neighborhoods/niches to consider (`niches.k` passed to BuildNicheAssay).
#' @param num_neighbors Numeric. Number of cell/spot neighbors to consider per cell (`neighbors.k` passed to BuildNicheAssay). Default 30.
#' @param group_by Character. Metadata column to group cells by when building niches (e.g., "SingleR_label").
#' @param fov Character. Column name containing the field-of-view within the spatial Seurat object. Default "fov".
#'
#' @return Seurat object with niche assay added.
#'
#' @examples
#' \dontrun{
#' obj <- run_niche_analysis(seurat_obj, 
#' optimal_k = 10, 
#' num_neighbors = 30, 
#' group_by = "SingleR_label", 
#' fov = "fov")
#' }
#'
#' @export
run_niche_analysis <- function(seurat_obj, optimal_k, num_neighbors = 30, group_by, fov = "fov"){
 seurat_obj <- BuildNicheAssay(object = seurat_obj, group.by = group_by, fov = fov,
    niches.k = optimal_k, neighbors.k = num_neighbors)
 return(seurat_obj)
}

#' Annotate cells using SingleR and celldex reference
#'
#' @description
#' Run SingleR annotation on expression data extracted from a spatial Seurat object using a celldex reference dataset.
#'
#' @param seurat_obj A Seurat object.
#' @param reference Character. One of "HumanPrimaryCellAtlasData", "BlueprintEncodeData", "ImmGenData", "MouseRNAseqData", "MonacoImmuneData".
#' @param assay Character or name. The assay to pull expression from (e.g., "Xenium", "SCT").
#' @param labels_name Character. Metadata column name to hold assigned labels. Default "SingleR_label".
#'
#' @return Seurat object with new metadata column containing SingleR labels.
#'
#' @examples
#' \dontrun{
#' obj <- run_annotation(seurat_obj, reference = "HumanPrimaryCellAtlasData", assay = "SCT")
#' }
#'
#' @export
#' @importFrom SingleR SingleR
#' @importFrom celldex HumanPrimaryCellAtlasData BlueprintEncodeData ImmGenData MouseRNAseqData MonacoImmuneData
#' @importFrom Seurat GetAssayData
run_annotation <- function(seurat_obj, reference = "HumanPrimaryCellAtlasData", assay = c("Xenium", "SCT"), labels_name = "SingleR_label"){
  check_pkgs <- function(pkgs){
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if(length(missing)){
    stop(sprintf("Missing required packages: %s.",
                 paste(missing, collapse=", ")))
  }
  invisible(TRUE)
}
  check_pkgs(c("SingleR","celldex","Seurat"))
  expr <- as.matrix(Seurat::GetAssayData(seurat_obj, assay, slot = "data"))
  ref_fun <- switch(reference,
                    HumanPrimaryCellAtlasData = celldex::HumanPrimaryCellAtlasData,
                    BlueprintEncodeData = celldex::BlueprintEncodeData,
                    ImmGenData = celldex::ImmGenData,
                    MouseRNAseqData = celldex::MouseRNAseqData,
                    MonacoImmuneData = celldex::MonacoImmuneData,
                    stop("Unknown reference"))
  ref <- ref_fun()
  pred <- SingleR::SingleR(test = expr, ref = ref, labels = ref$label.main)
  seurat_obj@meta.data[[labels_name]] <- pred$labels
  return(seurat_obj)
}
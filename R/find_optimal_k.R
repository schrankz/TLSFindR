#' Find optimal number of neighborhoods K using cell type mixing and transcriptional coherence
#'
#' @description
#' Evaluate k-nearest neighbor values by balancing transcriptional coherence with
#' cell type mixing (diversity). This approach is designed to identify biologically
#' meaningful niches where different cell types interact, such as T-B cell 
#' co-localization. The optimal k maximizes both transcriptional similarity and
#' cell type diversity within neighborhoods.
#'
#' @param seurat_obj A Seurat object with spatial coordinates accessible via `GetTissueCoordinates`.
#' @param k_seq Numeric vector. Sequence of k values to test. If NULL (default), will
#'   automatically start from the number of unique cell types in `group_by` and go up
#'   to `k_max` by `k_step`. If provided, this overrides the automatic calculation.
#' @param assay Character. Assay to use for expression data. Default "SCT".
#' @param n_features Integer. Number of highly variable features to use. Default 2000.
#' @param slot Character. Slot to pull expression data from. Default "data".
#' @param group_by Character. Metadata column containing cell/spot type labels 
#'   (e.g., "SingleR_label"). Required for calculating cell type mixing scores.
#'   Default "SingleR_label".
#' @param k_max Integer. Maximum k value to test when auto-generating k_seq.
#'   Only used if k_seq is NULL. Default 50.
#' @param k_step Integer. Step size for k sequence when auto-generating k_seq.
#'   Only used if k_seq is NULL. Default 1.
#' @param coherence_weight Numeric. Weight for transcriptional coherence score (0-1).
#'   Default 0.5 gives equal weight to coherence and mixing. Higher values (e.g., 0.7)
#'   emphasize transcriptional similarity; lower values (e.g., 0.3) emphasize cell type
#'   diversity.
#'
#' @return A list containing:
#' \itemize{
#'   \item optimal_k: The selected k value
#'   \item coherence_scores: Vector of transcriptional coherence scores for each k
#'   \item mixing_scores: Vector of cell type mixing scores for each k
#'   \item combined_scores: Vector of weighted combined scores for each k
#'   \item k_values: The k values tested
#' }
#'
#' @details
#' For each k value, computes spatial k-nearest neighbors and calculates two scores:
#' 
#' \strong{1. Transcriptional coherence:} Pearson correlation between each cell's 
#' expression profile and the mean expression of its spatial neighbors. Higher values
#' indicate spatially proximal cells are transcriptionally similar.
#' 
#' \strong{2. Cell type mixing:} Shannon entropy of cell type composition in each 
#' neighborhood. Higher entropy indicates more diverse cell type mixing, which can
#' reveal functional niches where multiple cell types interact.
#' 
#' The two scores are normalized to 0-1, then combined using `coherence_weight`:
#' 
#' \code{combined_score = coherence_weight * coherence + (1 - coherence_weight) * mixing}
#' 
#' The optimal k is chosen at the elbow point (maximum second derivative) of the
#' combined score curve, representing the point where adding more neighbors provides
#' diminishing returns.
#' 
#' \strong{When k_seq is NULL:} The function automatically determines the starting k
#' based on the number of unique cell types in `group_by`. For example, if you have
#' 8 cell types, k will range from 8 to `k_max` (default 50) by `k_step` (default 1).
#' This ensures neighborhoods are large enough to capture meaningful interactions
#' between multiple cell types.
#' 
#' This approach is particularly useful for:
#' \itemize{
#'   \item Identifying immune cell niches (e.g., T-B cell interactions)
#'   \item Detecting tumor microenvironments with diverse cell populations
#'   \item Finding spatial domains in low-resolution data (e.g., Visium)
#'   \item Characterizing tissue organization patterns
#' }
#'
#' @examples
#' \dontrun{
#' # Auto-start from number of unique cell types
#' result <- find_optimal_k(seurat_obj, group_by = "cell_type")
#' bestk <- result$optimal_k
#' 
#' # Use custom k range
#' result <- find_optimal_k(seurat_obj, 
#'                          group_by = "cell_type",
#'                          k_seq = seq(5, 30, 2))
#' 
#' # Auto-start but customize max and step
#' result <- find_optimal_k(seurat_obj,
#'                          group_by = "cell_type", 
#'                          k_max = 100,
#'                          k_step = 5)
#' 
#' # Emphasize cell type mixing over coherence
#' result <- find_optimal_k(seurat_obj,
#'                          group_by = "cell_type",
#'                          coherence_weight = 0.3)
#' 
#' # Visualize the optimization
#' plot(result$k_values, result$combined_scores, 
#'      type = "b", pch = 19,
#'      xlab = "k (neighborhood size)", 
#'      ylab = "Combined score",
#'      main = "Optimal k selection")
#' abline(v = result$optimal_k, col = "red", lty = 2)
#' 
#' # Compare individual scores
#' par(mfrow = c(1, 3))
#' plot(result$k_values, result$coherence_scores, type = "b",
#'      main = "Transcriptional Coherence", xlab = "k", ylab = "Coherence")
#' plot(result$k_values, result$mixing_scores, type = "b",
#'      main = "Cell Type Mixing", xlab = "k", ylab = "Mixing (entropy)")
#' plot(result$k_values, result$combined_scores, type = "b",
#'      main = "Combined Score", xlab = "k", ylab = "Score")
#' abline(v = result$optimal_k, col = "red", lty = 2)
#' 
#' # Use optimal k in BuildNicheAssay
#' seurat_obj <- BuildNicheAssay(seurat_obj,
#'                                fov = "fov",
#'                                niches.k = result$optimal_k,
#'                                neighbors.k = result$optimal_k)
#' }
#'
#' @export
#' @importFrom FNN get.knn
#' @importFrom Seurat GetTissueCoordinates GetAssayData VariableFeatures FindVariableFeatures
find_optimal_k <- function(seurat_obj, 
                           k_seq = NULL,
                           assay = "SCT",
                           n_features = 2000,
                           slot = "data",
                           group_by = "SingleR_label",
                           k_max = 50,
                           k_step = 1,
                           coherence_weight = 0.5) {
  
  # Check required packages
  check_pkgs <- function(pkgs) {
    missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing)) {
      stop(sprintf("Missing required packages: %s.",
                   paste(missing, collapse = ", ")))
    }
    invisible(TRUE)
  }
  check_pkgs(c("FNN", "Seurat"))
  
  # Check that group_by is provided
  if (is.null(group_by)) {
    stop("Cell/spot labels must be provided in 'group_by' for cell type mixing analysis.")
  }
  
  if (!group_by %in% colnames(seurat_obj@meta.data)) {
    stop(sprintf("Column '%s' not found in seurat_obj metadata.", group_by))
  }
  
  # Set k_seq based on number of unique cell types if not provided
  if (is.null(k_seq)) {
    n_cell_types <- length(unique(seurat_obj@meta.data[[group_by]]))
    k_seq <- seq(n_cell_types, k_max, k_step)
    message(sprintf("Starting k from %d unique cell types in '%s'", 
                    n_cell_types, group_by))
  }
  
  # Get spatial coordinates and cell types
  coords <- as.matrix(GetTissueCoordinates(seurat_obj)[, 1:2])
  cell_types <- seurat_obj@meta.data[[group_by]]
  
  # Get expression data for highly variable features
  all_hvf <- VariableFeatures(seurat_obj, assay = assay)
  if (length(all_hvf) == 0) {
    message("No variable features found. Running FindVariableFeatures()...")
    seurat_obj <- FindVariableFeatures(seurat_obj, assay = assay, verbose = FALSE)
    all_hvf <- VariableFeatures(seurat_obj, assay = assay)
  }
  hvf <- head(all_hvf, n_features)
  
  expr_data <- GetAssayData(seurat_obj, assay = assay, slot = slot)
  expr_data <- as.matrix(expr_data[hvf, ])
  
  # Compute coherence and mixing scores for each k
  coherence_scores <- numeric(length(k_seq))
  mixing_scores <- numeric(length(k_seq))
  
  for (i in seq_along(k_seq)) {
    k <- k_seq[i]
    
    # Get spatial k-nearest neighbors
    knn_result <- FNN::get.knn(coords, k = k)
    
    # Calculate transcriptional coherence
    coherence <- vapply(seq_len(nrow(coords)), function(j) {
      neighbors <- knn_result$nn.index[j, ]
      neighbor_mean <- rowMeans(expr_data[, neighbors, drop = FALSE])
      cor(expr_data[, j], neighbor_mean, use = "complete.obs")
    }, numeric(1))
    
    coherence_scores[i] <- mean(coherence, na.rm = TRUE)
    
    # Calculate cell type mixing (Shannon entropy)
    mixing <- vapply(seq_len(nrow(coords)), function(j) {
      neighbors <- knn_result$nn.index[j, ]
      neighbor_types <- cell_types[neighbors]
      type_counts <- table(neighbor_types)
      type_props <- type_counts / sum(type_counts)
      # Shannon entropy - higher means more mixing
      -sum(type_props * log(type_props + 1e-10))
    }, numeric(1))
    
    mixing_scores[i] <- mean(mixing, na.rm = TRUE)
  }
  
  # Normalize scores to 0-1 range
  coherence_norm <- (coherence_scores - min(coherence_scores)) / 
                    (max(coherence_scores) - min(coherence_scores) + 1e-10)
  mixing_norm <- (mixing_scores - min(mixing_scores)) / 
                 (max(mixing_scores) - min(mixing_scores) + 1e-10)
  
  # Combine scores with weights
  combined_scores <- coherence_weight * coherence_norm + 
                     (1 - coherence_weight) * mixing_norm
  
  # Find elbow point using second derivative
  optimal_k <- if (length(combined_scores) >= 3) {
    second_deriv <- diff(diff(combined_scores))
    k_seq[which.max(abs(second_deriv)) + 1]
  } else {
    # Fallback if not enough k values
    k_seq[which.max(combined_scores)]
  }
  
  # Return results
  list(
    optimal_k = optimal_k,
    coherence_scores = coherence_scores,
    mixing_scores = mixing_scores,
    combined_scores = combined_scores,
    k_values = k_seq
  )
}
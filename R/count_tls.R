# Internal helpers for counting physically separate TLS structures.
# Used by summarize_niches(count_TLS = TRUE). Not exported.

#' Get x/y coordinates and cell names for one image
#'
#' Handles the column layouts returned by `GetTissueCoordinates()` for
#' FOV-based images (Xenium, CosMx, MERFISH: columns x, y, cell) and
#' Visium images (x, y, cell in VisiumV2; imagerow, imagecol in VisiumV1).
#'
#' @keywords internal
#' @noRd
.get_image_coords <- function(obj, image) {
  co <- as.data.frame(Seurat::GetTissueCoordinates(obj, image = image))
  cells <- if ("cell" %in% colnames(co)) as.character(co$cell) else rownames(co)

  if (all(c("x", "y") %in% colnames(co))) {
    xy <- co[, c("x", "y")]
  } else if (all(c("imagecol", "imagerow") %in% colnames(co))) {
    xy <- co[, c("imagecol", "imagerow")]
  } else {
    num_cols <- which(vapply(co, is.numeric, logical(1)))
    if (length(num_cols) < 2) {
      stop(sprintf("Could not find x/y coordinates for image '%s'.", image))
    }
    xy <- co[, num_cols[1:2]]
  }

  out <- data.frame(cell = cells, x = as.numeric(xy[[1]]), y = as.numeric(xy[[2]]),
                    stringsAsFactors = FALSE)
  out[stats::complete.cases(out), , drop = FALSE]
}

#' Label connected components of points that lie within `radius` of each other
#'
#' Two points are linked when they are within `radius`; points linked through a
#' chain of such links share a component. DBSCAN with `minPts = 1` gives
#' exactly these components (every point is a core point), using a k-d tree so
#' it scales to large sections.
#'
#' @return Integer vector of component ids, one per row of `coords`.
#' @keywords internal
#' @noRd
.connected_components <- function(coords, radius) {
  n <- nrow(coords)
  if (n == 0) return(integer(0))
  if (n == 1) return(1L)
  as.integer(dbscan::dbscan(as.matrix(coords), eps = radius, minPts = 1)$cluster)
}

#' Count physically separate TLS structures within TLS niches
#'
#' @param obj Seurat object with spatial coordinates.
#' @param tls_cells Character vector of cells belonging to TLS niches.
#' @param cell_types Named character vector of cell labels (names = cells).
#' @param niches Named vector of niche ids (names = cells).
#' @param Bcell_label Label used for B cells.
#' @param min_cells Minimum cells/spots for a structure to be counted.
#' @param max_dist Linking distance in coordinate units, or NULL for automatic.
#' @param dist_multiplier Multiplier on the median nearest-neighbour distance
#'   used when `max_dist` is NULL.
#' @param images Image names to use, or NULL for all images.
#'
#' @return A list with `TLS_summary`, `TLS_counts` and `TLS_cells`.
#' @keywords internal
#' @noRd
.count_tls_structures <- function(obj, tls_cells, cell_types, niches, Bcell_label,
                                  min_cells, max_dist, dist_multiplier, images) {
  if (is.null(images)) images <- Seurat::Images(obj)
  if (length(images) == 0) {
    stop("count_TLS = TRUE requires spatial coordinates, but the object has no images/FOVs.")
  }

  seen <- character(0)
  structure_rows <- list()
  count_rows <- list()
  cell_rows <- list()

  for (img in images) {
    all_xy <- .get_image_coords(obj, img)

    # Avoid counting the same cells twice (e.g., a cropped FOV of the same section)
    dup <- all_xy$cell %in% seen
    if (any(dup)) {
      warning(sprintf(
        "%d cells in image '%s' also appear in an earlier image and were skipped. If this image is a crop of another, set `image` to the one you want counted.",
        sum(dup), img), call. = FALSE)
      all_xy <- all_xy[!dup, , drop = FALSE]
    }
    seen <- c(seen, all_xy$cell)
    if (nrow(all_xy) < 2) next

    # Linking distance: user-supplied, or scaled from the typical cell/spot spacing
    if (is.null(max_dist)) {
      nn1 <- FNN::get.knn(as.matrix(all_xy[, c("x", "y")]), k = 1)$nn.dist[, 1]
      radius <- dist_multiplier * stats::median(nn1)
    } else {
      radius <- max_dist
    }

    xy <- all_xy[all_xy$cell %in% tls_cells, , drop = FALSE]
    n_structures <- 0L
    n_excluded <- 0L

    if (nrow(xy) > 0) {
      comp <- .connected_components(xy[, c("x", "y")], radius)
      sizes <- table(comp)
      kept <- as.integer(names(sizes)[sizes >= min_cells])
      n_excluded <- sum(sizes < min_cells)

      if (length(kept) > 0) {
        # Number structures from largest to smallest
        kept <- kept[order(sizes[as.character(kept)], decreasing = TRUE)]
        ids <- sprintf("%s_TLS_%d", img, seq_along(kept))
        id_map <- stats::setNames(ids, kept)

        in_kept <- comp %in% kept
        cell_rows[[img]] <- data.frame(
          cell = xy$cell[in_kept],
          TLS_ID = unname(id_map[as.character(comp[in_kept])]),
          stringsAsFactors = FALSE
        )

        structure_rows[[img]] <- do.call(rbind, lapply(seq_along(kept), function(i) {
          sel <- comp == kept[i]
          cts <- cell_types[xy$cell[sel]]
          n_b <- sum(cts == Bcell_label, na.rm = TRUE)
          data.frame(
            TLS_ID = ids[i],
            image = img,
            n_cells = sum(sel),
            n_Bcells = n_b,
            Bcell_Prop = n_b / sum(sel),
            centroid_x = mean(xy$x[sel]),
            centroid_y = mean(xy$y[sel]),
            niches = paste(sort(unique(as.character(niches[xy$cell[sel]]))), collapse = ", "),
            stringsAsFactors = FALSE
          )
        }))
        n_structures <- length(kept)
      }
    }

    count_rows[[img]] <- data.frame(
      image = img,
      n_TLS = n_structures,
      n_fragments_excluded = n_excluded,
      min_TLS_cells = min_cells,
      linking_distance = radius,
      stringsAsFactors = FALSE
    )
  }

  empty_summary <- data.frame(
    TLS_ID = character(0), image = character(0), n_cells = integer(0),
    n_Bcells = integer(0), Bcell_Prop = numeric(0), centroid_x = numeric(0),
    centroid_y = numeric(0), niches = character(0), stringsAsFactors = FALSE
  )
  empty_cells <- data.frame(cell = character(0), TLS_ID = character(0),
                            stringsAsFactors = FALSE)

  TLS_summary <- if (length(structure_rows)) do.call(rbind, structure_rows) else empty_summary
  TLS_cells <- if (length(cell_rows)) do.call(rbind, cell_rows) else empty_cells
  TLS_counts <- do.call(rbind, count_rows)
  rownames(TLS_summary) <- NULL
  rownames(TLS_cells) <- NULL
  rownames(TLS_counts) <- NULL

  missing <- setdiff(tls_cells, seen)
  if (length(missing) > 0) {
    warning(sprintf(
      "%d TLS-niche cells had no spatial coordinates in the selected image(s) and were not counted.",
      length(missing)), call. = FALSE)
  }

  list(TLS_summary = TLS_summary, TLS_counts = TLS_counts, TLS_cells = TLS_cells)
}

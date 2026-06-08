#' Convert a B2C Seurat object to a SpatialExperiment for hoodscanR analysis
#'
#' @param obj A Seurat object (post-aggregated B2C object) with a \code{spatial}
#'   reduction containing SPATIAL_1 and SPATIAL_2 coordinates.
#' @param annot_col The name of the metadata column in \code{obj} containing
#'   cell type annotations to pass to \code{hoodscanR::readHoodData}.
#' @return A \code{SpatialExperiment} object ready for hoodscanR analysis.
#' @export
b2c_to_hood <- function(obj, annot_col) {
  if (!requireNamespace("SpatialExperiment", quietly = TRUE))
    stop("SpatialExperiment is required but not installed. Install via: BiocManager::install('SpatialExperiment')")
  if (!requireNamespace("hoodscanR", quietly = TRUE))
    stop("hoodscanR is required but not installed. Install via: BiocManager::install('hoodscanR')")
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE))
    stop("SingleCellExperiment is required but not installed. Install via: BiocManager::install('SingleCellExperiment')")

  coords <- obj@reductions$spatial@cell.embeddings
  obj@reductions <- list()
  sce <- Seurat::as.SingleCellExperiment(obj)
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = sce@assays@data$counts),
    colData = SingleCellExperiment::colData(sce),
    spatialCoords = coords
  )
  spe <- hoodscanR::readHoodData(spe, anno_col = annot_col)
  spe
}

#' Compute a neighborhood co-localization correlation matrix using hoodscanR
#'
#' @param spe A \code{SpatialExperiment} object prepared with \code{b2c_to_hood}.
#' @param k Integer. Number of nearest neighbors to use (default 100).
#' @return A correlation matrix of cell-type neighborhood co-localization probabilities.
#' @export
get_cor_hood <- function(spe, k = 100) {
  if (!requireNamespace("hoodscanR", quietly = TRUE))
    stop("hoodscanR is required but not installed. Install via: BiocManager::install('hoodscanR')")

  fnc <- hoodscanR::findNearCells(spe, k = k)           # get neighbors
  pm  <- hoodscanR::scanHoods(fnc$distance)              # infer probability based on distance
  hoods <- hoodscanR::mergeByGroup(pm, fnc$cells)        # group by celltypes
  spe <- hoodscanR::mergeHoodSpe(spe, hoods)             # add hoods to metadata
  spe <- hoodscanR::calcMetrics(spe, pm_cols = colnames(hoods))  # generate summaries
  cor <- hoodscanR::plotColocal(spe, pm_cols = colnames(hoods), return_matrix = TRUE)  # correlation matrix
  cor
}

#' Extract the upper triangle of a square matrix as a data frame
#'
#' @param mat A square numeric matrix with row and column names.
#' @return A data frame with columns \code{celltype1}, \code{celltype2}, and \code{r}
#'   representing pairs from the upper triangle (excluding the diagonal).
#' @export
get_upper <- function(mat) {
  idx <- upper.tri(mat, diag = FALSE)
  data.frame(
    celltype1 = rownames(mat)[row(mat)[idx]],
    celltype2 = colnames(mat)[col(mat)[idx]],
    r = mat[idx]
  )
}

#' Plot hoodscanR focal cell-type neighborhood data
#'
#' Harmonizes an edge table so that the focal cell type always appears in
#' \code{celltype1}, then plots a faceted point plot.
#'
#' @param edge_df A data frame with at least columns \code{celltype1},
#'   \code{celltype2}, \code{r}, and \code{z}, plus any additional columns
#'   referenced by \code{x}, \code{metric}, and \code{group}.
#' @param focus Character. The cell type to focus on (searched in both
#'   \code{celltype1} and \code{celltype2}).
#' @param x Character. Column name to use for the x-axis.
#' @param metric Character. Column name to use for the y-axis metric.
#' @param facet_x Character. Row facet variable (default \code{"."} for no row
#'   facet).
#' @param facet_y Character. Column facet variable (default \code{"celltype2"}).
#' @param group Character. Column name used for point color/group aesthetics.
#' @param sort Logical. If \code{TRUE} (default), order facet panels by
#'   descending mean of \code{metric}.
#' @return A \code{ggplot2} object.
#' @export
plot_hood_focal <- function(edge_df, focus, x, metric, facet_x = ".", facet_y = "celltype2", group, sort = TRUE) {
  df <-
    rbind(
      edge_df[grep(focus, edge_df$celltype1),],
      edge_df[grep(focus, edge_df$celltype2),] %>%
        apply(., 1, function(u) u[c(2, 1, 3:length(u))]) %>%
        t() %>% as.data.frame() %>%
        setNames(c("celltype1", "celltype2", colnames(edge_df)[3:length(edge_df)])) %>%
        dplyr::mutate(r = as.numeric(r),
               z = as.numeric(z))
    )
  if(sort) df$celltype2 <- factor(df$celltype2, levels = (dplyr::group_by(df, celltype2) %>% dplyr::summarize(mm = mean(!!rlang::sym(metric))) %>% dplyr::arrange(dplyr::desc(mm)))$celltype2 %>% unique())
  ggplot2::ggplot(df, ggplot2::aes(x = !!rlang::sym(x), y = !!rlang::sym(metric))) +
    ggplot2::geom_hline(yintercept = 0, col = "gray") +
    ggplot2::geom_point(ggplot2::aes(group = !!rlang::sym(group), col = !!rlang::sym(group))) +
    ggplot2::facet_grid(stats::as.formula(paste(facet_x, "~", facet_y))) +
    ggplot2::theme_light() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
          strip.text.x = ggplot2::element_text(angle = 90, hjust = 0, vjust = 0.5, colour = "black"),
          strip.text.y = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5, colour = "black"),
          strip.background = ggplot2::element_rect(fill = NA),
          plot.title = ggplot2::element_text(hjust = 0.5))
}

#' Convert a hoodscanR edge table back to a symmetric matrix
#'
#' @param edge_wide A data frame (edge table) with one row per cell-type pair.
#' @param value_col Character. Column name containing the numeric values to fill
#'   the matrix (default \code{"delta_z"}).
#' @param celltype1_col Character. Column name for the first cell type
#'   (default \code{"celltype1"}).
#' @param celltype2_col Character. Column name for the second cell type
#'   (default \code{"celltype2"}).
#' @param diag_value Numeric. Value to place on the diagonal (default \code{0}).
#' @param fill_missing Numeric. Value to use for missing pairs
#'   (default \code{NA_real_}).
#' @return A symmetric named numeric matrix.
#' @export
edge_table_to_matrix <- function(edge_wide,
                                 value_col = "delta_z",
                                 celltype1_col = "celltype1",
                                 celltype2_col = "celltype2",
                                 diag_value = 0,
                                 fill_missing = NA_real_) {
  stopifnot(all(c(celltype1_col, celltype2_col, value_col) %in% colnames(edge_wide)))

  celltypes <- sort(unique(c(edge_wide[[celltype1_col]], edge_wide[[celltype2_col]])))
  mat <- matrix(fill_missing, nrow = length(celltypes), ncol = length(celltypes),
                dimnames = list(celltypes, celltypes))

  diag(mat) <- diag_value

  i <- match(edge_wide[[celltype1_col]], celltypes)
  j <- match(edge_wide[[celltype2_col]], celltypes)
  v <- edge_wide[[value_col]]

  mat[cbind(i, j)] <- v
  mat[cbind(j, i)] <- v

  mat
}

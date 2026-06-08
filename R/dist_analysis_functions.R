#' @importFrom purrr map_dfr
#' @importFrom rlang sym .data
NULL

# =====================================================================
# VisiumHD / B2CViz get_dist() neighborhood pipeline
# =====================================================================
# Designed for a directed, radius-based edge table exported by get_dist().
# Assumptions for this version:
# - one row = one origin -> neighbor pair within a search radius
# - both A->B and B->A may be present
# - every cell in a ROI appears at least once as an origin
# - no separate metadata table is required
#
# Required columns in the combined table:
#   roi, group,
#   origin_name, origin_value,
#   neighbor_name, neighbor_value,
#   distance
# =====================================================================

# ---------------------------------------------------------------------
# 1) Standardize get_dist() input
# ---------------------------------------------------------------------

#' Standardize a combined get_dist() output table
#'
#' Renames and coerces columns to a canonical format required by the
#' downstream neighborhood analysis pipeline.
#'
#' @param df A data frame combining multiple get_dist() outputs, with ROI
#'   and group columns added.
#' @param roi_col Character. Name of the ROI identifier column
#'   (default \code{"roi"}).
#' @param group_col Character. Name of the group identifier column
#'   (default \code{"group"}).
#' @param origin_id_col Character. Column for the origin cell ID
#'   (default \code{"origin_name"}).
#' @param origin_type_col Character. Column for the origin cell type / value
#'   (default \code{"origin_value"}).
#' @param neighbor_id_col Character. Column for the neighbor cell ID
#'   (default \code{"neighbor_name"}).
#' @param neighbor_type_col Character. Column for the neighbor cell type / value
#'   (default \code{"neighbor_value"}).
#' @param distance_col Character. Column for the pairwise distance
#'   (default \code{"distance"}).
#' @return A standardized data frame with columns:
#'   \code{roi}, \code{group}, \code{origin_id}, \code{origin_type},
#'   \code{neighbor_id}, \code{neighbor_type}, \code{distance}.
#'   Rows with any \code{NA} in these columns are removed.
#' @export
standardize_get_dist <- function(df,
                                 roi_col = "roi",
                                 group_col = "group",
                                 origin_id_col = "origin_name",
                                 origin_type_col = "origin_value",
                                 neighbor_id_col = "neighbor_name",
                                 neighbor_type_col = "neighbor_value",
                                 distance_col = "distance") {
  df %>%
    dplyr::transmute(
      roi = as.character(.data[[roi_col]]),
      group = as.character(.data[[group_col]]),
      origin_id = as.character(.data[[origin_id_col]]),
      origin_type = as.character(.data[[origin_type_col]]),
      neighbor_id = as.character(.data[[neighbor_id_col]]),
      neighbor_type = as.character(.data[[neighbor_type_col]]),
      distance = as.numeric(.data[[distance_col]])
    ) %>%
    dplyr::filter(
      !is.na(roi), !is.na(group),
      !is.na(origin_id), !is.na(origin_type),
      !is.na(neighbor_id), !is.na(neighbor_type),
      !is.na(distance)
    )
}

# ---------------------------------------------------------------------
# 2) Derive cell table from unique origins
# ---------------------------------------------------------------------

#' Derive a cell-level table from unique origin entries in the edge table
#'
#' Valid because all cells are assumed to occur as origins in the edge table.
#'
#' @param edges A standardized edge data frame (output of
#'   \code{\link{standardize_get_dist}}).
#' @return A data frame with columns \code{roi}, \code{group},
#'   \code{cell_id}, and \code{cell_type} — one row per unique cell.
#' @export
derive_cells_from_origins <- function(edges) {
  edges %>%
    dplyr::distinct(roi, group, cell_id = origin_id, cell_type = origin_type)
}

# ---------------------------------------------------------------------
# 3) Basic QC
# ---------------------------------------------------------------------

#' Run basic quality-control checks on a standardized edge table
#'
#' @param edges A standardized edge data frame.
#' @return A named list with QC information:
#'   \describe{
#'     \item{n_edges}{Total number of edges.}
#'     \item{n_cells}{Total number of unique origin cells.}
#'     \item{n_rois}{Number of distinct ROIs.}
#'     \item{n_groups}{Number of distinct groups.}
#'     \item{self_edges}{Rows where origin and neighbor are the same cell.}
#'     \item{duplicate_edges}{Edges that appear more than once.}
#'     \item{missing_neighbor_ids}{Neighbor IDs not found among origin IDs.}
#'     \item{origins_with_multiple_labels}{Origin cells mapped to more than one
#'       cell type.}
#'   }
#' @export
check_get_dist_table <- function(edges) {
  cells <- derive_cells_from_origins(edges)

  missing_neighbor_ids <- dplyr::anti_join(
    edges %>% dplyr::distinct(roi, neighbor_id),
    cells %>% dplyr::distinct(roi, cell_id),
    by = c("roi", "neighbor_id" = "cell_id")
  )

  mismatched_origin_labels <- edges %>%
    dplyr::distinct(roi, origin_id, origin_type) %>%
    dplyr::count(roi, origin_id) %>%
    dplyr::filter(n > 1)

  list(
    n_edges = nrow(edges),
    n_cells = nrow(cells),
    n_rois = dplyr::n_distinct(edges$roi),
    n_groups = dplyr::n_distinct(edges$group),
    self_edges = edges %>% dplyr::filter(origin_id == neighbor_id),
    duplicate_edges = edges %>% dplyr::count(roi, group, origin_id, neighbor_id, distance) %>% dplyr::filter(n > 1),
    missing_neighbor_ids = missing_neighbor_ids,
    origins_with_multiple_labels = mismatched_origin_labels
  )
}

# ---------------------------------------------------------------------
# 4) Cleaning helpers
# ---------------------------------------------------------------------

#' Remove self-edges and optionally deduplicate an edge table
#'
#' @param edges A standardized edge data frame.
#' @param drop_self Logical. If \code{TRUE} (default), remove rows where
#'   \code{origin_id == neighbor_id}.
#' @param deduplicate Logical. If \code{TRUE}, remove exact duplicate rows
#'   (default \code{FALSE}).
#' @return A cleaned edge data frame.
#' @export
clean_edges <- function(edges, drop_self = TRUE, deduplicate = FALSE) {
  out <- edges

  if (drop_self) {
    out <- out %>% dplyr::filter(origin_id != neighbor_id)
  }

  if (deduplicate) {
    out <- out %>% dplyr::distinct()
  }

  out
}

#' Subset an edge table to a distance range
#'
#' @param edges A standardized edge data frame.
#' @param max_radius Numeric or \code{NULL}. Maximum distance (inclusive).
#'   If \code{NULL} (default), no upper limit is applied.
#' @param min_radius Numeric. Minimum distance (inclusive), default \code{0}.
#' @return A filtered edge data frame.
#' @export
subset_radius <- function(edges, max_radius = NULL, min_radius = 0) {
  out <- edges %>% dplyr::filter(distance >= min_radius)
  if (!is.null(max_radius)) out <- out %>% dplyr::filter(distance <= max_radius)
  out
}

# ---------------------------------------------------------------------
# 5) ROI-level cell abundance derived from origins
# ---------------------------------------------------------------------

#' Compute cell-type abundance per ROI from origin cells
#'
#' @param edges A standardized edge data frame.
#' @return A data frame with columns \code{group}, \code{roi},
#'   \code{cell_type}, \code{n_cells}, and \code{global_freq}.
#' @export
cell_abundance <- function(edges) {
  derive_cells_from_origins(edges) %>%
    dplyr::count(group, roi, cell_type, name = "n_cells") %>%
    dplyr::group_by(group, roi) %>%
    dplyr::mutate(global_freq = n_cells / sum(n_cells)) %>%
    dplyr::ungroup()
}

# ---------------------------------------------------------------------
# 6) ROI-level pairwise neighborhood summary
# ---------------------------------------------------------------------

#' Summarize pairwise neighborhood statistics per ROI
#'
#' For each origin-type / neighbor-type pair and each ROI, compute edge
#' counts, presence counts, distance statistics, and derived enrichment
#' metrics.
#'
#' @param edges A standardized edge data frame.
#' @param max_radius Numeric or \code{NULL}. Upper distance limit applied
#'   before summarizing.
#' @param min_radius Numeric. Lower distance limit (default \code{0}).
#' @return A data frame with one row per (group, roi, origin_type,
#'   neighbor_type) combination and columns including \code{n_edges},
#'   \code{mean_neighbors_per_origin}, \code{frac_origin_with_neighbor},
#'   \code{log2_enrichment}, and others.
#' @export
summarize_pairs <- function(edges, max_radius = NULL, min_radius = 0) {
  edges_use <- subset_radius(edges, max_radius = max_radius, min_radius = min_radius)
  cells <- derive_cells_from_origins(edges)

  origin_counts <- cells %>%
    dplyr::count(group, roi, origin_type = cell_type, name = "n_origin_cells")

  neighbor_abundance <- cells %>%
    dplyr::count(group, roi, neighbor_type = cell_type, name = "n_neighbor_cells") %>%
    dplyr::group_by(group, roi) %>%
    dplyr::mutate(global_freq = n_neighbor_cells / sum(n_neighbor_cells)) %>%
    dplyr::ungroup()

  edge_counts <- edges_use %>%
    dplyr::count(group, roi, origin_type, neighbor_type, name = "n_edges")

  presence_counts <- edges_use %>%
    dplyr::distinct(group, roi, origin_id, origin_type, neighbor_type) %>%
    dplyr::count(group, roi, origin_type, neighbor_type, name = "n_origin_with_neighbor")

  distance_stats <- edges_use %>%
    dplyr::group_by(group, roi, origin_type, neighbor_type) %>%
    dplyr::summarise(
      mean_distance = mean(distance, na.rm = TRUE),
      median_distance = median(distance, na.rm = TRUE),
      q25_distance = quantile(distance, 0.25, na.rm = TRUE),
      q75_distance = quantile(distance, 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  edge_counts %>%
    dplyr::left_join(presence_counts, by = c("group", "roi", "origin_type", "neighbor_type")) %>%
    dplyr::left_join(distance_stats, by = c("group", "roi", "origin_type", "neighbor_type")) %>%
    dplyr::left_join(origin_counts, by = c("group", "roi", "origin_type")) %>%
    dplyr::left_join(neighbor_abundance, by = c("group", "roi", "neighbor_type")) %>%
    dplyr::group_by(group, roi, origin_type) %>%
    dplyr::mutate(prop_of_neighbors = n_edges / sum(n_edges)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      radius_max = ifelse(is.null(max_radius), Inf, max_radius),
      radius_min = min_radius,
      mean_neighbors_per_origin = n_edges / n_origin_cells,
      frac_origin_with_neighbor = n_origin_with_neighbor / n_origin_cells,
      enrichment_vs_global = prop_of_neighbors / global_freq,
      log2_enrichment = log2(enrichment_vs_global)
    )
}

# ---------------------------------------------------------------------
# 7) Distance-bin radial profiles
# ---------------------------------------------------------------------

#' Compute radial distance-bin profiles for cell-type pairs
#'
#' Bins edges by distance and computes per-bin neighbor counts normalized by
#' origin cell count and shell geometry.
#'
#' @param edges A standardized edge data frame.
#' @param breaks Numeric vector of distance breakpoints defining the bins
#'   (default \code{c(0, 25, 50, 100, 200, 400, 800, Inf)}).
#' @param right Logical. If \code{FALSE} (default), intervals are left-closed
#'   \code{[a, b)}.
#' @return A data frame with one row per (group, roi, origin_type,
#'   neighbor_type, dist_bin) and columns including
#'   \code{mean_neighbors_per_origin_bin},
#'   \code{neighbors_per_origin_per_area}, and
#'   \code{neighbors_per_origin_per_width}.
#' @export
summarize_radial_bins <- function(edges,
                                  breaks = c(0, 25, 50, 100, 200, 400, 800, Inf),
                                  right = FALSE) {
  cells <- derive_cells_from_origins(edges)

  origin_counts <- cells %>%
    dplyr::count(group, roi, origin_type = cell_type, name = "n_origin_cells")

  bin_tbl <- tidyr::tibble(
    dist_bin = cut(breaks[-length(breaks)],
                   breaks = breaks,
                   include.lowest = TRUE,
                   right = right),
    r_inner = breaks[-length(breaks)],
    r_outer = breaks[-1]
  ) %>%
    dplyr::mutate(
      shell_area = pi * (r_outer^2 - r_inner^2),
      shell_width = r_outer - r_inner,
      shell_midpoint = (r_inner + r_outer) / 2
    )

  edges %>%
    dplyr::mutate(
      dist_bin = cut(distance, breaks = breaks, include.lowest = TRUE, right = right)
    ) %>%
    dplyr::filter(!is.na(dist_bin)) %>%
    dplyr::count(group, roi, origin_type, neighbor_type, dist_bin, name = "n_edges") %>%
    dplyr::left_join(origin_counts, by = c("group", "roi", "origin_type")) %>%
    dplyr::left_join(bin_tbl, by = "dist_bin") %>%
    dplyr::group_by(group, roi, origin_type, dist_bin) %>%
    dplyr::mutate(
      prop_in_bin = n_edges / sum(n_edges)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      mean_neighbors_per_origin_bin = n_edges / n_origin_cells,
      neighbors_per_origin_per_area = n_edges / (n_origin_cells * shell_area),
      neighbors_per_origin_per_width = n_edges / (n_origin_cells * shell_width)
    )
}

# ---------------------------------------------------------------------
# 8) Cumulative summaries within chosen radii
# ---------------------------------------------------------------------

#' Compute cumulative neighborhood statistics at multiple radii
#'
#' For each element of \code{radii}, counts edges within that radius and
#' computes per-origin-cell metrics.
#'
#' @param edges A standardized edge data frame.
#' @param radii Numeric vector of cumulative radii to evaluate
#'   (default \code{c(25, 50, 100, 200, 400)}).
#' @return A data frame combining results for all radii, with columns
#'   \code{radius}, \code{mean_neighbors_per_origin}, and
#'   \code{frac_origin_with_neighbor}.
#' @export
summarize_cumulative <- function(edges, radii = c(25, 50, 100, 200, 400)) {
  cells <- derive_cells_from_origins(edges)

  origin_counts <- cells %>%
    dplyr::count(group, roi, origin_type = cell_type, name = "n_origin_cells")

  purrr::map_dfr(radii, function(r) {
    edges_r <- edges %>% dplyr::filter(distance <= r)

    edge_counts <- edges_r %>%
      dplyr::count(group, roi, origin_type, neighbor_type, name = "n_edges_within_r")

    presence_counts <- edges_r %>%
      dplyr::distinct(group, roi, origin_id, origin_type, neighbor_type) %>%
      dplyr::count(group, roi, origin_type, neighbor_type, name = "n_origin_with_neighbor_within_r")

    edge_counts %>%
      dplyr::full_join(presence_counts, by = c("group", "roi", "origin_type", "neighbor_type")) %>%
      tidyr::replace_na(list(
        n_edges_within_r = 0,
        n_origin_with_neighbor_within_r = 0
      )) %>%
      dplyr::left_join(origin_counts, by = c("group", "roi", "origin_type")) %>%
      dplyr::mutate(
        radius = r,
        mean_neighbors_per_origin = n_edges_within_r / n_origin_cells,
        frac_origin_with_neighbor = n_origin_with_neighbor_within_r / n_origin_cells
      )
  })
}

# ---------------------------------------------------------------------
# 9) ROI-group summaries and group comparisons
# ---------------------------------------------------------------------

#' Summarize pair statistics across ROIs within each group
#'
#' @param pair_stats A data frame produced by \code{\link{summarize_pairs}}.
#' @return A data frame with one row per (group, origin_type, neighbor_type)
#'   and columns for mean and SD of key metrics across ROIs.
#' @export
summarize_groups <- function(pair_stats) {
  pair_stats %>%
    dplyr::group_by(group, origin_type, neighbor_type) %>%
    dplyr::summarise(
      n_rois = dplyr::n_distinct(roi),
      mean_neighbors_per_origin_mean = mean(mean_neighbors_per_origin, na.rm = TRUE),
      mean_neighbors_per_origin_sd = sd(mean_neighbors_per_origin, na.rm = TRUE),
      frac_origin_with_neighbor_mean = mean(frac_origin_with_neighbor, na.rm = TRUE),
      frac_origin_with_neighbor_sd = sd(frac_origin_with_neighbor, na.rm = TRUE),
      median_distance_mean = mean(median_distance, na.rm = TRUE),
      median_distance_sd = sd(median_distance, na.rm = TRUE),
      log2_enrichment_mean = mean(log2_enrichment, na.rm = TRUE),
      log2_enrichment_sd = sd(log2_enrichment, na.rm = TRUE),
      .groups = "drop"
    )
}

#' Compare a neighborhood metric between two ROI groups
#'
#' For each origin-type / neighbor-type pair, tests for a difference between
#' \code{group_a} and \code{group_b} using a Wilcoxon rank-sum test and
#' applies Benjamini-Hochberg correction across all pairs.
#'
#' @param pair_stats A data frame produced by \code{\link{summarize_pairs}}.
#' @param metric Character. Column name of the metric to compare
#'   (default \code{"mean_neighbors_per_origin"}).
#' @param group_a Character or \code{NULL}. Label of the first group. If
#'   \code{NULL}, inferred automatically when exactly two groups are present.
#' @param group_b Character or \code{NULL}. Label of the second group.
#' @param min_rois_per_group Integer. Minimum number of ROIs with non-missing
#'   values required in each group to perform the test (default \code{2}).
#' @return A data frame with columns \code{origin_type}, \code{neighbor_type},
#'   \code{n_roi_group_a}, \code{n_roi_group_b}, \code{mean_group_a},
#'   \code{mean_group_b}, \code{delta}, \code{p_value}, and \code{p_adj},
#'   sorted by adjusted p-value and absolute delta.
#' @export
compare_groups <- function(pair_stats,
                           metric = "mean_neighbors_per_origin",
                           group_a = NULL,
                           group_b = NULL,
                           min_rois_per_group = 2) {
  if (is.null(group_a) || is.null(group_b)) {
    gs <- unique(pair_stats$group)
    if (length(gs) != 2) stop("Please provide group_a and group_b explicitly.")
    group_a <- gs[1]
    group_b <- gs[2]
  }

  pair_stats %>%
    dplyr::filter(group %in% c(group_a, group_b)) %>%
    dplyr::select(group, roi, origin_type, neighbor_type, value = dplyr::all_of(metric)) %>%
    dplyr::group_by(origin_type, neighbor_type) %>%
    dplyr::summarise(
      n_roi_group_a = sum(group == group_a & !is.na(value)),
      n_roi_group_b = sum(group == group_b & !is.na(value)),
      mean_group_a = mean(value[group == group_a], na.rm = TRUE),
      mean_group_b = mean(value[group == group_b], na.rm = TRUE),
      delta = mean_group_b - mean_group_a,
      p_value = tryCatch(
        if ((sum(group == group_a & !is.na(value)) >= min_rois_per_group) &&
            (sum(group == group_b & !is.na(value)) >= min_rois_per_group)) {
          wilcox.test(value[group == group_a], value[group == group_b])$p.value
        } else {
          NA_real_
        },
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    dplyr::arrange(p_adj, dplyr::desc(abs(delta)))
}

# ---------------------------------------------------------------------
# 10) Optional edge-aware filtering if you have cell coordinates
# ---------------------------------------------------------------------

#' Compute each cell's distance to the nearest ROI border
#'
#' @param coords A data frame with columns \code{roi}, \code{cell_id},
#'   \code{x}, and \code{y}.
#' @param roi_windows A data frame with columns \code{roi}, \code{xmin},
#'   \code{xmax}, \code{ymin}, and \code{ymax} defining the bounding box of
#'   each ROI.
#' @return \code{coords} augmented with columns \code{d_left}, \code{d_right},
#'   \code{d_bottom}, \code{d_top}, and \code{dist_to_border}.
#' @export
compute_border_distance <- function(coords, roi_windows) {
  coords %>%
    dplyr::left_join(roi_windows, by = "roi") %>%
    dplyr::mutate(
      d_left = x - xmin,
      d_right = xmax - x,
      d_bottom = y - ymin,
      d_top = ymax - y,
      dist_to_border = pmin(d_left, d_right, d_bottom, d_top)
    )
}

#' Keep only origin cells whose distance to the ROI border exceeds a radius
#'
#' Ensures that every retained origin cell has a complete neighborhood sphere
#' within the ROI, reducing edge effects.
#'
#' @param edges A standardized edge data frame.
#' @param coords A data frame with columns \code{roi}, \code{cell_id},
#'   \code{x}, and \code{y} (e.g. from the \code{locations} element of
#'   \code{\link{get_dist}} output).
#' @param roi_windows A data frame with columns \code{roi}, \code{xmin},
#'   \code{xmax}, \code{ymin}, \code{ymax}.
#' @param radius Numeric. Minimum required distance to all borders.
#' @return A filtered edge data frame containing only edges whose origin cell
#'   is at least \code{radius} away from all ROI borders.
#' @export
trim_origins_by_border <- function(edges, coords, roi_windows, radius) {
  eligible <- compute_border_distance(coords, roi_windows) %>%
    dplyr::filter(dist_to_border >= radius) %>%
    dplyr::select(roi, origin_id = cell_id)

  edges %>%
    dplyr::semi_join(eligible, by = c("roi", "origin_id"))
}

# ---------------------------------------------------------------------
# 11) Full wrapper pipeline
# ---------------------------------------------------------------------

#' Run the full neighborhood analysis pipeline
#'
#' A convenience wrapper that sequentially calls QC, optional border trimming,
#' \code{\link{summarize_pairs}}, \code{\link{summarize_radial_bins}},
#' \code{\link{summarize_cumulative}}, and \code{\link{summarize_groups}}.
#'
#' @param edges A standardized edge data frame (output of
#'   \code{\link{standardize_get_dist}}).
#' @param mean_step Numeric. The median nearest-neighbor step size (in the same
#'   units as \code{distance}) used to rescale the distance breakpoints and
#'   radii. Typically obtained as the \code{step} element of the \code{plot_b2c}
#'   return value.
#' @param distance_breaks Numeric vector. Bin breakpoints \emph{before}
#'   multiplication by \code{mean_step}
#'   (default \code{c(0, 25, 50, 100, 200, 400, 800, Inf)}).
#' @param radii Numeric vector. Cumulative radii \emph{before} multiplication
#'   by \code{mean_step} (default \code{c(25, 50, 100, 200, 400)}).
#' @param border_trim Logical. If \code{TRUE}, remove origin cells closer to
#'   the ROI border than \code{max(radii)} (requires \code{coords} and
#'   \code{roi_windows}).
#' @param coords A data frame of cell coordinates required when
#'   \code{border_trim = TRUE}.
#' @param roi_windows A data frame of ROI bounding boxes required when
#'   \code{border_trim = TRUE}.
#' @return A named list with elements: \code{qc}, \code{pair_stats},
#'   \code{radial_stats}, \code{cumulative_stats}, \code{group_pair_stats},
#'   and \code{cells}.
#' @export
run_full_pipeline <- function(edges,
                              mean_step,
                              distance_breaks = c(0, 25, 50, 100, 200, 400, 800, Inf),
                              radii = c(25, 50, 100, 200, 400),
                              border_trim = FALSE,
                              coords = NULL,
                              roi_windows = NULL) {
  distance_breaks <- distance_breaks * mean_step
  radii <- radii * mean_step

  qc <- check_get_dist_table(edges)

  edges_use <- edges
  if (border_trim) {
    if (is.null(coords) || is.null(roi_windows)) {
      stop("border_trim = TRUE requires coords and roi_windows.")
    }
    edges_use <- trim_origins_by_border(edges_use, coords = coords, roi_windows = roi_windows, radius = max(radii, na.rm = TRUE))
  }

  pair_stats      <- summarize_pairs(edges_use)
  radial_stats    <- summarize_radial_bins(edges_use, breaks = distance_breaks)
  cumulative_stats <- summarize_cumulative(edges_use, radii = radii)
  group_pair_stats <- summarize_groups(pair_stats)

  list(
    qc               = qc,
    pair_stats       = pair_stats,
    radial_stats     = radial_stats,
    cumulative_stats = cumulative_stats,
    group_pair_stats = group_pair_stats,
    cells            = derive_cells_from_origins(edges_use)
  )
}

# ---------------------------------------------------------------------
# 12) Scaling back to microns
# ---------------------------------------------------------------------

#' Rescale distance columns in a pipeline result object to physical units
#'
#' After running \code{\link{run_full_pipeline}}, use this function to divide
#' all distance-valued columns by \code{factor} (converting from internal
#' coordinate units to microns, for example).
#'
#' @param res A list returned by \code{\link{run_full_pipeline}}.
#' @param factor Positive numeric. The number of coordinate units per micron
#'   (default \code{10}).
#' @return The same list structure with distance columns rescaled and
#'   \code{dist_bin} labels updated. An attribute
#'   \code{"distance_rescale_factor"} is attached to the returned list.
#' @export
rescale_res_distances <- function(res, factor = 10) {
  if (is.null(factor) || !is.numeric(factor) || length(factor) != 1 || factor <= 0) {
    stop("factor must be one positive numeric value.")
  }

  out <- res

  scale_if_present <- function(df, cols, factor) {
    if (is.null(df)) return(df)
    cols_present <- intersect(cols, colnames(df))
    if (length(cols_present) > 0) {
      df <- df %>%
        dplyr::mutate(dplyr::across(dplyr::all_of(cols_present), ~ .x / factor))
    }
    df
  }

  relabel_dist_bin <- function(x, factor) {
    x_chr <- as.character(x)
    x_chr <- gsub("\\[|\\]|\\(|\\)", "", x_chr)
    parts <- strsplit(x_chr, ",", fixed = TRUE)

    out_chr <- vapply(parts, function(p) {
      if (length(p) != 2) return(NA_character_)
      left <- trimws(p[1])
      right <- trimws(p[2])

      left_num <- suppressWarnings(as.numeric(left))
      right_num <- suppressWarnings(as.numeric(right))

      left_new <- if (is.na(left_num)) left else format(round(left_num / factor), trim = TRUE, scientific = FALSE)
      right_new <- if (is.na(right_num)) right else format(round(right_num / factor), trim = TRUE, scientific = FALSE)

      paste0("[", left_new, ", ", right_new, ")")
    }, character(1))

    factor(out_chr, levels = unique(out_chr))
  }

  if (!is.null(out$pair_stats)) {
    out$pair_stats <- scale_if_present(
      out$pair_stats,
      cols = c("mean_distance", "median_distance", "q25_distance", "q75_distance",
               "radius_min", "radius_max"),
      factor = factor
    )
  }

  if (!is.null(out$radial_stats)) {
    out$radial_stats <- scale_if_present(
      out$radial_stats,
      cols = c("r_inner", "r_outer", "shell_width", "shell_midpoint"),
      factor = factor
    )

    if ("shell_area" %in% colnames(out$radial_stats)) {
      out$radial_stats <- out$radial_stats %>%
        dplyr::mutate(shell_area = shell_area / (factor^2))
    }

    if ("neighbors_per_origin_per_area" %in% colnames(out$radial_stats)) {
      out$radial_stats <- out$radial_stats %>%
        dplyr::mutate(neighbors_per_origin_per_area = neighbors_per_origin_per_area * (factor^2))
    }

    if ("neighbors_per_origin_per_width" %in% colnames(out$radial_stats)) {
      out$radial_stats <- out$radial_stats %>%
        dplyr::mutate(neighbors_per_origin_per_width = neighbors_per_origin_per_width * factor)
    }

    if ("dist_bin" %in% colnames(out$radial_stats)) {
      out$radial_stats$dist_bin <- relabel_dist_bin(out$radial_stats$dist_bin, factor = factor)
    }
  }

  if (!is.null(out$cumulative_stats)) {
    out$cumulative_stats <- scale_if_present(
      out$cumulative_stats,
      cols = c("radius"),
      factor = factor
    )
  }

  if (!is.null(out$group_pair_stats)) {
    out$group_pair_stats <- scale_if_present(
      out$group_pair_stats,
      cols = c("median_distance_mean", "median_distance_sd"),
      factor = factor
    )
  }

  attr(out, "distance_rescale_factor") <- factor
  out
}

# ---------------------------------------------------------------------
# 13) Plotting functions
# ---------------------------------------------------------------------

#' Plot neighborhood statistics from a pipeline result
#'
#' Creates a faceted point (and optionally line) plot for a chosen focal
#' origin cell type, with flexible axis and facet mappings.
#'
#' @param res A list returned by \code{\link{run_full_pipeline}}.
#' @param table Character. Name of the result list element to plot (e.g.
#'   \code{"pair_stats"}, \code{"radial_stats"}, \code{"cumulative_stats"}).
#' @param focus Character. The origin cell type to focus on.
#' @param x Character. Column name for the x-axis.
#' @param y Character. Column name for the y-axis metric.
#' @param facet_x Character. Row facet formula term (default \code{"."}).
#' @param facet_y Character. Column facet formula term (default \code{"."}).
#' @param group Character. Column name used for color / group aesthetics.
#' @param plot.line Logical. If \code{TRUE} (default), add connecting lines
#'   between points grouped by \code{id}.
#' @return A \code{ggplot2} object.
#' @export
plot_stats <- function(res, table, focus, x, y, facet_x = ".", facet_y = ".", group, plot.line = TRUE) {
  df <- res[[table]]
  df$id <- paste0(df$neighbor_type, "_", df$roi)
  definitions <-
    rbind(
      c("n_edges", "Number of directed origin-to-neighbor pairs"),
      c("n_origin_with_neighbor", "Number of distinct origin cells with neighbor"),
      c("mean_distance", "Mean distance between origin and neighbor"),
      c("median_distance", "Median distance between origin and neighbor"),
      c("q25_distance", "Q25 distance between origin and neighbor"),
      c("q75_distance", "Q75 distance between origin and neighbor"),
      c("n_origin_cells", "Total number of origin cells"),
      c("n_neighbor_cells", "Total number of neighbor cells"),
      c("global_freq", "Fraction of neighbor in the ROI globally"),
      c("prop_of_neighbors", "Fraction of origin-neighbor pair vs all origin pairs"),
      c("mean_neighbors_per_origin", "Average number of neighbors around one origin cell"),
      c("frac_origin_with_neighbor", "Fraction of origin cells paired with this neighbor"),
      c("enrichment_vs_global", "Neighbor enrichment vs ROI-wide abundance"),
      c("log2_enrichment", "Neighbor enrichment vs ROI-wide abundance"),
      c("dist_bin", "Distance bin"),
      c("prop_in_bin", "Fraction of origin-neighbor pair vs all origin pairs"),
      c("mean_neighbors_per_origin_bin", "Average number of neighbors per bin"),
      c("neighbors_per_origin_per_area", "Average number of neighbors per bin area"),
      c("neighbors_per_origin_per_width", "Average number of neighbors per bin width"),
      c("n_edges_within_r", "Number of directed origin-to-neighbor pairs"),
      c("n_origin_with_neighbor_within_r", "Number of distinct origin cells with neighbor")
    )
  matched.def <- match(y, definitions[,1])
  def <- ifelse(is.na(matched.def), y, definitions[matched.def, 2])
  p <-
    ggplot2::ggplot(dplyr::filter(df, origin_type == focus), ggplot2::aes(!!rlang::sym(x), !!rlang::sym(y))) +
    ggplot2::geom_hline(yintercept = 0, col = "gray") +
    ggplot2::geom_point(ggplot2::aes(group = id, col = !!rlang::sym(group))) +
    ggplot2::facet_grid(stats::as.formula(paste(facet_x, "~", facet_y))) +
    ggplot2::theme_light() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
          strip.text.x = ggplot2::element_text(angle = 90, hjust = 0, vjust = 0.5, colour = "black"),
          strip.text.y = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5, colour = "black"),
          strip.background = ggplot2::element_rect(fill = NA),
          plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::ggtitle(paste0(def, " (in ", table, " data)"))
  if(plot.line) p <- p + ggplot2::geom_line(ggplot2::aes(group = id, col = !!rlang::sym(group)))
  p
}

#' Plot group-level heatmap for a focal cell type
#'
#' Requires the \code{circlize} and \code{ComplexHeatmap} packages
#' (Bioconductor). Creates a heatmap where rows are groups and columns are
#' neighbor cell types.
#'
#' @param res A list returned by \code{\link{run_full_pipeline}}.
#' @param focus Character. The origin cell type to focus on.
#' @param metric Character. The column in \code{res$group_pair_stats} to
#'   visualize (e.g. \code{"mean_neighbors_per_origin_mean"}).
#' @param clustered Logical. If \code{TRUE}, cluster columns by hierarchical
#'   clustering (columns without any \code{NA} values are required).
#'   Default \code{FALSE}.
#' @return A \code{ComplexHeatmap::Heatmap} object.
#' @export
plot_hm <- function(res, focus, metric, clustered = FALSE) {
  if (!requireNamespace("circlize", quietly = TRUE))
    stop("circlize is required but not installed. Install via: BiocManager::install('circlize')")
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE))
    stop("ComplexHeatmap is required but not installed. Install via: BiocManager::install('ComplexHeatmap')")

  df <- dplyr::filter(res$group_pair_stats, origin_type == focus) %>%
    dplyr::select(group, neighbor_type, !!rlang::sym(metric))

  mat <- t(sapply(split(df, df$group), function(u) {
    u[[metric]][match(unique(sort(df$neighbor_type)), u$neighbor_type)]
  }))
  colnames(mat) <- unique(sort(df$neighbor_type))

  if(clustered) mat <- mat[, apply(mat, 2, function(u) !any(is.na(u)))]

  mat.color <- NULL
  if(min(mat, na.rm = TRUE) < 0) {
    mat.color <- circlize::colorRamp2(c(-max(abs(mat), na.rm = TRUE), 0, max(abs(mat), na.rm = TRUE)), c("lightblue", "white", "orangered"))
  } else {
    mat.color <- circlize::colorRamp2(c(0, max(mat, na.rm = TRUE)/2, max(mat, na.rm = TRUE)), c("#ffffcc", "gold", "#800026"))
  }

  ComplexHeatmap::Heatmap(mat,
          column_title = metric,
          name = metric,
          cluster_rows = FALSE,
          cluster_columns = clustered,
          na_col = "gray30",
          row_names_side = "left",
          col = mat.color)
}

#' Plot a 2-D neighbor density map centered on a focal cell type
#'
#' Samples \code{n} origin cells of the specified type, centers their
#' neighbor coordinates at the origin, and visualizes the resulting spatial
#' density using \code{geom_bin_2d}.
#'
#' @param output_data A named list where each element corresponds to a ROI
#'   and contains, at index 1, a distances data frame and, at index 2, a
#'   locations data frame (matching the format returned by
#'   \code{\link{get_dist}}).
#' @param roi Character. Name of the ROI to use (must be a key in
#'   \code{output_data}).
#' @param origin Character. The origin cell type (matched against
#'   \code{origin_value}).
#' @param neighbor Character. The neighbor cell type (matched against
#'   \code{neighbor_value}).
#' @param n Integer. Number of origin cells to sample (default \code{100}).
#' @param seed Integer. Random seed for sampling (default \code{123}).
#' @param ... Additional arguments passed to
#'   \code{ggplot2::scale_fill_gradientn}.
#' @return A \code{ggplot2} object showing the aggregated neighbor density.
#' @export
plot_density <- function(output_data, roi, origin, neighbor, n = 100, seed = 123, ...) {
  sub.cells <- dplyr::filter(output_data[[roi]][[1]],
                    origin_value == origin,
                    neighbor_value == neighbor)
  set.seed(seed)
  origins <- sample(sub.cells$origin_name, n)
  cell.data <- list()
  for(elmt in origins) {
    selected.cells <- dplyr::filter(output_data[[roi]][[1]], origin_name == elmt, neighbor_value == neighbor)
    cells.location <- dplyr::filter(output_data[[roi]][[2]],
                             cell_id %in% c(selected.cells$neighbor_name, selected.cells$origin_name))
    # center data
    cells.location$x <- cells.location$x - dplyr::filter(cells.location, cell_id == elmt)$x
    cells.location$y <- cells.location$y - dplyr::filter(cells.location, cell_id == elmt)$y
    cell.data[[elmt]] <- cells.location
  }
  final.data <- do.call("rbind", cell.data)
  final.data <- dplyr::filter(final.data, !cell_id %in% origins)
  ggplot2::ggplot(data.frame()) +
    ggplot2::geom_bin_2d(data = final.data, ggplot2::aes(x, y), size = 1) +
    ggplot2::scale_fill_gradientn(colors = c("#4575B4", "#ABD9E9", "#FEE090", "#F46D43", "#A50026"), ...) +
    ggplot2::coord_fixed() +
    ggplot2::theme_light() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          axis.text  = ggplot2::element_blank(),
          axis.title = ggplot2::element_blank())
}

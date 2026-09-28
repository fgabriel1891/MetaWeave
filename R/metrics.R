#' Simulate an interaction-count matrix
#' @param network Ecological network or probability matrix.
#' @param size Total number of interactions.
#' @param seed Optional random seed.
#' @return An integer matrix with the input dimensions and names. Counts sum to `size` when positive probability mass exists; otherwise all entries are zero.
#' @export
simulate_interaction_counts <- function(network, size = 100L, seed = NULL) {
  mat <- if (inherits(network, "ecological_network")) network$matrix else network
  if (!is.null(seed)) set.seed(seed)
  probs <- as.numeric(mat)
  if (!length(probs) || sum(probs, na.rm = TRUE) <= 0) return(matrix(0L, nrow(mat), ncol(mat), dimnames = dimnames(mat)))
  probs[is.na(probs)] <- 0
  out <- matrix(as.integer(stats::rmultinom(1, size, probs / sum(probs))), nrow(mat), ncol(mat), dimnames = dimnames(mat))
  out
}

#' Simulate a network
#' @param network Probability network.
#' @param size Total interaction count.
#' @param seed Optional seed.
#' @return A weighted `ecological_network` containing simulated interaction counts and the source network in its metadata.
#' @export
simulate_network <- function(network, size = 100L, seed = NULL) {
  new_ecological_network(simulate_interaction_counts(network, size, seed), "weighted", list(source = network),
                         network$row_group, network$column_group,
                         directed = isTRUE(network$directed),
                         self_links = if (is.null(network$self_links)) NA else network$self_links)
}

#' Summarize one ecological network
#' @param network Ecological network or numeric matrix.
#' @param threshold Optional threshold for realized links.
#' @return A one-row data frame with `row_richness`, `column_richness`,
#'   `possible_links` (number of eligible matrix entries), `expected_links`
#'   (sum of eligible entries, or number at least `threshold`), and
#'   `mean_probability` (mean of eligible nonmissing entries).
#'   For networks with `self_links = FALSE`, pairs with the same species name
#'   are excluded. Unnamed same-group square networks use the diagonal;
#'   unnamed networks across different groups require species names.
#'   Other inputs use all matrix entries, including missing entries in the
#'   possible-link count. Missing values are omitted from sums and means.
#'   Without eligible entries the mean is `NA`; if all are missing it is `NaN`.
#'   For weighted/observed networks these fields summarize weights, not
#'   link probabilities; use `threshold` to count entries at or above a cutoff.
#' @examples
#' p <- matrix(c(0.8, 0.3), 2)
#' summarize_network(p)
#' summarize_network(p, threshold = 0.5)
#' @export
summarize_network <- function(network, threshold = NULL) {
  mat <- if (inherits(network, "ecological_network")) network$matrix else network
  allowed <- matrix(TRUE, nrow(mat), ncol(mat))
  if (inherits(network, "ecological_network") && identical(network$self_links, FALSE)) {
    if (!is.null(rownames(mat)) && !is.null(colnames(mat))) {
      allowed <- outer(rownames(mat), colnames(mat), `!=`)
    } else if (identical(network$row_group, network$column_group)) {
      if (nrow(mat) != ncol(mat)) stop("Unnamed same-group networks must be square.")
      diag(allowed) <- FALSE
    } else {
      stop("Species names are required to identify forbidden self-links across groups.")
    }
  }
  entries <- mat[allowed]
  possible <- length(entries)
  expected <- if (is.null(threshold)) sum(entries, na.rm = TRUE) else sum(entries >= threshold, na.rm = TRUE)
  data.frame(row_richness = nrow(mat), column_richness = ncol(mat), possible_links = possible,
             expected_links = expected, mean_probability = if (possible) mean(entries, na.rm = TRUE) else NA_real_)
}

#' Summarize a collection of networks
#' @param x Network collection.
#' @param threshold Optional probability threshold.
#' @return The input `network_collection` with summary columns appended to its `index`; see [summarize_network()] for column definitions.
#' @export
summarize_networks <- function(x, threshold = NULL) {
  metrics <- do.call(rbind, lapply(x$networks, summarize_network, threshold = threshold))
  x$index <- cbind(x$index, metrics)
  x
}

#' Placeholder adapter for rarefaction packages
#' @description Simulates a standardized count matrix and evaluates a user-supplied metric function.
#' @param prob_mat Probability matrix.
#' @param target_size Standardized number of interactions.
#' @param metric A function accepting a count matrix.
#' @param n_per_level Number of replicates.
#' @param seed Optional seed.
#' @return A numeric mean of the supplied metric across simulated count matrices, with missing metric values omitted.
#' @export
rarefied_network_metric <- function(prob_mat, target_size = 100L, metric, n_per_level = 10L, seed = NULL) {
  if (!is.function(metric)) stop("`metric` must be a function accepting a count matrix.")
  seeds <- if (is.null(seed)) rep(NA_integer_, n_per_level) else seed + seq_len(n_per_level) - 1L
  vals <- vapply(seeds, function(s) metric(simulate_interaction_counts(prob_mat, target_size, if (is.na(s)) NULL else s)), numeric(1))
  mean(vals, na.rm = TRUE)
}

#' Add standard metrics to network collection
#' @inheritParams summarize_networks
#' @return A `network_collection` with summary columns appended to `index`; see [summarize_networks()].
#' @export
add_network_metrics <- function(x, threshold = NULL) summarize_networks(x, threshold)

#' Map cell metrics back to rasters
#' @param x Network collection with metrics in its index.
#' @param metrics Metric columns; inferred when omitted.
#' @return A `spatial_result` with a multilayer `terra::SpatRaster` in `data`, one layer per requested metric, and metric names in `metrics`.
#' @export
map_metrics <- function(x, metrics = NULL) {
  if (is.null(x$template)) stop("A raster template is required.")
  reserved <- c("cell_id", "x", "y")
  if (is.null(metrics)) metrics <- setdiff(names(x$index)[vapply(x$index, is.numeric, logical(1))], reserved)
  layers <- lapply(metrics, function(nm) {
    r <- terra::rast(x$template)
    terra::values(r) <- NA_real_
    r[x$index$cell_id] <- x$index[[nm]]
    names(r) <- nm
    r
  })
  new_spatial_result(terra::rast(layers), metrics)
}

#' Compatibility alias returning metric rasters
#' @inheritParams map_metrics
#' @return A `terra::SpatRaster` with one layer per requested metric.
#' @export
network_metrics_raster <- function(x, metrics = NULL) map_metrics(x, metrics)$data

#' Complete compatibility workflow
#' @param ... Arguments passed to `downscale_networks()`.
#' @return A list with `result` (a summarized `network_collection`) and `rasters` (a `terra::SpatRaster` of metrics).
#' @export
run_network_downscaling <- function(...) {
  result <- add_network_metrics(downscale_networks(...))
  list(result = result, rasters = network_metrics_raster(result))
}

#' Run generic spatial ecological-network inference
#' @param distributions Named distribution groups.
#' @param model Inference model.
#' @param min_species Minimum local richness in each group.
#' @param threshold Optional probability threshold for network summaries.
#' @return A list with `result` (a summarized `network_collection`) and `spatial` (a `spatial_result` containing metric rasters).
#' @export
run_spatial_inference <- function(distributions, model, min_species = 1L, threshold = NULL) {
  result <- summarize_networks(
    reconstruct_networks(distributions, model, min_species),
    threshold = threshold
  )
  list(result = result, spatial = map_metrics(result))
}

#' Convert an interaction table to a matrix
#' @param data Data frame.
#' @param row_col,column_col,value_col Column names.
#' @param fill Missing-pair value.
#' @return A numeric matrix with row and column species names. Unspecified pairs receive `fill`; repeated pairs use the last supplied value.
#' @examples
#' records <- data.frame(plant = c("a", "b"), animal = c("x", "x"), p = c(0.8, 0.3))
#' interaction_table_to_matrix(records, "plant", "animal", "p")
#' @export
interaction_table_to_matrix <- function(data, row_col, column_col, value_col, fill = 0) {
  rows <- unique(as.character(data[[row_col]])); cols <- unique(as.character(data[[column_col]]))
  out <- matrix(fill, length(rows), length(cols), dimnames = list(rows, cols))
  out[cbind(match(data[[row_col]], rows), match(data[[column_col]], cols))] <- data[[value_col]]
  out
}

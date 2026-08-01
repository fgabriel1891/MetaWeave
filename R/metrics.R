#' Simulate an interaction-count matrix
#' @param network Ecological network or probability matrix.
#' @param size Total number of interactions.
#' @param seed Optional random seed.
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
#' @export
simulate_network <- function(network, size = 100L, seed = NULL) {
  new_ecological_network(simulate_interaction_counts(network, size, seed), "weighted", list(source = network),
                         network$row_group, network$column_group)
}

#' Summarize one ecological network
#' @param network Ecological network or numeric matrix.
#' @param threshold Optional threshold for realized links.
#' @export
summarize_network <- function(network, threshold = NULL) {
  mat <- if (inherits(network, "ecological_network")) network$matrix else network
  possible <- length(mat)
  expected <- if (is.null(threshold)) sum(mat, na.rm = TRUE) else sum(mat >= threshold, na.rm = TRUE)
  data.frame(row_richness = nrow(mat), column_richness = ncol(mat), possible_links = possible,
             expected_links = expected, mean_probability = if (possible) mean(mat, na.rm = TRUE) else NA_real_)
}

#' Summarize a collection of networks
#' @param x Network collection.
#' @param threshold Optional probability threshold.
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
#' @export
rarefied_network_metric <- function(prob_mat, target_size = 100L, metric, n_per_level = 10L, seed = NULL) {
  if (!is.function(metric)) stop("`metric` must be a function accepting a count matrix.")
  seeds <- if (is.null(seed)) rep(NA_integer_, n_per_level) else seed + seq_len(n_per_level) - 1L
  vals <- vapply(seeds, function(s) metric(simulate_interaction_counts(prob_mat, target_size, if (is.na(s)) NULL else s)), numeric(1))
  mean(vals, na.rm = TRUE)
}

#' Add standard metrics to network collection
#' @export
add_network_metrics <- function(x, threshold = NULL) summarize_networks(x, threshold)

#' Map cell metrics back to rasters
#' @param x Network collection with metrics in its index.
#' @param metrics Metric columns; inferred when omitted.
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
#' @export
network_metrics_raster <- function(x, metrics = NULL) map_metrics(x, metrics)$data

#' Complete compatibility workflow
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
#' @export
interaction_table_to_matrix <- function(data, row_col, column_col, value_col, fill = 0) {
  rows <- unique(as.character(data[[row_col]])); cols <- unique(as.character(data[[column_col]]))
  out <- matrix(fill, length(rows), length(cols), dimnames = list(rows, cols))
  out[cbind(match(data[[row_col]], rows), match(data[[column_col]], cols))] <- data[[value_col]]
  out
}

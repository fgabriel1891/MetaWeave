#' Construct a distribution collection
#'
#' @param data A named `terra::SpatRaster` whose layers represent species.
#' @param guild Optional guild label for the collection.
#' @param metadata Optional species metadata data frame.
#' @export
new_distribution_collection <- function(data, guild = NULL, metadata = NULL) {
  if (!inherits(data, "SpatRaster")) stop("`data` must be a SpatRaster.")
  if (is.null(names(data)) || any(!nzchar(names(data)))) stop("Raster layers must be named.")
  structure(list(data = data, guild = guild, metadata = metadata), class = "distribution_collection")
}

#' @export
print.distribution_collection <- function(x, ...) {
  cat("<distribution_collection>\n", length(x$data), " species\n", sep = "")
  invisible(x)
}

#' Construct a local assemblage
#' @param species Named list of character vectors, one per ecological group.
#' @param cell_id Optional raster cell identifier.
#' @param coordinates Optional named numeric vector containing x and y.
#' @export
new_assemblage <- function(species, cell_id = NA_integer_, coordinates = NULL) {
  if (!is.list(species) || is.null(names(species))) stop("`species` must be a named list.")
  species <- lapply(species, as.character)
  structure(list(species = species, cell_id = cell_id, coordinates = coordinates), class = "assemblage")
}

#' @export
print.assemblage <- function(x, ...) {
  cat("<assemblage>", if (!is.na(x$cell_id)) paste(" cell", x$cell_id), "\n")
  for (nm in names(x$species)) cat(nm, ": ", length(x$species[[nm]]), " species\n", sep = "")
  invisible(x)
}

#' Construct a pluggable inference model
#' @param predict A function accepting `(model, assemblage)` and returning an ecological network.
#' @param name Model name.
#' @param parameters Model-specific parameters.
#' @param subclass Optional additional class name.
#' @export
new_inference_model <- function(predict, name = "custom", parameters = list(), subclass = NULL) {
  if (!is.function(predict)) stop("`predict` must be a function.")
  structure(list(name = name, predict = predict, parameters = parameters),
            class = c(subclass, "inference_model"))
}

#' @export
print.inference_model <- function(x, ...) {
  cat("<inference_model>", x$name, "\n")
  invisible(x)
}

#' Construct an ecological network
#' @param matrix Numeric interaction matrix.
#' @param type One of probability, binary, weighted, or observed.
#' @param metadata Optional metadata.
#' @param row_group,column_group Ecological groups represented by rows and columns.
#' @param directed Whether matrix orientation represents directed interactions.
#' @param self_links Whether self-links are permitted; `NA` when not applicable.
#' @export
new_ecological_network <- function(matrix, type = c("probability", "binary", "weighted", "observed"), metadata = list(),
                                   row_group = "rows", column_group = "columns",
                                   directed = FALSE, self_links = NA) {
  type <- match.arg(type)
  if (!is.matrix(matrix) || !is.numeric(matrix)) stop("`matrix` must be a numeric matrix.")
  if (type == "probability" && any(matrix < 0 | matrix > 1, na.rm = TRUE)) stop("Probabilities must lie in [0, 1].")
  if (!is.logical(directed) || length(directed) != 1L || is.na(directed)) stop("`directed` must be TRUE or FALSE.")
  if (!is.logical(self_links) || length(self_links) != 1L) stop("`self_links` must be TRUE, FALSE, or NA.")
  structure(list(matrix = matrix, type = type, metadata = metadata,
                 row_group = as.character(row_group), column_group = as.character(column_group),
                 directed = directed, self_links = self_links),
            class = c(paste0(type, "_network"), "ecological_network"))
}

#' @export
print.ecological_network <- function(x, ...) {
  cat("<ecological_network:", x$type, ">", nrow(x$matrix), x$row_group,
      "x", ncol(x$matrix), x$column_group,
      if (isTRUE(x$directed)) "directed" else "undirected", "\n")
  invisible(x)
}

#' Construct a collection of spatial ecological networks
#' @param networks List of ecological networks.
#' @param index Cell-level data frame.
#' @param template Spatial raster template.
#' @param settings Optional workflow settings.
#' @export
new_network_collection <- function(networks, index, template = NULL, settings = list()) {
  if (!is.list(networks) || !is.data.frame(index) || length(networks) != nrow(index)) {
    stop("`networks` and rows of `index` must have equal length.")
  }
  structure(list(networks = networks, index = index, template = template, settings = settings),
            class = "network_collection")
}

#' @export
print.network_collection <- function(x, ...) {
  cat("<network_collection>", length(x$networks), " local networks\n")
  invisible(x)
}

#' Construct a spatial result
#' @param data A spatial object.
#' @param metrics Names of mapped metrics.
#' @export
new_spatial_result <- function(data, metrics = names(data)) {
  structure(list(data = data, metrics = metrics), class = "spatial_result")
}

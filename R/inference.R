#' Infer a network for one assemblage
#' @param assemblage An assemblage.
#' @param model An inference model.
#' @export
infer_network <- function(assemblage, model) {
  if (!inherits(assemblage, "assemblage")) stop("`assemblage` must be an assemblage.")
  if (!inherits(model, "inference_model")) stop("`model` must be an inference model.")
  out <- model$predict(model, assemblage)
  if (!inherits(out, "ecological_network")) stop("Model must return an ecological network.")
  out
}

#' Infer networks for many assemblages
#' @param assemblages List of assemblages.
#' @param model Inference model.
#' @param template Optional spatial template.
#' @export
infer_networks <- function(assemblages, model, template = NULL) {
  nets <- lapply(assemblages, infer_network, model = model)
  idx <- do.call(rbind, lapply(assemblages, function(a) {
    coordinates <- as.numeric(a$coordinates)
    x <- if (length(coordinates) >= 1L) coordinates[[1L]] else NA_real_
    y <- if (length(coordinates) >= 2L) coordinates[[2L]] else NA_real_
    data.frame(
      cell_id = unname(a$cell_id),
      x = unname(x),
      y = unname(y),
      row.names = NULL
    )
  }))
  rownames(idx) <- NULL
  new_network_collection(nets, idx, template, settings = list(model = model$name))
}

#' Reconstruct spatial networks from arbitrary distribution groups
#' @param distributions Named list of aligned distribution rasters or collections.
#' @param model Inference model whose group names match `distributions`.
#' @param min_species Minimum number of locally present species required in every group.
#' @export
reconstruct_networks <- function(distributions, model, min_species = 1L) {
  if (!is.list(distributions) || !length(distributions)) stop("`distributions` must be a non-empty named list.")
  template <- distributions[[1L]]
  if (inherits(template, "distribution_collection")) template <- template$data
  assemblages <- assemble_communities(distributions, min_species = min_species)
  infer_networks(assemblages, model, template = template[[1L]])
}

#' Construct a species-level probability-matrix model
#'
#' Use this backend when an external model returns probabilities for named
#' species pairs directly, including ensemble-averaged predictions for which
#' no single coherent block partition exists.
#'
#' @param probability_matrix Numeric species-by-species probability matrix.
#'   Row and column names are required and must be unique.
#' @param row_group Assemblage group corresponding to matrix rows.
#' @param column_group Assemblage group corresponding to matrix columns.
#' @param name Descriptive model name stored in outputs.
#' @export
probability_matrix_model <- function(probability_matrix, row_group, column_group,
                                     name = "species probability matrix") {
  if (!is.matrix(probability_matrix) || !is.numeric(probability_matrix)) {
    stop("`probability_matrix` must be a numeric matrix.")
  }
  if (is.null(rownames(probability_matrix)) || is.null(colnames(probability_matrix))) {
    stop("`probability_matrix` requires species row and column names.")
  }
  if (any(!nzchar(rownames(probability_matrix))) || any(!nzchar(colnames(probability_matrix)))) {
    stop("Probability-matrix species names cannot be empty.")
  }
  if (anyDuplicated(rownames(probability_matrix)) || anyDuplicated(colnames(probability_matrix))) {
    stop("Probability-matrix species names must be unique within each group.")
  }
  if (any(!is.finite(probability_matrix))) {
    stop("`probability_matrix` cannot contain missing or infinite values.")
  }
  if (any(probability_matrix < 0 | probability_matrix > 1)) {
    stop("Species-pair probabilities must lie in [0, 1].")
  }
  if (!is.character(row_group) || length(row_group) != 1L || !nzchar(row_group) ||
      !is.character(column_group) || length(column_group) != 1L || !nzchar(column_group)) {
    stop("`row_group` and `column_group` must each be one non-empty name.")
  }
  if (identical(row_group, column_group)) stop("Bipartite group names must differ.")

  pred <- function(model, assemblage) {
    par <- model$parameters
    if (!all(c(par$row_group, par$column_group) %in% names(assemblage$species))) {
      stop("Assemblage lacks model groups: ", par$row_group, " and ", par$column_group, ".")
    }
    rows <- assemblage$species[[par$row_group]]
    columns <- assemblage$species[[par$column_group]]
    rows <- rows[rows %in% rownames(par$probability_matrix)]
    columns <- columns[columns %in% colnames(par$probability_matrix)]
    mat <- par$probability_matrix[rows, columns, drop = FALSE]
    new_ecological_network(
      mat,
      type = "probability",
      metadata = list(cell_id = assemblage$cell_id, model = model$name),
      row_group = par$row_group,
      column_group = par$column_group
    )
  }

  new_inference_model(
    pred,
    name = name,
    parameters = list(
      probability_matrix = probability_matrix,
      row_group = row_group,
      column_group = column_group
    ),
    subclass = "probability_matrix_model"
  )
}

#' Construct a generic stochastic block model backend
#' @param row_lookup Data frame mapping species in the matrix rows to blocks.
#' @param column_lookup Data frame mapping species in the matrix columns to blocks.
#' @param theta Guild-by-guild probability matrix.
#' @param row_group Name of the assemblage group placed in matrix rows.
#' @param column_group Name of the assemblage group placed in matrix columns.
#' @param species_col Species column in lookup tables.
#' @param guild_col Guild column in lookup tables.
#' @export
block_model <- function(row_lookup, column_lookup, theta, row_group, column_group,
                        species_col = "species", guild_col = "guild") {
  validate_lookup <- function(x, label) {
    if (!all(c(species_col, guild_col) %in% names(x))) stop(label, " lookup lacks required columns.")
    if (anyDuplicated(x[[species_col]])) stop(label, " species must be unique.")
  }
  validate_lookup(row_lookup, "Row-group")
  validate_lookup(column_lookup, "Column-group")
  if (!is.character(row_group) || length(row_group) != 1L || !nzchar(row_group)) stop("`row_group` must be one non-empty name.")
  if (!is.character(column_group) || length(column_group) != 1L || !nzchar(column_group)) stop("`column_group` must be one non-empty name.")
  if (identical(row_group, column_group)) stop("Bipartite row and column groups must have different names.")
  if (!is.matrix(theta) || any(theta < 0 | theta > 1, na.rm = TRUE)) stop("`theta` must be a probability matrix.")
  pred <- function(model, assemblage) {
    par <- model$parameters
    if (!all(c(par$row_group, par$column_group) %in% names(assemblage$species))) {
      stop("Assemblage lacks model groups: ", par$row_group, " and ", par$column_group, ".")
    }
    rows <- assemblage$species[[par$row_group]]
    columns <- assemblage$species[[par$column_group]]
    ri <- match(rows, par$row_lookup[[par$species_col]])
    ci <- match(columns, par$column_lookup[[par$species_col]])
    keep_r <- !is.na(ri); keep_c <- !is.na(ci)
    rows <- rows[keep_r]; columns <- columns[keep_c]; ri <- ri[keep_r]; ci <- ci[keep_c]
    rg <- as.character(par$row_lookup[[par$guild_col]][ri])
    cg <- as.character(par$column_lookup[[par$guild_col]][ci])
    if (!length(rows) || !length(columns)) {
      mat <- matrix(numeric(), length(rows), length(columns), dimnames = list(rows, columns))
    } else {
      if (is.null(rownames(par$theta)) || is.null(colnames(par$theta))) stop("`theta` requires guild row and column names.")
      missing_rows <- setdiff(unique(rg), rownames(par$theta))
      missing_columns <- setdiff(unique(cg), colnames(par$theta))
      if (length(missing_rows) || length(missing_columns)) stop("Lookup guilds are missing from `theta` dimnames.")
      mat <- par$theta[rg, cg, drop = FALSE]
      dimnames(mat) <- list(rows, columns)
    }
    new_ecological_network(mat, "probability",
      list(cell_id = assemblage$cell_id, model = "SBM"), par$row_group, par$column_group)
  }
  new_inference_model(pred, "stochastic block model",
    list(row_lookup = row_lookup, column_lookup = column_lookup, theta = theta,
         row_group = row_group, column_group = column_group,
         species_col = species_col, guild_col = guild_col), "block_model")
}

#' Construct the original palm-mammal SBM backend
#' @description Compatibility wrapper around `block_model()`.
#' @param palm_lookup,mammal_lookup Species-to-guild lookup tables.
#' @inheritParams block_model
#' @export
sbm_model <- function(palm_lookup, mammal_lookup, theta, species_col = "species", guild_col = "guild") {
  block_model(palm_lookup, mammal_lookup, theta, "palms", "mammals", species_col, guild_col)
}

#' Build one local generic block-model network directly
#' @param rows,columns Species names in each network group.
#' @param cell_id Optional raster cell identifier.
#' @inheritParams block_model
#' @export
local_block_network <- function(rows, columns, row_lookup, column_lookup, theta,
                                row_group = "rows", column_group = "columns",
                                species_col = "species", guild_col = "guild", cell_id = NA_integer_) {
  species <- setNames(list(rows, columns), c(row_group, column_group))
  a <- new_assemblage(species, cell_id)
  infer_network(a, block_model(row_lookup, column_lookup, theta, row_group, column_group,
                               species_col, guild_col))$matrix
}

#' Build one local palm-mammal probability network directly
#' @description Compatibility wrapper around `local_block_network()`.
#' @param palms,mammals Species names in each network group.
#' @param cell_id Optional raster cell identifier.
#' @inheritParams sbm_model
#' @export
local_probability_network <- function(palms, mammals, palm_lookup, mammal_lookup, theta,
                                      species_col = "species", guild_col = "guild", cell_id = NA_integer_) {
  local_block_network(palms, mammals, palm_lookup, mammal_lookup, theta,
                      "palms", "mammals", species_col, guild_col, cell_id)
}

#' Compatibility wrapper for spatial network downscaling
#' @param palm_stack,mammal_stack Aligned species distribution rasters.
#' @param min_palms,min_mammals Minimum local richness in each group.
#' @inheritParams sbm_model
#' @export
downscale_networks <- function(palm_stack, mammal_stack, palm_lookup, mammal_lookup, theta,
                               species_col = "species", guild_col = "guild", min_palms = 1L, min_mammals = 1L) {
  assemblages <- assemble_communities(list(palms = palm_stack, mammals = mammal_stack), min_species = 1L)
  assemblages <- Filter(function(a) length(a$species$palms) >= min_palms && length(a$species$mammals) >= min_mammals, assemblages)
  infer_networks(assemblages, sbm_model(palm_lookup, mammal_lookup, theta, species_col, guild_col), palm_stack[[1]])
}

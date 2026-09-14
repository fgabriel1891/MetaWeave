#' Create a standard raster grid
#' @param extent Numeric vector `c(xmin, xmax, ymin, ymax)`.
#' @param resolution Cell resolution.
#' @param crs Coordinate reference system.
#' @export
create_standard_grid <- function(extent = c(-180, 180, -90, 90), resolution = 1, crs = "EPSG:4326") {
  terra::rast(xmin = extent[1], xmax = extent[2], ymin = extent[3], ymax = extent[4], resolution = resolution, crs = crs)
}

#' Create a bounding-box polygon
#' @param extent Numeric vector `c(xmin, xmax, ymin, ymax)`.
#' @param crs Coordinate reference system.
#' @export
make_bbox <- function(extent, crs = "EPSG:4326") {
  terra::as.polygons(terra::ext(extent), crs = crs)
}

#' List ESRI shapefiles
#' @param path Directory.
#' @param recursive Search recursively.
#' @export
list_shapefiles <- function(path, recursive = FALSE) {
  list.files(path, pattern = "\\.shp$", full.names = TRUE, recursive = recursive, ignore.case = TRUE)
}

#' Rasterize one range file
#' @param shp_file Shapefile path.
#' @param template Raster template.
#' @param touches Count all touched cells.
#' @export
rasterize_range_file <- function(shp_file, template, touches = TRUE) {
  v <- terra::project(terra::vect(shp_file), terra::crs(template))
  r <- terra::rasterize(v, template, field = 1, background = NA, touches = touches)
  names(r) <- tools::file_path_sans_ext(basename(shp_file))
  r
}

#' Rasterize rows grouped by species
#' @param x An sf or SpatVector object.
#' @param template Raster template.
#' @param species_col Species-name column.
#' @param species_keep Optional species subset.
#' @param touches Count all touched cells.
#' @export
rasterize_range_rows <- function(x, template, species_col, species_keep = NULL, touches = TRUE) {
  v <- terra::vect(x)
  vals <- terra::values(v)
  if (!species_col %in% names(vals)) stop("Species column not found.")
  spp <- as.character(vals[[species_col]])
  if (!is.null(species_keep)) spp <- intersect(unique(spp), species_keep) else spp <- unique(spp)
  out <- lapply(spp, function(sp) {
    one <- v[as.character(terra::values(v)[[species_col]]) == sp, ]
    one <- terra::project(one, terra::crs(template))
    r <- terra::rasterize(one, template, field = 1, background = NA, touches = touches)
    names(r) <- sp
    r
  })
  if (!length(out)) stop("No species remained after filtering.")
  terra::rast(out)
}

#' Build a species range stack
#' @param source Directory, vector file, sf, or SpatVector.
#' @param template Raster template.
#' @param species_col Species column for multi-species data.
#' @param species_keep Optional species subset.
#' @param region Optional crop region.
#' @param touches Count all touched cells.
#' @export
build_range_stack <- function(source, template, species_col = NULL, species_keep = NULL, region = NULL, touches = TRUE) {
  if (is.character(source) && length(source) == 1L && dir.exists(source)) {
    files <- list_shapefiles(source)
    if (!length(files)) stop("No shapefiles found.")
    out <- terra::rast(lapply(files, rasterize_range_file, template = template, touches = touches))
  } else {
    if (is.character(source)) source <- terra::vect(source)
    if (is.null(species_col)) stop("`species_col` is required for multi-species data.")
    out <- rasterize_range_rows(source, template, species_col, species_keep, touches)
  }
  if (!is.null(region)) out <- crop_filter_stack(out, region)
  out
}

#' Crop a stack and remove species absent from the region
#' @param stack Species raster stack.
#' @param region Crop region.
#' @param return_cropped Return cropped rather than original layers.
#' @export
crop_filter_stack <- function(stack, region, return_cropped = TRUE) {
  cropped <- terra::crop(stack, region)
  present <- as.numeric(terra::global(!is.na(cropped), "sum", na.rm = TRUE)[, 1]) > 0
  if (return_cropped) cropped[[present]] else stack[[present]]
}

#' Prepare palm and mammal range stacks
#' @param palm_source,mammal_source Range sources accepted by `build_range_stack()`.
#' @param mammal_species_col Species-name column in the mammal source.
#' @param mammal_species_keep Optional mammal species subset.
#' @param output_dir Optional directory in which to save the range stacks.
#' @inheritParams build_range_stack
#' @export
prepare_range_stacks <- function(palm_source, mammal_source, mammal_species_col, mammal_species_keep = NULL,
                                 template = create_standard_grid(), region = NULL, output_dir = NULL, touches = TRUE) {
  palms <- build_range_stack(palm_source, template, region = region, touches = touches)
  mammals <- build_range_stack(mammal_source, template, mammal_species_col, mammal_species_keep, region, touches)
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(palms, file.path(output_dir, "palm_ranges_grid_stack.RDS"))
    saveRDS(mammals, file.path(output_dir, "mammal_ranges_grid_stack.RDS"))
  }
  list(palms = palms, mammals = mammals)
}

#' Assemble local communities from aligned distribution rasters
#' @param distributions Named list of SpatRaster or distribution_collection objects.
#' @param min_species Minimum species per group.
#' @export
assemble_communities <- function(distributions, min_species = 1L) {
  if (!is.list(distributions) || is.null(names(distributions))) stop("`distributions` must be a named list.")
  rasters <- lapply(distributions, function(x) if (inherits(x, "distribution_collection")) x$data else x)
  ref <- rasters[[1]]
  ok <- vapply(rasters[-1], function(x) terra::compareGeom(ref, x, stopOnError = FALSE), logical(1))
  if (length(ok) && !all(ok)) stop("All distribution rasters must share geometry.")
  values <- lapply(rasters, terra::values, mat = TRUE)
  valid <- Reduce(`&`, lapply(values, function(v) rowSums(!is.na(v) & v != 0) >= min_species))
  cells <- which(valid)
  xy <- terra::xyFromCell(ref, cells)
  lapply(seq_along(cells), function(k) {
    i <- cells[k]
    spp <- Map(function(v, r) names(r)[which(!is.na(v[i, ]) & v[i, ] != 0)], values, rasters)
    new_assemblage(spp, cells[k], c(x = xy[k, 1], y = xy[k, 2]))
  })
}

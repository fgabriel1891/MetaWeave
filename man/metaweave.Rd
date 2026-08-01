\name{metaweave}
\alias{metaweave}
\alias{new_distribution_collection}
\alias{new_assemblage}
\alias{new_inference_model}
\alias{new_ecological_network}
\alias{new_network_collection}
\alias{new_spatial_result}
\alias{create_standard_grid}
\alias{make_bbox}
\alias{list_shapefiles}
\alias{rasterize_range_file}
\alias{rasterize_range_rows}
\alias{build_range_stack}
\alias{crop_filter_stack}
\alias{prepare_range_stacks}
\alias{probability_matrix_model}
\alias{reconstruct_networks}
\alias{assemble_communities}
\alias{block_model}
\alias{infer_network}
\alias{infer_networks}
\alias{sbm_model}
\alias{local_probability_network}
\alias{local_block_network}
\alias{downscale_networks}
\alias{simulate_interaction_counts}
\alias{simulate_network}
\alias{summarize_network}
\alias{summarize_networks}
\alias{rarefied_network_metric}
\alias{add_network_metrics}
\alias{map_metrics}
\alias{network_metrics_raster}
\alias{run_network_downscaling}
\alias{run_spatial_inference}
\alias{interaction_table_to_matrix}
\title{Spatial ecological network inference}
\description{
MetaWeave reconstructs local assemblages from distributions, applies pluggable
ecological inference models, simulates observations, summarizes networks, and
maps cell-level results. See the package vignette and README for the public API.
}
\usage{
new_distribution_collection(data, guild = NULL, metadata = NULL)
new_assemblage(species, cell_id = NA_integer_, coordinates = NULL)
new_inference_model(predict, name = "custom", parameters = list(), subclass = NULL)
new_ecological_network(matrix, type = c("probability", "binary", "weighted", "observed"), metadata = list(), row_group = "rows", column_group = "columns")
create_standard_grid(extent = c(-180, 180, -90, 90), resolution = 1, crs = "EPSG:4326")
assemble_communities(distributions, min_species = 1L)
infer_network(assemblage, model)
infer_networks(assemblages, model, template = NULL)
block_model(row_lookup, column_lookup, theta, row_group, column_group, species_col = "species", guild_col = "guild")
summarize_network(network, threshold = NULL)
summarize_networks(x, threshold = NULL)
map_metrics(x, metrics = NULL)
}
\details{
The stochastic block model is the first included inference backend. Custom
backends are created with \code{new_inference_model()} and must return an
\code{ecological_network}.
}
\keyword{package}

# metaweave 0.3.0

* Renames the package and framework from rangerast to MetaWeave.
* Positions the package around weaving regional ecological information into local spatial networks.
* Adds `probability_matrix_model()` for direct use of named species-pair probability matrices, including ensemble predictions from external models.

# rangerast 0.2.0

* Generalizes bipartite inference from hard-coded palm and mammal groups to arbitrary named ecological groups.
* Adds `block_model()` and `local_block_network()` as the generic block-model API.
* Adds `reconstruct_networks()` and `run_spatial_inference()` as generic end-to-end workflows.
* Ecological networks now record their row and column group identities.
* Network summaries now expose generic `row_richness` and `column_richness` metrics.
* Keeps `sbm_model()`, `local_probability_network()`, and the original downscaling wrappers for backward compatibility.
* Fixes network-index construction when coordinate values carry names.

# rangerast 0.1.0

* Introduces model-agnostic objects for distributions, assemblages, inference models, ecological networks, network collections, and spatial results.
* Adds spatial assembly and range-rasterization utilities.
* Adds the stochastic block model inference backend.
* Adds simulation, network summaries, spatial mapping, and compatibility wrappers for the original downscaling workflow.

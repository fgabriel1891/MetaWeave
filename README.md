# MetaWeave

<!-- badges: start -->
[![R-CMD-check](https://github.com/USERNAME/metaweave/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/USERNAME/metaweave/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`metaweave` is a model-agnostic framework for weaving regional distributions and ecological models into spatially explicit local networks.

The core workflow is:

```text
DistributionCollection -> Assemblage -> InferenceModel -> EcologicalNetwork
                       -> Simulation -> Metrics -> SpatialResult
```

## Installation

```r
# install.packages("remotes")
remotes::install_github("fgabriel1891/metaweave")
```

## Generic producer-consumer example

```r
library(metaweave)

producer_lookup <- data.frame(species = c("Producer_a", "Producer_b"), guild = c("P1", "P2"))
consumer_lookup <- data.frame(species = c("Consumer_x", "Consumer_y"), guild = c("C1", "C2"))
theta <- matrix(c(.8, .2, .3, .7), 2, dimnames = list(c("P1", "P2"), c("C1", "C2")))

assemblage <- new_assemblage(list(
  producers = c("Producer_a", "Producer_b"),
  consumers = c("Consumer_x", "Consumer_y")
))

model <- block_model(
  producer_lookup,
  consumer_lookup,
  theta,
  row_group = "producers",
  column_group = "consumers"
)
network <- infer_network(assemblage, model)
summarize_network(network)
```

For aligned distribution rasters, the complete generic workflow is:

```r
out <- run_spatial_inference(
  distributions = list(producers = producer_rasters, consumers = consumer_rasters),
  model = model,
  min_species = 2
)

out$result$index
plot(out$spatial$data)
```

Group names are user-defined; producers and consumers are only examples. Any
backend made with `new_inference_model()` can replace the block model, provided
it accepts an assemblage and returns an `ecological_network`.

External models that already return named species-pair probabilities, such as
ensemble-averaged SBM predictions, can be used without inventing a single
block-level `theta`:

```r
model <- probability_matrix_model(
  probability_matrix = ensemble_probabilities,
  row_group = "producers",
  column_group = "consumers",
  name = "cassandRa SBM ensemble"
)
```

## Scope

`metaweave` integrates spatial distributions, ecological inference, simulation, network analysis, and mapping. It is not intended to replace GIS software, species distribution models, or specialist network-analysis packages.

## Development

```r
devtools::document()
devtools::test()
devtools::check()
```

Before publishing, replace `USERNAME` and the placeholder maintainer email in `DESCRIPTION`, this README, and `CITATION.cff`.

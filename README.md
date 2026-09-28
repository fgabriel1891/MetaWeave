# metaweave

Spatial ecological network inference from species distributions and supplied
interaction models. Supports species-pair probabilities, stochastic block models,
and directed constrained network ensembles.

## Installation

Install the release from CRAN when available:

```r
install.packages("metaweave")
```

## A minimal example

```r
library(metaweave)
p <- matrix(c(0.8, 0.2, 0.3, 0.7), 2, byrow = TRUE,
            dimnames = list(c("plant_a", "plant_b"), c("animal_x", "animal_y")))
model <- probability_matrix_model(p, "plants", "animals")
site <- new_assemblage(list(plants = c("plant_a", "plant_b"), animals = "animal_x"))
network <- infer_network(site, model)
summarize_network(network)
```

For a complete runnable spatial example, use
`vignette("getting-started", package = "metaweave")`.
The tutorial uses synthetic data; publication case studies are distributed
separately. Use `citation("metaweave")` for the software citation.

Author, maintainer and copyright holder: Gabriel Munoz, Concordia University.
Contact: gmunozacevedo@proton.me. License: MIT.

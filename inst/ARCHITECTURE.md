# ADR-001: Model-agnostic spatial ecological inference

Status: accepted for v0.1.0

## Vision

MetaWeave is a framework for reconstructing and analysing spatially explicit ecological structure from incomplete ecological information. Rasterization and downscaling are workflows within that framework, not its defining abstractions.

## Decision

The package uses a staged object model:

`DistributionCollection -> Assemblage -> InferenceModel -> EcologicalNetwork -> Simulation -> Metrics -> SpatialResult`

An inference backend receives an assemblage and returns an ecological network. Downstream functions depend on the network contract, not the backend. `block_model()` is the first generic backend; trait, phylogenetic, Bayesian, machine-learning, and null models can implement the same interface. Ecological group names are supplied by users and are never assumed by the core inference machinery.

Ecological networks are not restricted to bipartite probability matrices. The base object records a network type and supports probability, binary, weighted, and observed representations.

## Responsibilities

- `distribution_collection`: where species may occur.
- `assemblage`: which species co-occur in a local unit.
- `inference_model`: how ecological relationships are inferred.
- `ecological_network`: the inferred or observed relationship structure.
- Simulation functions: observation and uncertainty processes.
- Metric functions: network properties independent of inference origin.
- `spatial_result`: attachment of cell-level results to geography.

## Boundaries

MetaWeave does not replace GIS software, species distribution models, specialist network-analysis packages, or plotting systems. It integrates their outputs through stable ecological objects and transformations.

## Consequences

The architecture requires slightly more explicit objects than a single analysis script. In return, inference methods become replaceable, metrics become comparable across models, and spatial workflows remain independent of any one study system.

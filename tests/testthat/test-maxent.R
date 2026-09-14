test_that("fixed-link ensembles preserve constraints and produce marginals", {
  assemblage <- new_assemblage(
    list(web = c("a", "b", "c", "d")),
    cell_id = 7L
  )
  model <- maxent_model(
    row_group = "web",
    constraint = "links",
    links = 3,
    ensemble_size = 40,
    self_links = FALSE,
    seed = 10,
    keep_ensemble = TRUE
  )

  network <- infer_network(assemblage, model)
  ensemble <- network$metadata$ensemble

  expect_s3_class(network, "probability_network")
  expect_true(network$directed)
  expect_false(network$self_links)
  expect_equal(diag(network$matrix), setNames(rep(0, 4), letters[1:4]))
  expect_true(all(vapply(ensemble, sum, numeric(1)) == 3))
  expect_equal(network$matrix, Reduce(`+`, ensemble) / length(ensemble))
})

test_that("connectance is converted to the nearest feasible local link count", {
  assemblage <- new_assemblage(list(
    consumers = c("c1", "c2"),
    resources = c("r1", "r2", "r3")
  ))
  model <- maxent_model(
    "consumers", "resources",
    constraint = "connectance",
    connectance = 0.5,
    ensemble_size = 20,
    seed = 2,
    keep_ensemble = TRUE
  )

  network <- infer_network(assemblage, model)

  expect_true(all(vapply(network$metadata$ensemble, sum, numeric(1)) == 3))
  expect_equal(sum(network$matrix), 3)
})

test_that("degree sampler preserves directed margins and forbids self-links", {
  species <- c("a", "b", "c", "d")
  degrees <- setNames(rep(1, 4), species)
  model <- maxent_model(
    row_group = "web",
    constraint = "degree_sequence",
    out_degree = degrees,
    in_degree = degrees,
    ensemble_size = 25,
    self_links = FALSE,
    seed = 3,
    burn_in = 50,
    thin = 10,
    keep_ensemble = TRUE
  )

  network <- infer_network(new_assemblage(list(web = species)), model)

  for (adjacency in network$metadata$ensemble) {
    expect_equal(rowSums(adjacency), degrees)
    expect_equal(colSums(adjacency), degrees)
    expect_equal(diag(adjacency), setNames(rep(0, 4), letters[1:4]))
  }
})

test_that("seeded maximum-entropy inference preserves global RNG state", {
  set.seed(811)
  state_before <- .Random.seed
  model <- maxent_model("web", constraint = "links", links = 2,
                        ensemble_size = 10, seed = 12)
  assemblage <- new_assemblage(list(web = c("a", "b", "c")), cell_id = 2L)

  first <- infer_network(assemblage, model)
  state_after <- .Random.seed
  second <- infer_network(assemblage, model)

  expect_identical(state_after, state_before)
  expect_equal(first$matrix, second$matrix)
})

test_that("different cell identifiers receive reproducible distinct streams", {
  model <- maxent_model("web", constraint = "links", links = 4,
                        ensemble_size = 30, seed = 15)
  a1 <- new_assemblage(list(web = letters[1:5]), cell_id = 1L)
  a2 <- new_assemblage(list(web = letters[1:5]), cell_id = 2L)

  n1 <- infer_network(a1, model)
  n2 <- infer_network(a2, model)

  expect_false(isTRUE(all.equal(n1$matrix, n2$matrix)))
  expect_equal(n1$matrix, infer_network(a1, model)$matrix)
})

test_that("infeasible maximum-entropy constraints fail informatively", {
  assemblage <- new_assemblage(list(web = c("a", "b")))

  too_many <- maxent_model("web", constraint = "links", links = 3,
                           self_links = FALSE)
  expect_error(infer_network(assemblage, too_many), "only 2 dyads")

  inconsistent <- maxent_model(
    "web",
    constraint = "degree_sequence",
    out_degree = c(a = 1, b = 1),
    in_degree = c(a = 1, b = 0)
  )
  expect_error(infer_network(assemblage, inconsistent), "sums.*differ")

  impossible <- maxent_model(
    "web",
    constraint = "degree_sequence",
    out_degree = c(a = 1, b = 1),
    in_degree = c(a = 2, b = 0),
    self_links = FALSE
  )
  expect_error(infer_network(assemblage, impossible), "Infeasible `in_degree`")
})

test_that("maximum-entropy backend works with infer_networks", {
  model <- maxent_model("web", constraint = "connectance", connectance = 0.5,
                        ensemble_size = 10, seed = 5)
  assemblages <- list(
    new_assemblage(list(web = c("a", "b", "c")), 1L, c(x = 0, y = 0)),
    new_assemblage(list(web = c("b", "c", "d")), 2L, c(x = 1, y = 0))
  )

  result <- infer_networks(assemblages, model)

  expect_s3_class(result, "network_collection")
  expect_length(result$networks, 2)
  expect_true(all(vapply(result$networks, function(x) x$directed, logical(1))))
})

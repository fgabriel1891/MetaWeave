test_that("SBM model maps guilds to probabilities", {
  p <- data.frame(species = c("p1", "p2"), guild = c("a", "b"))
  m <- data.frame(species = c("m1", "m2"), guild = c("x", "y"))
  theta <- matrix(c(.9, .2, .3, .8), 2, dimnames = list(c("a", "b"), c("x", "y")))
  a <- new_assemblage(list(palms = c("p1", "p2"), mammals = c("m1", "m2")))
  n <- infer_network(a, sbm_model(p, m, theta))
  expect_s3_class(n, "probability_network")
  expect_equal(unname(n$matrix), unname(theta))
})

test_that("generic block model uses arbitrary assemblage group names", {
  producers <- data.frame(species = c("p1", "p2"), guild = c("a", "b"))
  consumers <- data.frame(species = c("c1", "c2"), guild = c("x", "y"))
  theta <- matrix(c(.9, .2, .3, .8), 2,
                  dimnames = list(c("a", "b"), c("x", "y")))
  assemblage <- new_assemblage(list(
    producers = c("p1", "p2"),
    consumers = c("c1", "c2")
  ))
  model <- block_model(producers, consumers, theta, "producers", "consumers")

  network <- infer_network(assemblage, model)

  expect_equal(network$row_group, "producers")
  expect_equal(network$column_group, "consumers")
  expect_equal(unname(network$matrix), unname(theta))
})

test_that("species probability matrices can be used directly", {
  probabilities <- matrix(
    c(.8, .2, .4, .9), 2, byrow = TRUE,
    dimnames = list(c("p1", "p2"), c("c1", "c2"))
  )
  model <- probability_matrix_model(
    probabilities,
    row_group = "producers",
    column_group = "consumers"
  )
  assemblage <- new_assemblage(list(
    producers = c("p2", "absent_producer"),
    consumers = c("c2", "c1")
  ))

  network <- infer_network(assemblage, model)

  expect_equal(rownames(network$matrix), "p2")
  expect_equal(colnames(network$matrix), c("c2", "c1"))
  expect_equal(unname(network$matrix), matrix(c(.9, .4), 1))
  expect_equal(network$row_group, "producers")
  expect_equal(network$column_group, "consumers")
})

test_that("invalid species probability matrices are rejected", {
  invalid <- matrix(1.2, 1, 1, dimnames = list("p1", "c1"))
  expect_error(
    probability_matrix_model(invalid, "producers", "consumers"),
    "\\[0, 1\\]"
  )
})

test_that("network summary distinguishes possible and expected links", {
  n <- new_ecological_network(matrix(c(.25, .75, .5, .5), 2), "probability")
  x <- summarize_network(n)
  expect_equal(x$possible_links, 4)
  expect_equal(x$expected_links, 2)
  expect_equal(x$row_richness, 2)
  expect_equal(x$column_richness, 2)
})

test_that("custom inference backends obey the common contract", {
  pred <- function(model, assemblage) new_ecological_network(matrix(.5, 1, 1), "probability")
  model <- new_inference_model(pred, "test")
  out <- infer_network(new_assemblage(list(a = "x")), model)
  expect_s3_class(out, "ecological_network")
})

test_that("network collection index ignores coordinate value names", {
  pred <- function(model, assemblage) {
    new_ecological_network(matrix(.5, 1, 1), "probability")
  }
  model <- new_inference_model(pred, "test")
  assemblages <- list(
    new_assemblage(list(a = "x"), cell_id = 1L,
                    coordinates = c(x = 0.5, y = 1.5))
  )

  out <- infer_networks(assemblages, model)

  expect_equal(out$index$cell_id, 1L)
  expect_equal(out$index$x, 0.5)
  expect_equal(out$index$y, 1.5)
})

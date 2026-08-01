test_that("standard grid has requested geometry", {
  skip_if_not_installed("terra")
  r <- create_standard_grid(c(0, 2, 0, 2), 1)
  expect_equal(terra::ncell(r), 4)
})

test_that("assemblages are extracted from aligned rasters", {
  skip_if_not_installed("terra")
  r <- create_standard_grid(c(0, 2, 0, 1), 1)
  p <- c(r, r); names(p) <- c("p1", "p2"); terra::values(p) <- matrix(c(1, NA, 1, 1), 2, byrow = TRUE)
  m <- r; names(m) <- "m1"; terra::values(m) <- 1
  a <- assemble_communities(list(producers = p, consumers = m))
  expect_length(a, 2)
  expect_equal(a[[1]]$species$producers, "p1")
})

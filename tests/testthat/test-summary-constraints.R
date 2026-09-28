test_that("forbidden self-links use species identity, not diagonal position", {
  mat <- matrix(0.5, 2, 2, dimnames = list(c("a", "b"), c("b", "a")))
  net <- new_ecological_network(mat, row_group = "s", column_group = "s",
                                directed = TRUE, self_links = FALSE)
  s <- summarize_network(net)
  expect_equal(s$possible_links, 2)
  expect_equal(s$expected_links, 1)
  expect_equal(s$mean_probability, 0.5)
  expect_equal(summarize_network(net, threshold = 0)$expected_links, 2)
})
test_that("maxent connectance and summary use the same eligible pairs", {
  a <- new_assemblage(list(species = letters[1:4]))
  m <- maxent_model("species", connectance = 0.25, ensemble_size = 20,
                    self_links = FALSE, seed = 42)
  s <- summarize_network(infer_network(a, m))
  expect_equal(s$possible_links, 12)
  expect_equal(s$expected_links, 3)
  expect_equal(s$mean_probability, 0.25)
})
test_that("plain matrices and unrestricted networks retain their summaries", {
  mat <- matrix(c(0.2, 0.8, 0.4, 0.6), 2)
  expect_equal(summarize_network(mat)$possible_links, 4)
  expect_equal(summarize_network(mat)$expected_links, 2)
  expect_equal(summarize_network(new_ecological_network(mat)), summarize_network(mat))
})
test_that("unnamed same-group networks handle an empty eligible domain", {
  n <- new_ecological_network(matrix(0, 1, 1), row_group = "s", column_group = "s",
                              directed = TRUE, self_links = FALSE)
  s <- summarize_network(n)
  expect_equal(s$possible_links, 0)
  expect_equal(s$expected_links, 0)
  expect_true(is.na(s$mean_probability))
  n$row_group <- "other"
  expect_error(summarize_network(n), "Species names")
})
test_that("self-link matching supports rectangular overlapping groups", {
  mat <- matrix(0.5, 2, 3, dimnames = list(c("a", "b"), c("b", "c", "a")))
  n <- new_ecological_network(mat, self_links = FALSE)
  expect_equal(summarize_network(n)$possible_links, 4)
  expect_equal(summarize_network(n)$expected_links, 2)
})

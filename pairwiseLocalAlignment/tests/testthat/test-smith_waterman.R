test_that("returns error if not numeric", {
  input <- "string example"
  expect_error(smith_waterman(input, -5, -6, 'GAATC', 'CATACG'))
  })

test_that("returns error if not numeric", {
  input <- "string example"
  expect_error(smith_waterman(10, input, -6, 'GAATC', 'CATACG'))
})

test_that("returns error if not numeric", {
  input <- "string example"
  expect_error(smith_waterman(10, -5, input, 'GAATC', 'CATACG'))
})

test_that("returns error if not character", {
  input <- 12
  expect_error(smith_waterman(10, -5, -6, input, 'CATACG'))
})

test_that("returns error if not character", {
  input <- 12
  expect_error(smith_waterman(10, -5, -6, 'GAATC', input))
})

test_that("return warning if match < mismatch", {
  expect_warning(smith_waterman(-5, 5, -6, 'GAATC', 'CATACG'))
})

test_that("computing alignment of 'GAATC' and 'CATACG'", {
  expect_equal(smith_waterman(10, -5, -6, 'GAATC', 'CATACG'), list("AT-C", "ATAC"))
})

test_that("length of the resulting list = 2", {
  expect_length(smith_waterman(10, -5, -6, 'GAATC', 'CATACG'), 2)
})

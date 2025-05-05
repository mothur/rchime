# test rchime function

test_that("test uchime parameters ", {
  options <- rchime_options(
    minh = 0.6, abskew = 1.9, mindiv = 0.505,
    xn = 8.01, dn = 1.401, xa = 1.001,
    minchunk = 63, chunks = 5, idsmoothwindow = 33,
    maxp = 3,
    skipgaps = FALSE, skipgaps2 = FALSE
  )

  expect_equal(options$minh, 0.6)
  expect_equal(options$abskew, 1.9)
  expect_equal(options$mindiv, 0.505)
  expect_equal(options$xn, 8.01)
  expect_equal(options$dn, 1.401)
  expect_equal(options$xa, 1.001)
  expect_equal(options$minchunk, 63)
  expect_equal(options$chunks, 5)
  expect_equal(options$idsmoothwindow, 33)
  expect_equal(options$maxp, 3)
  expect_false(options$skipgaps)
  expect_false(options$skipgaps2)
})

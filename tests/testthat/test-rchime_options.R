# test rchime function

test_that("test uchime parameters ", {
  options <- rchime_options(
    minh = 0.6, abskew = 1.9, mindiv = 0.505,
    xn = 8.01, dn = 1.401, maxp = 4
  )

  expect_equal(options$minh, 0.6)
  expect_equal(options$abskew, 1.9)
  expect_equal(options$mindiv, 0.505)
  expect_equal(options$xn, 8.01)
  expect_equal(options$dn, 1.401)
  expect_equal(options$maxp, 4)
})

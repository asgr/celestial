test_that("cosmological distances have expected identities and ordering", {
  distances <- cosdist(z = c(0, 0.5, 1), H0 = 70, OmegaM = 0.3)

  expect_equal(distances$z, c(0, 0.5, 1))
  expect_equal(distances$CoDist[1], 0)
  expect_true(all(diff(distances$CoDist) > 0))
  expect_equal(distances$LumDist, (1 + distances$z) * distances$CoDist)
  expect_equal(distances$AngDist, distances$CoDist / (1 + distances$z))
})

test_that("cosmological distance functions validate redshifts", {
  expect_error(cosdist(z = -1), "z must be > -1")
  expect_error(cosdist(z = Inf), "z must be finite")
})

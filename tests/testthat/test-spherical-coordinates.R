test_that("spherical and Cartesian coordinates round-trip", {
  spherical <- cbind(long = c(0, 90, 180), lat = c(0, 30, -45),
                     radius = c(1, 2, 3))

  cartesian <- sph2car(spherical, deg = TRUE)
  expect_equal(car2sph(cartesian, deg = TRUE), spherical, tolerance = 1e-10)
})

test_that("equatorial and galactic coordinates round-trip", {
  equatorial <- cbind(RA = c(0, 123.4, 359), Dec = c(-30, 0, 70))

  expect_equal(gal2eq(eq2gal(equatorial)), equatorial, tolerance = 1e-10)
})

test_that("sky area handles full-sky and unit conversions", {
  full_sky <- skyarea(long = c(0, 360), lat = c(-90, 90), outunit = "sr")

  expect_equal(unname(full_sky["area"]), 4 * pi, tolerance = 1e-10)
  expect_equal(unname(full_sky["areafrac"]), 1)

  degrees <- skyarea(long = c(0, 10), lat = c(0, 10), outunit = "deg2")
  steradians <- skyarea(long = c(0, 10), lat = c(0, 10), outunit = "sr")
  expect_equal(unname(steradians["area"]),
               unname(degrees["area"]) * (pi / 180)^2,
               tolerance = 1e-10)
})

test_that("sky area rejects unsupported units and non-physical regions", {
  expect_error(skyarea(inunit = "hour"), "inunit must be")
  expect_error(skyarea(outunit = "degree"), "outunit must be")
  expect_error(skyarea(long = c(10, 0)), "non-physical")
})

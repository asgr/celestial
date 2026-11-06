test_that("decimal degrees convert to formatted DMS and HMS", {
  expect_equal(
    deg2dms(c(-27.345, 12.5), type = "cat", sep = ":"),
    c("-27:20:42.00", "+12:30:00.00")
  )
  expect_equal(
    deg2hms(c(182.5, 0), type = "cat", sep = ":"),
    c("12:10:00.00", "00:00:00.00")
  )
})

test_that("DMS values and strings convert back to decimal degrees", {
  expect_equal(
    dms2deg(d = c(27, 12), m = c(20, 30), s = c(42, 0), sign = c(-1, 1)),
    c(-27.345, 12.5)
  )
  expect_equal(
    dms2deg(c("-27:20:42", "+12:30:00")),
    c(-27.345, 12.5)
  )
  expect_equal(hms2deg(c("12:10:00", "00:00:00")), c(182.5, 0))
})

test_that("coordinate conversions reject values outside their domains", {
  expect_error(deg2dms(91), "All deg values")
  expect_error(dms2deg(d = 12, m = 60), "All m values")
  expect_error(hms2deg(h = 25), "All h values")
})

test_that("single-coordinate matching chooses the nearest in-radius source", {
  matches <- coordmatchsing(
    RAref = 0, Decref = 0,
    coordcompare = rbind(c(1, 0), c(0, 0)),
    rad = 2, radunit = "deg"
  )

  expect_equal(matches$ID, c(2, 1))
  expect_equal(unname(matches$bestmatch["compareID"]), 2)
  expect_equal(unname(matches$bestmatch["sep"]), 0, tolerance = 1e-10)
})

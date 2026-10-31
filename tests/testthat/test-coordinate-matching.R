test_that("coordinate matching returns matches and separations in requested units", {
  reference <- rbind(c(10, 0), c(30, 0))
  comparison <- rbind(c(10 + 1 / 3600, 0), c(80, 0))

  matches <- coordmatch(reference, comparison, rad = 2, radunit = "asec")

  expect_equal(matches$Nmatch, c(1, 0))
  expect_equal(matches$bestmatch$refID, 1)
  expect_equal(matches$bestmatch$compareID, 1)
  expect_equal(matches$bestmatch$sep, 1, tolerance = 1e-5)
})

test_that("single-coordinate matching chooses the nearest in-radius source", {
  matches <- coordmatchsing(
    RAref = 0, Decref = 0,
    coordcompare = rbind(c(1, 0), c(0, 0)),
    rad = 2, radunit = "deg"
  )

  expect_equal(matches$ID, c(2, 1))
  expect_equal(matches$bestmatch["compareID"], 2)
  expect_equal(matches$bestmatch["sep"], 0, tolerance = 1e-10)
})

test_that("linked coordinate groups combine connected objects", {
  groups <- group_links(rbind(c(10, 20), c(20, 30), c(40, 50)))

  expect_equal(groups, list(c(10, 20, 30), c(40, 50)))
})

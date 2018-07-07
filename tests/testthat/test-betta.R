test_that("betta isn't stupid", {
  n <- 25
  
  bb <- betta(c(rep(7, n), rep(8, n)),
              c(rep(2, n), rep(2, n)),
              cbind(1, c(rep(0, n), rep(1, n))))
  expect_is(bb, "list")
  
  expect_equal(bb$table[2,1] %>% unname, 1)
  
  expect_is(betta_pic(c(rep(7, n), rep(8, n)),
                      c(rep(2, n), rep(2, n)),
                      cbind(1, c(rep(0, n), rep(1, n)))), "ggplot")
  
})
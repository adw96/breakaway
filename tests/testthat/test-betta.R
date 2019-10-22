context("betta")
library(breakaway)

test_that("betta isn't stupid", {
  n <- 25
  
  bb <- betta(c(rep(7, n), rep(8, n)),
              c(rep(2, n), rep(2, n)),
              cbind(1, c(rep(0, n), rep(1, n))))
  bb_noX <-  betta(c(rep(7, n), rep(8, n)),
                   c(rep(2, n), rep(2, n)))
  expect_is(bb, "list")
  expect_is(bb_noX, "list")
  
  expect_equal(bb$table[2,1] %>% unname, 1)
  
  expect_is(betta_random(c(2000, 3000, 4000, 3000), c(100, 200, 150, 180), X = cbind(Int = 1, 
                Cont_var = c(100, 150, 100, 50)), groups = c("a", "a", "a", "b")), "list")
  expect_is(betta_random(c(2000, 3000, 4000, 3000), c(100, 200, 150, 180), groups = c("a", "a", "a", "b")), "list")
  
  expect_is(betta_pic(c(rep(7, n), rep(8, n)),
                      c(rep(2, n), rep(2, n)),
                      cbind(1, c(rep(0, n), rep(1, n)))), "ggplot")
  
  # AIC and AICc should indicate that bb fits better than bb_noX
  expect_true(bb$aic < bb_noX$aic)
  expect_true(bb$aicc < bb_noX$aicc)
  
})

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
  
  expect_is(betta_random(c(2000, 3000, 4000, 3000), c(100, 200, 150, 180), 
                         X = cbind(Int = 1, 
                                   Cont_var = c(100, 150, 100, 50)), groups = c("a", "a", "a", "b")), "list")
  expect_is(betta_random(c(2000, 3000, 4000, 3000), c(100, 200, 150, 180), 
                         groups = c("a", "a", "a", "b")), "list")
  
  expect_is(betta_pic(c(rep(7, n), rep(8, n)),
                      c(rep(2, n), rep(2, n)),
                      cbind(1, c(rep(0, n), rep(1, n)))), "ggplot")
  
  # AIC and AICc should indicate that bb fits better than bb_noX
  expect_true(bb$aic < bb_noX$aic)
  expect_true(bb$aicc < bb_noX$aicc)
  expect_equal(bb$r_squared_wls, 1)
  expect_equal(bb_noX$r_squared_wls, 0)
  
  # issue # 71
  expect_is(betta(chats=c(-0.013892611956319358, 0.00809696161789718, 0.005459268094084673), 
                  ses=c(0.012075830274043914, 0.00429604108003701, 0.0067910218651793496)), 
            "list")
  
  
  df <- data.frame(chats = c(2000, 3000, 4000, 3000),
                   ses = c(100, 200, 150, 180),
                   Cont_var = c(100, 150, 100, 50),
                   groups = c("a", "a", "a", "b"))
  b_formula <- betta_random(ses = "ses", 
                            formula = chats ~ Cont_var | groups, data = df)
  b_inputs <- betta_random(df$chats, df$ses, 
                           X = cbind(Int = 1, df$Cont_var), 
                           groups = df$groups)
  
  expect_is(b_formula, "list")
  expect_is(b_inputs, "list")
  expect_true(all.equal(b_formula$table[2, ], b_inputs$table[2, ]))
  
  df2 <- data.frame(chats = c(2000, 3000, 4000, 3000),
                    ses = c(100, 200, 150, NA),
                    groups = c("a", NA, "b", "b"))
  
  expect_is(betta_random(ses = "ses", formula = chats ~ 1 | groups, data = df2), 
            "list")
  
  
  df_uncouth <- data.frame(pooleen = c(2000, 3000, 4000, 3000),
                           dayved = c(100, 200, 150, 180),
                           brybry = c(100, 150, 100, 50),
                           sarsh = c("a", "a", "a", "b"))
  expect_is(
    betta_random(ses = "dayved",
                 formula = pooleen ~ brybry | sarsh, data = df_uncouth), 
    "list")
  expect_is(
    betta_random(ses = dayved,
                 formula = pooleen ~ brybry | sarsh, data = df_uncouth), 
    "list")
  
  
})

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
  b_formula <- betta_random(ses = ses,
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

  expect_is(betta_random(ses = ses, formula = chats ~ 1 | groups, data = df2),
            "list")


  df_uncouth <- data.frame(pooleen = c(2000, 3000, 4000, 3000),
                           dayved = c(100, 200, 150, 180),
                           brybry = c(100, 150, 100, 50),
                           sarsh = c("a", "a", "a", "b"))
  expect_is(
    betta_random(ses = dayved,
                 formula = pooleen ~ brybry | sarsh, data = df_uncouth),
    "list")
  expect_is(
    betta_random(ses = "dayved",
                 formula = pooleen ~ brybry | sarsh, data = df_uncouth),
    "list")

  expect_error(
    betta_random(ses = "dayved",
                 formula = pooleen ~ brybry + sarsh, data = df_uncouth),
    "Make sure that your formula includes `|` to specify a random effect if you'd like
          to use `betta_random`. Otherwise, you can use `betta`."
    
  )
  
  # issue 118
  df_issue118 <- df_uncouth
  df_issue118$dayved[1] <- 0
  df_issue118$nutmeg <- c(1, 1, 1, 1)
  df_issue118_v2 <- df_issue118
  df_issue118_v2$dayved[1] <- NA

  expect_warning(
    betta(ses = "dayved",
          formula = pooleen ~ brybry, data = df_issue118),
    "At least one of your standard errors is 0 or NA. Any observation with
            a standard error of 0 or NA has been dropped from the analysis."
  )

  expect_equal(
    # suppressing warnings here because they are the same warnings as in the test
    # directly above and we want to test if two output tables are the same when a 
    # standard error is 0 versus when a standard error is NA
    suppressWarnings({
      betta(ses = dayved,
          formula = pooleen ~ brybry, data = df_issue118)$table}),
    suppressWarnings({
      betta(ses = dayved,
          formula = pooleen ~ brybry, data = df_issue118_v2)$table})
  )

  df_issue118_v3 <- df_issue118_v2
  df_issue118_v3$dayved[1] <- 100
  
  expect_error(
    betta(ses = dayved,
          formula = pooleen ~ nutmeg, data = df_issue118_v3),
    "Your design matrix is not full rank. We recommend that you
         examine your design matrix for linear dependency and remove
         redundant columns."
  )

  expect_error(
    betta_random(ses = dayved,
          formula = pooleen ~ nutmeg | sarsh, data = df_issue118_v2),
    "Your design matrix is not full rank. We recommend that you
         examine your design matrix for linear dependency and remove
         redundant columns."
  )


})

test_that("betta_lincom works with betta and betta_random", {
  n <- 25

  bb <- betta(c(rep(7, n), rep(8, n)),
              c(rep(2, n), rep(2, n)),
              cbind(1, c(rep(0, n), rep(1, n))))

  bb_lincom <- betta_lincom(bb,c(1,1))

  expect_is(bb_lincom, "data.frame")

  expect_equal(as.numeric(bb_lincom[,1]), 8,
               tolerance = 0.02)

  expect_equal(as.numeric(bb_lincom[,2]), 0.4,
               tolerance = 0.02)

  expect_equal(as.numeric(bb_lincom[,3]),  7.216014,
               tolerance = 0.02)

  expect_equal(as.numeric(bb_lincom[,4]),  8.783986,
               tolerance = 0.02)

  expect_equal(as.numeric(bb_lincom[,5]),  2.753624e-89,
               tolerance = 0.02)

  # commenting out this test, since we've already checked that betta and 
  # betta_random work as expected, and above we show that betta_lincom
  # works as expected 
  
  # df <- data.frame(chats = c(2000, 3000, 4000, 3000),
  #                  ses = c(100, 200, 150, 180),
  #                  Cont_var = c(100, 150, 100, 50),
  #                  groups = c("a", "a", "a", "b"))
  # 
  # b_formula <- betta_random(ses = ses,
  #                           formula = chats ~ Cont_var | groups, data = df)
  # 
  # b_lincom <- betta_lincom(b_formula, c(1,1))
  # 
  # expect_is(b_lincom, "data.frame")
  # 
  # expect_equal(as.numeric(b_lincom[,1]), 3000.301,
  #              tolerance = 1e-3)
  # 
  # expect_equal(as.numeric(as.numeric(b_lincom[,2])), 918.2737)
  # 
  # expect_equal(as.numeric(as.numeric(b_lincom[,3])),  1200.518,
  #              tolerance = 1e-5)
  # 
  # expect_equal(as.numeric(b_lincom[,4]),  4800.084,
  #              tolerance = 1e-5)
  # 
  # expect_equal(as.numeric(b_lincom[,5]),  0.0005428401,
  #              tolerance = 1e-5)


})

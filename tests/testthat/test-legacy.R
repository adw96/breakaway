context("legacy")
library(breakaway)

test_that("Updated S3 outputs consistent", {
  skip_on_cran()

  wlrm_t <- wlrm_transformed(apples)
  wlrm_ut <- wlrm_untransformed(apples)
  cb <- chao_bunge(apples)
  chao <- chao1(apples, answers = T)
  chao_bc <- chao1_bc(apples, answers = T)
  ba <- breakaway(apples, answers = T, output = F, plot = F)
  kp <- kemp(apples, answers = T, output = F, plot = F)
  ba_nof1 <- breakaway_nof1(apples[-1, ], answers = T, output = F, plot = F)

  ### New format
  expect_equal(wlrm_t$estimate, 1178.603, tolerance = 1e-3)
  expect_equal(wlrm_t$error, 28.17513, tolerance = 1e-3)

  expect_equal(wlrm_ut$estimate, 1329.659, tolerance = 1e-3)
  expect_equal(wlrm_ut$error, 77.00265, tolerance = 1e-3)

  expect_equal(cb$estimate, 1376.942, tolerance = 1e-3)
  expect_equal(cb$error, 64.43054, tolerance = 1e-3)

  expect_equal(chao$estimate, 1241.286, tolerance = 1e-3)
  expect_equal(chao$error, 38.05458, tolerance = 1e-3)

  expect_equal(chao_bc$estimate, 1238.912, tolerance = 1e-3)
  expect_equal(chao_bc$error, 37.69171, tolerance = 1e-3)

  expect_equal(ba$estimate, 1552.416, tolerance = 1e-1)
  expect_equal(ba$error, 304.7071, tolerance = 1e-1)

  expect_equal(ba_nof1$estimate, 1500.401, tolerance = 1e-1)
  expect_equal(ba_nof1$error, 1341.351, tolerance = 1e-1)

  ### Kemp format
  expect_equal(kp$estimate, 1552.416, tolerance = 1e-1)
  expect_equal(kp$error, 304.7071, tolerance = 1e-1)

})

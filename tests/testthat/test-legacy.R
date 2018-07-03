test_that("Updated S3 outputs consistent", {
  
  wlrm_t <- wlrm_transformed(apples)
  wlrm_ut <- wlrm_untransformed(apples)
  cb <- chao_bunge(apples)
  chao <- chao1(apples, answers = T)
  chao_bc <- chao1_bc(apples, answers = T)
  ba <- breakaway(apples, answers = T, output = F, plot = F)
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
  
  ### Old format
  expect_equal(ba$est, 1552.416, tolerance = 1e-1)
  expect_equal(ba$seest, 304.7071, tolerance = 1e-1)
  
  expect_equal(ba_nof1$est %>% unname, 1500.401, tolerance = 1e-1)
  expect_equal(ba_nof1$seest, 1341.351, tolerance = 1e-1)
  
})

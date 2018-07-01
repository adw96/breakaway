test_that("Updated S3 outputs consistent", {
  wlrm_t <- wlrm_transformed(apples)
  wlrm_ut <- wlrm_untransformed(apples)
  
  expect_equal(wlrm_t$estimate, 1178.603, tolerance = 1e-3)
  expect_equal(wlrm_t$error, 28.17513, tolerance = 1e-3)
  
  # expect_equal(wlrm_ut$estimate, 1329.659)
  # expect_equal(wlrm_ut$error, 77.00265)
  # 
})

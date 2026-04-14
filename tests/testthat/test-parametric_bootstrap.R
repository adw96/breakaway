set.seed(1)
test_that("test_submodel returns a p-value", {
 #generate data
  predictor <- rnorm(20)
  ses <- rexp(20,.01)
  b0 <- 500
  b1 <- 100
  outcome = b0 + b1*predictor + rnorm(20,0,ses) + rnorm(20,0,10)

  betta_fit <- betta(chats = outcome,
                                   ses = ses,
                                   formula = chats ~ predictor,
                                   data = tibble::tibble(predictor = predictor,
                                                 chats = outcome,
                                                 ses = ses))

  parametric_bootstrap_test <- test_submodel(betta_fit,
                                             ~1,
                                             nboot = 100)

  expect_is(parametric_bootstrap_test$pval, "numeric")

})

test_that("p-value are at least approximately uniform under null (univariate fixed-effects model)", {
  #generate data
  pvals <-numeric(100)

  for(i in 1:100){
  predictor <- rnorm(20)
  ses <- rexp(20,.01)
  b0 <- 500
  b1 <- 0
  outcome = b0 + b1*predictor + rnorm(20,0,ses) + rnorm(20,0,10)

  betta_fit <- betta(chats = outcome,
                     ses = ses,
                     formula = chats ~ predictor,
                     data = tibble(predictor = predictor,
                                   chats = outcome,
                                   ses = ses))

  parametric_bootstrap_test <- suppressWarnings(test_submodel(betta_fit,
                                             ~1,
                                             nboot = 10))
  pvals[i] <- parametric_bootstrap_test$pval
}

expect_true(max(abs(sapply(pvals, function(x) mean(pvals <= x)) - pvals)) <=
             0.25)

})

# Omitting this test from automatic testing because it takes too long to run 
# enough iterations, but manually check when updating the betta_random or 
# test_submodel code

# test_that("p-value are at least approximately uniform under null (univariate mixed-effects model)", {
#   pvals <- as.numeric(rep(NA, 50))
# 
#   groups <- rep(1:5,each = 4)
# 
#   for(i in 1:50){
#     predictor <- rnorm(20)
#     ses <- rexp(20,.01)
#     group_effects <- rnorm(5,sd = 50)
#     group_effects <- rep(group_effects,each = 4)
#     b0 <- 500
#     b1 <- 0
#     outcome = b0 + b1*predictor + rnorm(20,0,ses) + group_effects +  rnorm(20,0,1)
# 
#     betta_fit <- betta_random(chats = outcome,
#                        ses = ses,
#                        groups = groups,
#                        formula = chats ~ predictor|group,
#                        data = tibble(predictor = predictor,
#                                      group = groups,
#                                      chats = outcome,
#                                      ses = ses))
# 
#     parametric_bootstrap_test <- suppressWarnings(test_submodel(betta_fit,
#                                                ~1,
#                                                nboot = 10))
#     pvals[i] <- parametric_bootstrap_test$pval
#   }
# 
# 
#   expect_equal(max(abs(sapply(pvals, function(x) mean(pvals <= x)) - pvals)),
#                0.3)
# 
# 
# })

test_that("F-statistic is numeric", {
  predictor <- rnorm(20)
  ses <- rep(.01,20)
  b0 <- 500
  b1 <- 0
  outcome = b0 + b1*predictor + rnorm(20,0,ses) + rnorm(20,0,100)

  betta_fit <- betta(chats = outcome,
                                   ses = ses,
                                   formula = chats ~ predictor,
                                   data = tibble::tibble(predictor = predictor,
                                                 chats = outcome,
                                                 ses = ses))

  L <- matrix(c(0,1),nrow = 1)

  expect_true( is.numeric(get_F_stat(fitted_betta = betta_fit,
             L)))

})


test_that("F-test returns a list", {
  #generate data
  predictor <- rnorm(20)
  ses <- rep(.01,20)
  b0 <- 500
  b1 <- 0
  outcome = b0 + b1*predictor + rnorm(20,0,ses) + rnorm(20,0,100)

  betta_fit <- betta(chats = outcome,
                     ses = ses,
                     formula = chats ~ predictor,
                     data = tibble::tibble(predictor = predictor,
                                   chats = outcome,
                                   ses = ses))

  L <- matrix(c(0,1),nrow = 1)

  expect_true( is.list(F_test(betta_fit,
                                     L,
                              method = "asymptotic")))
  expect_true(is.list(F_test(betta_fit,
                             L,
                             method = "bootstrap",
                             nboot = 10)))

})

test_that("F-test returns a list for betta random", {
  #generate data
  groups <- rep(1:5,each = 4)
  predictor <- rnorm(20)
  group_effects <- rnorm(5,sd = 5)
  group_effects <- rep(group_effects,each = 4)
  ses <- rexp(20,.01)
  b0 <- 500
  b1 <- 100
  outcome = b0 + b1*predictor + group_effects + rnorm(20,0,ses) + rnorm(20,0,10)
  
  betta_random_fit <- betta_random(chats = outcome,
                                   ses = ses,
                                   formula = chats ~ predictor|group,
                                   data = tibble(predictor = predictor,
                                                 chats = outcome,
                                                 group = groups,
                                                 ses = ses))
  
  
  L <- matrix(c(0,1),nrow = 1)
  
  expect_true( is.list(F_test(betta_random_fit,
                              L,
                              method = "asymptotic")))
  expect_true(is.list(F_test(betta_random_fit,
                             L,
                             method = "bootstrap",
                             nboot = 10)))
  
})

test_that("simulate_betta returns a reasonable thing", {
  predictor <- rnorm(20)
  ses <- rexp(20,.01)
  b0 <- 500
  b1 <- 100
  outcome = b0 + b1*predictor + rnorm(20,0,ses) + rnorm(20,0,10)
  
  betta_fit <- betta(chats = outcome,
                     ses = ses,
                     formula = chats ~ predictor,
                     data = tibble::tibble(predictor = predictor,
                                           chats = outcome,
                                           ses = ses))
  
  simulations <- simulate_betta(betta_fit,100)
  
  expect_equal(length(simulations), 100)
  expect_equal(length(simulations[[1]]), 20)
})

test_that("simulate_betta does a reasonable thing", {
  #generate data
  set.seed(345)
  predictor <- rnorm(20)
  ses <- rexp(20,.01)
  b0 <- 500
  b1 <- 100
  outcome = b0 + b1*predictor + rnorm(20,0,ses) + rnorm(20,0,10)

  betta_fit <- betta(chats = outcome,
                     ses = ses,
                     formula = chats ~ predictor,
                     data = tibble::tibble(predictor = predictor,
                                   chats = outcome,
                                   ses = ses))

  simulations <- simulate_betta(betta_fit,1e5)

  simulations <- do.call(rbind,simulations)
  
  expect_equal(betta_fit$table[1, 1] + betta_fit$table[2, 1]*predictor, 
               apply(simulations, 2, mean),
               tolerance = 0.05)

})

test_that("simulate_betta_random does a reasonable thing", {
  groups <- rep(1:5,each = 4)
  predictor <- rnorm(20)
  group_effects <- rnorm(5,sd = 5)
  group_effects <- rep(group_effects,each = 4)
  ses <- rexp(20,.01)
  b0 <- 500
  b1 <- 100
  outcome = b0 + b1*predictor + group_effects + rnorm(20,0,ses) + rnorm(20,0,10)

  betta_random_fit <- betta_random(chats = outcome,
                                   ses = ses,
                                   formula = chats ~ predictor|group,
                                   data = tibble(predictor = predictor,
                                                 chats = outcome,
                                                 group = groups,
                                                 ses = ses))


  simulations <- simulate_betta_random(betta_random_fit,1e5)

  simulations <- do.call(rbind,simulations)

  expect_equal(betta_random_fit$table[1,1] + betta_random_fit$table[2,1]*predictor,
               apply(simulations,2,mean),
               tolerance = 0.05)

})

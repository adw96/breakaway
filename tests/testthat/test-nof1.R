context("breakaway nof1")
library(breakaway)

# nof1 can't currently take phyloseq!

test_that("breakaway nof1 works", {
  data("apples")
  apples_nof1 <- apples[-1,]
  
  ### Sample 111 was a particular problem after change
  expect_is(apples_nof1 %>% breakaway_nof1, "alpha_estimate")
  
  issue55 <- tibble::tibble("freq" = c(2, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 21, 22, 23, 24, 25, 34, 38, 42, 43, 55, 57, 59, 70, 72, 73, 80, 89, 93, 94, 96, 100, 111, 139, 149, 188, 232, 254, 295, 366, 436, 509, 551, 1296, 1450, 1777, 1836, 2324, 2643, 6468),
                    "counts" = c(5, 7, 7, 7, 5, 6, 1, 2, 2, 4, 2, 4, 2, 1, 4, 1, 1, 3, 1, 1, 1, 2, 1, 2, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
  
  b55 <- issue55 %>% breakaway_nof1
  
  b55_add_3 <- issue55 %>% rbind(data.frame("freq" = 3, "counts" = 2)) %>% as.data.frame
  b55_add_3 <- b55_add_3[order(b55_add_3$freq), ]

  expect_false(b55$estimate %>% is.na)
  expect_false((b55_add_3 %>% breakaway_nof1)$estimate %>% is.na)
  
})

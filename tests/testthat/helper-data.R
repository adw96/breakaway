# data for testing
library(dplyr)
random_samples <- sample(x = 1:143, size = 5, replace = F)
tables <- apply(toy_otu_table[,random_samples], 2, make_frequency_count_table)
df1 <- cbind(c(1,2,3,6), c(8,4,2,2))
df2 <- cbind(c(3,6,9,12,15,16), c(3,2,1,15,1,5))
df3 <- cbind(c(1,2,5,6,7,9), c(7,8,2,4,1,2))

datasets <- c(list(hawaii,
                 apples,
                 df1,
                 df2,
                 df3), tables)
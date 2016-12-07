## A primitive development version of sample size calculations for
## betta (alpha diversity differences). 

## Is almost certainly possible to find a closed form expression for this. 
## This is a really, really primitive way of doing it. I will come back to
## the math approach ASAP. (Sorry, I have a thesis to write!)

## read in some toy data
test_otu_table <- read.table("/Users/adw96/Dropbox/ForAmy/otutable")
save_test <- test_otu_table
names(test_otu_table) <- paste("Seq", 1:length(names(save_test)), sep="")

test_otu_table <- t(test_otu_table)
head(test_otu_table)

## read in the covariate info
metadata <- read.csv("/Users/adw96/Dropbox/ForAmy/EHM_Metadata_161019.csv")
head(metadata)

table(metadata$Day)

## Suppose we just want to compare species richness on Day 2
## between babies exposed to Antibiotics and those that aren't
metadata[metadata$Day == 2, "Antibiotics_MomOrBaby_new"]
metadata[metadata$Day == 2, "Subject"]

## Subset the data
colnames(test_otu_table)
day2s <- which(unlist(lapply(strsplit(colnames(test_otu_table), "-d"), function(x) x[2])) == "002")
day_2_otus <- test_otu_table[, day2s]
head(day_2_otus)

## How diverse are these datasets? Let's look at frequency count tables
day_2_frequency_counts <- build_frequency_count_tables(day_2_otus)
day_2_richnesses <- unlist(lapply(day_2_frequency_counts, function(x) sum(x[,2])))
day_2_simpsons <- unlist(lapply(day_2_frequency_counts, breakaway::simpson))

day_2_df <- data.frame("Richness"=day_2_richnesses, "Simpson"=day_2_simpsons, "Antibiotic"=metadata[metadata$Day == 2, "Antibiotics_MomOrBaby_new"])

## Any visible difference 
plot(day_2_simpsons, pch=ifelse(day_2_df$Antibiotic, 16, 1))
plot(day_2_richnesses, pch=ifelse(day_2_df$Antibiotic, 16, 1))


## look at sample sizes
sample_size_figure <- function(...) {
  stop('Function not yet implemented!')
  return(2)
}

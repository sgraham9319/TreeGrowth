library(dplyr)

# Load neighborhood data
all_nbhds <- read.csv("Data/Output_data/neighborhoods.csv",
                      stringsAsFactors = F)

# Set seed to ensure consistent test set specification
set.seed(100)

# Assign each tree to a test set
test_assign <- all_nbhds %>%
  group_by(tree_id) %>%
  summarize(stand_id = stand_id[1],
            species = species[1]) %>%
  group_by(stand_id) %>%
  mutate(test_set = sample(rep(1:4, length.out = n()), n())) %>%
  ungroup()

# Confirm equal number of trees in each test set per stand
table(test_assign$stand_id, test_assign$test_set)

# Confirm relatively equal test set distribution per species
table(test_assign$species, test_assign$test_set)

# Join to neighborhood data
rand_nbhds <- all_nbhds %>%
  left_join(test_assign %>% select(tree_id, test_set), by = "tree_id")

# Write random training and test sets to .csv
for(set in 1:4){
  
  # Isolate this training and test set
  test <- rand_nbhds %>%
    filter(test_set == set)
  training <- rand_nbhds %>%
    filter(test_set != set)
  
  # Write to .csv
  write.csv(test, paste("Data/Output_data/rand_test", set, ".csv",
                        sep = ""), row.names = F)
  write.csv(training, paste("Data/Output_data/rand_training", set, ".csv",
                            sep = ""), row.names = F)
}

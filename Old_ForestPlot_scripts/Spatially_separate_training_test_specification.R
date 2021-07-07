library(dplyr)

# Load neighborhood data
all_nbhds <- read.csv("Data/Output_data/neighborhoods.csv",
                      stringsAsFactors = F)

#=============================================
# Specifying spatially separated training sets
#=============================================

# Assign quadrants of each stand to a test set
set.seed(500)
stand <- rep(unique(all_nbhds$stand_id), each = 4)
quadrant <- rep_len(1:4, length.out = length(stand))
test_id <- vector()
while(length(test_id) < length(stand)){
  test_id <- c(test_id, sample(1:4, 4))
}
test_sets <- data.frame(stand, quadrant, test_id)

# Create vector of unique forest stand ids to loop through
stands <- unique(all_nbhds$stand_id)

# Create each training/test split and write to .csv (this will take < 1 minute)
for(set in 1:4){
  
  # Create empty dataframes for new training and test set
  test <- all_nbhds[F,]
  training <- all_nbhds[F,]
  
  # Loop through stands
  for(i in 1:length(stands)){
    
    # Subset neighborhood data to one stand
    nbhds <- all_nbhds %>% filter(stand_id == stands[i])
    
    # Identify quadrant assigned to test
    test_quad <- test_sets %>%
      filter(stand == stands[i] & test_id == set) %>%
      pull(quadrant)
    
    # Extract training and test data
    if(test_quad == 1){
      new_test <- nbhds %>% filter(between(x_coord, 15, 42.5) &
                                     between(y_coord, 15, 42.5))
      new_train_a <- nbhds %>% filter(between(x_coord, 57.5, 85) &
                                        between(y_coord, 15, 85))
      new_train_b <- nbhds %>% filter(between(x_coord, 15, 57.49) &
                                        between(y_coord, 57.5, 85))
    } else if(test_quad == 2){
      new_test <- nbhds %>% filter(between(x_coord, 57.5, 85) &
                                     between(y_coord, 15, 42.5))
      new_train_a <- nbhds %>% filter(between(x_coord, 15, 85) &
                                        between(y_coord, 57.5, 85))
      new_train_b <- nbhds %>% filter(between(x_coord, 15, 42.5) &
                                        between(y_coord, 15, 57.49))
    } else if(test_quad == 3){
      new_test <- nbhds %>% filter(between(x_coord, 57.5, 85) &
                                     between(y_coord, 57.5, 85))
      new_train_a <- nbhds %>% filter(between(x_coord, 15, 42.5) &
                                        between(y_coord, 15, 85))
      new_train_b <- nbhds %>% filter(between(x_coord, 42.51, 85) &
                                        between(y_coord, 15, 42.5))
    } else if(test_quad == 4){
      new_test <- nbhds %>% filter(between(x_coord, 15, 42.5) &
                                     between(y_coord, 57.5, 85))
      new_train_a <- nbhds %>% filter(between(x_coord, 15, 85) &
                                        between(y_coord, 15, 42.5))
      new_train_b <- nbhds %>% filter(between(x_coord, 57.5, 85) &
                                        between(y_coord, 42.51, 85))
    }
    
    # Add to cumulative training and test
    test <- bind_rows(test, new_test)
    training <- bind_rows(training, new_train_a, new_train_b)
    
  }
  
  # Write to .csv
  write.csv(test, paste("Data/Output_data/test", set, ".csv", sep = ""),
            row.names = F)
  write.csv(training, paste("Data/Output_data/training", set, ".csv", sep = ""),
            row.names = F)
  
}
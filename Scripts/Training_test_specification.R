library(dplyr)


# Assign quadrants of each stand to a test set
set.seed(500)
stand <- rep(unique(complete_nbhds$stand_id), each = 4)
quadrant <- rep_len(1:4, length.out = length(stand))
test_id <- vector()
while(length(test_id) < length(stand)){
  test_id <- c(test_id, sample(1:4, 4))
}
test_sets <- data.frame(stand, quadrant, test_id)

# Create empty dataframes for training and test sets
train_id <- 1
training <- complete_nbhds[F,]
test <- complete_nbhds[F,]
stands <- unique(complete_nbhds$stand_id)

for(i in 1:length(stands)){
  
  # Subset neighborhood data to one stand
  nbhds <- complete_nbhds %>% filter(stand_id == stands[i])
  
  # Identify quadrant assigned to test
  test_quad <- test_sets %>%
    filter(stand == stands[i] & test_id == train_id) %>%
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




# Add quadrant variable
mapping <- mapping %>%
  mutate(quadrant = case_when(
    x_coord < 50 & y_coord >= 50 ~ 4,
    x_coord >= 50 & y_coord >= 50 ~ 3,
    x_coord >= 50 & y_coord < 50 ~ 2,
    TRUE ~ 1
    ))

# Create data frame with one row for each quadrant in each stand
test_sets <- mapping %>%
  group_by(stand_id, quadrant) %>%
  summarize() %>%
  ungroup()

# Randomly assign each quadrant in each stand to a test set
test_id <- vector()
while(length(test_id) < nrow(test_sets)){
  test_id <- c(test_id, sample(1:4, 4))
}
test_sets <- test_sets %>%
  mutate(test_id = test_id)

# Add test set designation to mapping
mapping <- mapping %>%
  left_join(test_sets, by = c("stand_id", "quadrant"))

# Isolate training set
train1 <- mapping %>%
  filter(test_id != 4)

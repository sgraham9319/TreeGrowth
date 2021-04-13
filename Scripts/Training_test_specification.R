library(dplyr)
library(ForestPlot)

#==================
# Part 1. Load data
#==================

# Load mapping data
mapping <- read.csv("Data/Raw_data/mapping_2017.csv", stringsAsFactors = F)

# Load tree measurement data
tree <- read.csv("Data/Raw_data/tree_growth_2017.csv", stringsAsFactors = F)

# Load plot abiotic data
abio <- read.csv("Data/Raw_data/stand_abiotic_data.csv", stringsAsFactors = F)

#=============================
# Part 2. Create neighborhoods
#=============================

# Define neighborhood radius
nb_rad <- 15

# Obtain all neighborhood data
neighbors <- neighborhoods(mapping = mapping, radius = nb_rad, densities = T)

# Remove focal trees whose neighborhood overlaps stand boundary
neighbors <- neighbors %>%
  filter(x_coord >= nb_rad & x_coord <= 100 - nb_rad &
           y_coord >= nb_rad & y_coord <= 100 - nb_rad)

# Remove small competitors
neighbors <- neighbors %>%
  filter(size_cat_comp == "regular")

#================================
# Part 3. Calculate annual growth
#================================

# Calculate annual growth for all trees
growth <- growth_summary(tree)

# Remove trees that were only measured once and/or had negative growth
growth <- growth %>%
  filter(growth$first_record != growth$last_record &
           annual_growth >= 0)

# Check for any remaining NA growth values
sum(is.na(growth$annual_growth))
sum(is.na(growth$size_corr_growth))

#========================================
# Part 4. Join datasets and write to .csv
#========================================

# Combine neighbors, growth, and plot abiotic data
complete_nbhds <- neighbors %>%
  # Inner join here because not all trees in neighbors have growth data
  inner_join(growth %>% select(tree_id, midpoint_size, annual_growth,
                               size_corr_growth),
             by = "tree_id") %>%
  left_join(abio, by = "stand_id")

# Write to .csv
write.csv(complete_nbhds, "Data/Output_data/neighborhoods.csv", row.names = F)



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

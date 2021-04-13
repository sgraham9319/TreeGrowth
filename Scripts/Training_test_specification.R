library(dplyr)
library(ForestPlot)

#=====================
# Part 1. Loading data
#=====================

# Load mapping data
mapping <- read.csv("Data/mapping_2017.csv", stringsAsFactors = F)

# Load tree measurement data
tree <- read.csv("Data/tree_growth_2017.csv", stringsAsFactors = F)

# Load plot abiotic data
abio <- read.csv("Data/stand_abiotic_data.csv", stringsAsFactors = F)

#===============================
# Part 2. Creating neighborhoods
#===============================

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

#==================================
# Part 3. Calculating annual growth
#==================================

# Calculate annual growth for all trees
growth <- growth_summary(tree)

# Remove trees that were only measured once and/or had negative growth
growth <- growth %>%
  filter(growth$first_record != growth$last_record &
           annual_growth >= 0)

# Check for any remaining NA growth values
sum(is.na(growth$annual_growth))
sum(is.na(growth$size_corr_growth))

#===========================================================
# Part 4. Combining neighbors, growth, and plot abiotic data
#===========================================================

complete_nbhds <- neighbors %>%
  inner_join(growth %>% select(tree_id, midpoint_size, annual_growth,
                               size_corr_growth),
             by = "tree_id") %>%
  left_join(abio, by = "stand_id")







# Add quadrant variable
mapping <- mapping %>%
  mutate(quadrant = case_when(
    x_coord >= 50 & y_coord >= 50 ~ 4,
    x_coord < 50 & y_coord >= 50 ~ 3,
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

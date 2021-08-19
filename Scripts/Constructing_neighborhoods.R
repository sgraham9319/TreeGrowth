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
neighbors <- neighborhoods(mapping = mapping, radius = nb_rad)

# Calculate densities
nbhd_summ <- neighborhood_summary(neighbors, id_column = "tree_id", radius = 15)

# Join neighborhoods and densities
neighbors <- neighbors %>%
  left_join(nbhd_summ, by = "tree_id")

# Remove focal trees whose neighborhood overlaps stand boundary
neighbors <- neighbors %>%
  filter(x_coord >= nb_rad & x_coord <= 100 - nb_rad &
           y_coord >= nb_rad & y_coord <= 100 - nb_rad)

# Add size category of competitors
neighbors <- neighbors %>%
  left_join(mapping %>% select(tree_id, size_cat),
            by = c("id_comp" = "tree_id"))

# Remove small competitors
neighbors <- neighbors %>%
  filter(size_cat == "regular")

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
  left_join(abio, by = "stand_id") %>%
  select(-size_cat)

# Write to .csv
write.csv(complete_nbhds, "Data/Output_data/neighborhoods.csv", row.names = F)

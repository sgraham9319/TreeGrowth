library(dplyr)

# Load mapping data
mapping <- read.csv("Data/mapping_2017.csv")

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

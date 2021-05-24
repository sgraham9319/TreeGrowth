library(dplyr)

# Load required functions
#source("/gscratch/stf/sgraham3/Rscripts/lkhd_fitting_functions.R")
#source("/gscratch/stf/sgraham3/Rscripts/nci.R")
source("Functions/lkhd_fitting_functions.R")
source("Functions/nci.R")

# Create mean square error calculation function
mse <- function(x){
  sum((x$observations - x$predictions) ^ 2) / nrow(x)
}

# Define focal species
focal_sps <- "TSME"

# Define number of folds for cross-validation
nfolds <- 10

# Load training data
training <- read.csv("Data/Output_data/training1.csv")

# Subset to focal species and remove unneeded columns
sing_sp <- training %>%
  arrange(tree_id) %>%
  filter(species == focal_sps) %>%
  select(tree_id, stand_id, species, dbh, prox, sps_comp, dbh_comp,
         annual_growth, pet_mm)

# Change units of variables to give them similar ranges
sing_sp <- sing_sp %>%
  mutate(
    dbh = dbh / 10,                     # cm to dm
    dbh_comp = dbh_comp / 10,           # cm to dm
    annual_growth = annual_growth * 10, # cm to mm
    pet_dm = pet_mm / 100               # mm to dm
  ) %>%
  select(-pet_mm)

# Extract annual growth of each focal individual
focals <- sing_sp %>%
  group_by(tree_id) %>%
  summarize(dbh = dbh[1], annual_growth = annual_growth[1])

# Assign each tree to a validation set
set.seed(10)
focals <- focals %>%
  mutate(val_set = sample(rep(1:nfolds, length.out = n()), n()))

# Add validation set ids to sing_sp
sing_sp <- sing_sp %>%
  left_join(focals %>% select(tree_id, val_set), by = "tree_id")

# Create set of starting values
X0 <- c(5, 15, 25)
Xb <- 2
gmax <- 1
pet_a <- 3
pet_b <- 2
sigma <- 5
starting_vals <- expand.grid(X0 = X0, Xb = Xb, gmax = gmax, pet_a = pet_a,
                             pet_b = pet_b, sigma = sigma)

# Define cross-validation set
cv <- 1

# Define training and validation sets
sing_sp_val <- sing_sp %>%
  filter(val_set == cv)
sing_sp_t <- setdiff(sing_sp, sing_sp_val)

# Define training and validation focals
focals_val <- focals %>%
  filter(val_set == cv)
focals_t <- setdiff(focals, focals_val)

# Try optimizing one time for TSME - takes about 7 minutes
#par <- as.vector(starting_vals[1,])
#fit <- optim(par, neg_log_lkhd, method = "SANN")

# Convert starting values data frame to list format
start_vals <- split(starting_vals, 1:nrow(starting_vals))

# Run optimization with mclapply - this will not work on a Windows machine
optim_output <- mclapply(start_vals, nll_opt)

# Format output as data frame
optim_vals <- data.frame(matrix(unlist(optim_output),
                                nrow = nrow(starting_vals),
                                byrow = T))
optim_vals <- optim_vals[, 1:(ncol(starting_vals) + 1)]
names(optim_vals) <- c(paste(names(starting_vals), "_opt", sep = ""), "NLL")

# Combine starting and optimized values
output <- cbind(starting_vals, optim_vals)

# Make growth predictions
growth_test <- growth_pred(sing_sp_val,
                           output[fit_comp$index[i], "X0_opt"],
                           output[fit_comp$index[i], "Xb_opt"],
                           output[fit_comp$index[i], "gmax_opt"],
                           output[fit_comp$index[i], "pet_a_opt"],
                           output[fit_comp$index[i], "pet_b_opt"])

# Combine predicted and observed growth
obs_pred <- focals_val %>%
  left_join(growth_test, by = c("tree_id" = "ids")) %>%
  rename(observations = annual_growth,
         predictions = pred_grow)

# Calculate and store mean square error
mse(obs_pred)




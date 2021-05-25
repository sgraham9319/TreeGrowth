library(dplyr)
library(parallel)

# Load required functions
source("/gscratch/stf/sgraham3/Rscripts/lkhd_fitting_functions.R")
source("/gscratch/stf/sgraham3/Rscripts/nci.R")

# Create mean square error calculation function
mse <- function(x){
  sum((x$observations - x$predictions) ^ 2) / nrow(x)
}

# Define focal species
focal_sps <- "TSME"

# Define number of folds for cross-validation
nfolds <- 2

# Load training data
training <- read.csv("/gscratch/stf/sgraham3/data/training1.csv",
                     stringsAsFactors = F)
#training <- read.csv("Data/Output_data/training1.csv",
#                     stringsAsFactors = F)

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

# Define cross-validation parameters
cv <- 1:nfolds
X0 <- 5
cv_par_df <- expand.grid(cv = cv, X0 = X0)

# Convert cross-validation parameters data frame to list
cv_par_list <- split(cv_par_df, 1:nrow(cv_par_df))

# Add required data frames to each list element
for(i in cv_par_df$cv){
  
  # Define training and validation sets
  sing_sp_val <- sing_sp %>%
    filter(val_set == i)
  sing_sp_t <- setdiff(sing_sp, sing_sp_val)
  
  # Define training and validation focal trees
  focals_val <- focals %>%
    filter(val_set == i)
  focals_t <- setdiff(focals, focals_val)
  
  # Add to cv_par_list
  cv_par_list[[i]] <- list(cv_par_list[[i]], sing_sp_val, sing_sp_t,
                           focals_val, focals_t)
  
  # Name list elements
  names(cv_par_list[[i]]) <- c("params", "sing_sp_val", "sing_sp_t",
                               "focals_val", "focals_t")
  
}

#Testing
#cv_par <- cv_par_list[[1]]

# Define function for running optimization
test_func <- function(cv_par){
  
  # Run optimization
  optim_output <- optim(par = c(cv_par[["params"]]["X0"], c(2, 1, 3, 2, 5)),
                        fn = neg_log_lkhd, method = "SANN")
  
  # Extract optimization output
  optim_vals <- c(cv_par[["params"]]["cv"], cv_par[["params"]]["X0"],
                  optim_output[["par"]], optim_output[["value"]], NA)
  
  # Convert list output to vector
  optim_vals <- unlist(optim_vals)
  
  # Add names to output
  names(optim_vals) <- c("cv_set", "X0_start", "X0", "Xb", "gmax", "pet_a",
                         "pet_b", "sigma", "NLL", "mse")
  
  # Make growth predictions
  growth_test <- growth_pred(sing_sp_val,
                             optim_vals["X0"],
                             optim_vals["Xb"],
                             optim_vals["gmax"],
                             optim_vals["pet_a"],
                             optim_vals["pet_b"])
  
  # Combine predicted and observed growth
  obs_pred <- focals_val %>%
    left_join(growth_test, by = c("tree_id" = "ids")) %>%
    rename(observations = annual_growth,
           predictions = pred_grow)
  
  # Calculate and store mean square error
  optim_vals["mse"] <- mse(obs_pred)
  
  # Return optimization results
  return(optim_vals)
  
}

# Run optimization with mclapply - this will not work on a Windows machine
optim_res_list <- mclapply(cv_par_list, test_func)

# Combine listed results into a data frame
results <- bind_rows(optim_res_list)

# Save results
write.csv(results, paste("/gscratch/stf/sgraham3/output/no_comp_cv_practice",
                           focal_sps, ".csv", sep = ""), row.names = F)
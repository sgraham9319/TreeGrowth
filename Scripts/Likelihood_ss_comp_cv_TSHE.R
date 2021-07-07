library(dplyr)
library(parallel)

# Load required functions
source("Functions/lkhd_fitting_functions.R")
source("Functions/nci.R")
#source("/gscratch/stf/sgraham3/Rscripts/lkhd_fitting_functions.R")
#source("/gscratch/stf/sgraham3/Rscripts/nci.R")

# Create mean square error calculation function
mse <- function(x){
  sum((x$observations - x$predictions) ^ 2) / nrow(x)
}

# Define focal species and training set
focal_sps <- "TSHE"
set <- 1

# Define number of folds for cross-validation
nfolds <- 10

# Load common competitors data and extract for focal species
comm_comp <- read.csv("Data/Output_data/common_comps.csv", stringsAsFactors = F)
#comm_comp <- read.csv("/gscratch/stf/sgraham3/data/common_comps.csv",
#                      stringsAsFactors = F)
comm_comp <- comm_comp[, focal_sps]

# Load training data
training <- read.csv(paste("Data/Output_data/rand_training", set, ".csv",
                           sep = ""), stringsAsFactors = F)
#training <- read.csv(paste("/gscratch/stf/sgraham3/data/rand_training", set,
#                           ".csv", sep = ""), stringsAsFactors = F)

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

# Change rare competitors to OTHR
sing_sp <- sing_sp %>%
  mutate(sps_comp = if_else(sps_comp %in% comm_comp, sps_comp, "OTHR"))

# Create vector of competitor species
comps <- sort(unique(sing_sp$sps_comp))

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
X0 <- c(25, 30)
cv_par_df <- expand.grid(cv = cv, X0 = X0)

# Convert cross-validation parameters data frame to list
cv_par_list <- split(cv_par_df, 1:nrow(cv_par_df))

# Add required data frames to each list element
for(i in 1:length(cv_par_list)){
  
  # Define training and validation sets
  sing_sp_val <- sing_sp %>%
    filter(val_set == cv_par_df$cv[i])
  sing_sp_t <- setdiff(sing_sp, sing_sp_val)
  
  # Define training and validation focal trees
  focals_val <- focals %>%
    filter(val_set == cv_par_df$cv[i])
  focals_t <- setdiff(focals, focals_val)
  
  # Add to cv_par_list
  cv_par_list[[i]] <- list(cv_par_list[[i]], sing_sp_val, sing_sp_t,
                           focals_val, focals_t)
  
  # Name list elements
  names(cv_par_list[[i]]) <- c("params", "sing_sp_val", "sing_sp_t",
                               "focals_val", "focals_t")
  
}

# Define function for running optimization
opt_func <- function(cv_par){
  
  # Run optimization
  optim_output <- optim(par = c(unlist(cv_par[["params"]]["X0"]), 
                                c(2, 1, 3, 2, 1, 1, 1, 0.5, 0.5, 0.5, 0.5,
                                  0.5, 0.5, 0.5, 0.5, 0.5, 5)),
                        fn = neg_log_lkhd, method = "SANN")
  
  # Extract optimization output
  optim_vals <- c(cv_par[["params"]]["cv"], cv_par[["params"]]["X0"],
                  optim_output[["par"]], optim_output[["value"]], NA)
  
  # Convert list output to vector
  optim_vals <- unlist(optim_vals)
  
  # Add names to output
  names(optim_vals) <- c("cv_set", "X0_start", "X0", "Xb", "gmax", "pet_a",
                         "pet_b", "C", "alpha", "beta", comps, 
                         "sigma", "NLL", "mse")
  
  # Make growth predictions
  growth_test <- growth_pred(nbhd_data = sing_sp_val,
                             X0 = optim_vals["X0"],
                             Xb = optim_vals["Xb"],
                             gmax = optim_vals["gmax"],
                             pet_a = optim_vals["pet_a"],
                             pet_b = optim_vals["pet_b"],
                             C = optim_vals["C"],
                             alpha = optim_vals["alpha"],
                             beta = optim_vals["beta"],
                             lmd1 = optim_vals[comps[1]],
                             lmd2 = optim_vals[comps[2]],
                             lmd3 = optim_vals[comps[3]],
                             lmd4 = optim_vals[comps[4]],
                             lmd5 = optim_vals[comps[5]],
                             lmd6 = optim_vals[comps[6]],
                             lmd7 = optim_vals[comps[7]],
                             lmd8 = optim_vals[comps[8]],
                             lmd9 = optim_vals[comps[9]])
  
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

# Try optimizing one time - will take > 30 minutes
#par_list_sub <- cv_par_list[[1]]
#fit <- opt_func(par_list_sub)

# Run optimization with mclapply - this will not work on a Windows machine
optim_res_list <- mclapply(cv_par_list, opt_func)

# Combine listed results into a data frame
results <- bind_rows(optim_res_list)

# Save results
write.csv(results, paste("Data/Output_data/ss_comp_cv", set, "_",
                         focal_sps, ".csv", sep = ""), row.names = F)
#write.csv(results, paste("/gscratch/stf/sgraham3/output/ss_comp_cv", set, "_",
#                         focal_sps, ".csv", sep = ""), row.names = F)
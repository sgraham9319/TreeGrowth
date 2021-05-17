library(dplyr)
library(parallel)

# Load required functions
source("/gscratch/stf/sgraham3/Rscripts/lkhd_fitting_functions.R")
source("/gscratch/stf/sgraham3/Rscripts/nci.R")

# Define AICc calculation function
AICc_calc <- function(k, NLL, n){
  (2 * (k + NLL)) + ((2 * ((k^2) + k)) / (n - k - 1))
}

# Define coefficient of determination function
coef_det <- function(x){
  1 - (sum((x$observations - x$predictions)^2) / 
         sum((x$observations - mean(x$observations))^2))
}

# Specify training set and focal species
set <- 1
focal_sps <- "THPL"

# Load common competitors data and extract for focal species
comm_comp <- read.csv("/gscratch/stf/sgraham3/data/common_comps.csv",
                      stringsAsFactors = F)
comm_comp <- comm_comp[, focal_sps]

# Load training data
training <- read.csv(paste("/gscratch/stf/sgraham3/data/rand_training", set,
                           ".csv", sep = ""), stringsAsFactors = F)

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

# Separate focal trees into training and cv sets
set.seed(50)
focals_cv <- focals[sample(1:nrow(focals), round(0.2 * nrow(focals))), ]
focals_t <- setdiff(focals, focals_cv)

# Separate sing_sp into training and cv sets
sing_sp_cv <- sing_sp %>%
  filter(tree_id %in% focals_cv$tree_id)
sing_sp_t <- sing_sp %>%
  filter(tree_id %in% focals_t$tree_id)

# Create empty data frame to store fit data
fit_data <- data.frame("X0" = numeric())

#----------------------------------
# Model structure 1: No competition
#----------------------------------

# Create set of starting values
X0 <- c(5, 15)
Xb <- 2
gmax <- 1
pet_a <- 3
pet_b <- 2
sigma <- 5
starting_vals <- expand.grid(X0 = X0, Xb = Xb, gmax = gmax, pet_a = pet_a,
                             pet_b = pet_b, sigma = sigma)
starting_vals <- bind_rows(starting_vals, starting_vals, starting_vals)

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

# Add AICc column
output$AICc <- AICc_calc(ncol(starting_vals), output$NLL, nrow(focals_t))

# Add empty train_r2 and cv_r2 columns
output$train_r2 <- rep(NA, times = nrow(output))
output$cv_r2 <- rep(NA, times = nrow(output))

# Make growth predictions and calculate coefficient of determination
for(i in 1:nrow(output)){
  
  # Predict training growth
  growth_t <- growth_pred(sing_sp_t, output[i, "X0_opt"],
                          output[i, "Xb_opt"], output[i, "gmax_opt"],
                          output[i, "pet_a_opt"], output[i, "pet_b_opt"])
  
  # Predict test growth
  growth_cv <- growth_pred(sing_sp_cv, output[i, "X0_opt"],
                           output[i, "Xb_opt"], output[i, "gmax_opt"],
                           output[i, "pet_a_opt"], output[i, "pet_b_opt"])
  
  # Combine predicted and observed growth
  obs_pred_t <- focals_t %>%
    left_join(growth_t, by = c("tree_id" = "ids")) %>%
    rename(observations = annual_growth,
           predictions = pred_grow)
  obs_pred_cv <- focals_cv %>%
    left_join(growth_cv, by = c("tree_id" = "ids")) %>%
    rename(observations = annual_growth,
           predictions = pred_grow)
  
  # Calculate coefficients of determination
  output$train_r2[i] <- coef_det(obs_pred_t)
  output$cv_r2[i] <- coef_det(obs_pred_cv)
  
}

# Add results to fit data
fit_data <- bind_rows(fit_data, output)

#------------------------------------------
# Model structure 2: Equivalent competition
#------------------------------------------

# Add additional parameters to starting_vals
starting_vals <- starting_vals %>%
  select(-sigma) %>%
  mutate(
    C = rep(1, times = nrow(starting_vals)),
    alpha = rep(1, times = nrow(starting_vals)),
    beta = rep(1, times = nrow(starting_vals)),
    sigma = rep(5, times = nrow(starting_vals)))

# Remake starting values list
start_vals <- split(starting_vals, 1:nrow(starting_vals))

# Run optimization with mclapply
optim_output <- mclapply(start_vals, nll_opt)

# Format output as data frame
optim_vals <- data.frame(matrix(unlist(optim_output),
                                nrow = nrow(starting_vals),
                                byrow = T))
optim_vals <- optim_vals[, 1:(ncol(starting_vals) + 1)]
names(optim_vals) <- c(paste(names(starting_vals), "_opt", sep = ""), "NLL")

# Combine starting and optimized values
output <- cbind(starting_vals, optim_vals)

# Add AICc column
output$AICc <- AICc_calc(ncol(starting_vals), output$NLL, nrow(focals_t))

# Add empty train_r2 and cv_r2 columns
output$train_r2 <- rep(NA, times = nrow(output))
output$cv_r2 <- rep(NA, times = nrow(output))

# Make growth predictions and calculate coefficient of determination
for(i in 1:nrow(output)){
  
  # Predict training growth
  growth_t <- growth_pred(sing_sp_t, output[i, "X0_opt"],
                          output[i, "Xb_opt"], output[i, "gmax_opt"],
                          output[i, "pet_a_opt"], output[i, "pet_b_opt"],
                          output[i, "C_opt"], output[i, "alpha_opt"],
                          output[i, "beta_opt"])
  
  # Predict test growth
  growth_cv <- growth_pred(sing_sp_cv, output[i, "X0_opt"],
                           output[i, "Xb_opt"], output[i, "gmax_opt"],
                           output[i, "pet_a_opt"], output[i, "pet_b_opt"],
                           output[i, "C_opt"], output[i, "alpha_opt"],
                           output[i, "beta_opt"])
  
  # Combine predicted and observed growth
  obs_pred_t <- focals_t %>%
    left_join(growth_t, by = c("tree_id" = "ids")) %>%
    rename(observations = annual_growth,
           predictions = pred_grow)
  obs_pred_cv <- focals_cv %>%
    left_join(growth_cv, by = c("tree_id" = "ids")) %>%
    rename(observations = annual_growth,
           predictions = pred_grow)
  
  # Calculate coefficients of determination
  output$train_r2[i] <- coef_det(obs_pred_t)
  output$cv_r2[i] <- coef_det(obs_pred_cv)
  
}

# Add results to fit data
fit_data <- bind_rows(fit_data, output)

#-----------------------------------------------
# Model structure 3: Intra vs. inter competition
#-----------------------------------------------

# Add additional parameters to starting_vals
starting_vals <- starting_vals %>%
  select(-sigma) %>%
  mutate(
    intra = rep(0.5, times = nrow(starting_vals)),
    inter = rep(0.5, times = nrow(starting_vals)),
    sigma = rep(5, times = nrow(starting_vals)))

# Remake starting values list
start_vals <- split(starting_vals, 1:nrow(starting_vals))

# Run optimization with mclapply
optim_output <- mclapply(start_vals, nll_opt)

# Format output as data frame
optim_vals <- data.frame(matrix(unlist(optim_output),
                                nrow = nrow(starting_vals),
                                byrow = T))
optim_vals <- optim_vals[, 1:(ncol(starting_vals) + 1)]
names(optim_vals) <- c(paste(names(starting_vals), "_opt", sep = ""), "NLL")

# Combine starting and optimized values
output <- cbind(starting_vals, optim_vals)

# Add AICc column
output$AICc <- AICc_calc(ncol(starting_vals), output$NLL, nrow(focals_t))

# Add empty train_r2 and cv_r2 columns
output$train_r2 <- rep(NA, times = nrow(output))
output$cv_r2 <- rep(NA, times = nrow(output))

# Make growth predictions and calculate coefficient of determination
for(i in 1:nrow(output)){
  
  # Predict training growth
  growth_t <- growth_pred(sing_sp_t, output[i, "X0_opt"],
                          output[i, "Xb_opt"], output[i, "gmax_opt"],
                          output[i, "pet_a_opt"], output[i, "pet_b_opt"],
                          output[i, "C_opt"], output[i, "alpha_opt"],
                          output[i, "beta_opt"], output[i, "intra_opt"],
                          output[i, "inter_opt"])
  
  # Predict test growth
  growth_cv <- growth_pred(sing_sp_cv, output[i, "X0_opt"],
                           output[i, "Xb_opt"], output[i, "gmax_opt"],
                           output[i, "pet_a_opt"], output[i, "pet_b_opt"],
                           output[i, "C_opt"], output[i, "alpha_opt"],
                           output[i, "beta_opt"], output[i, "intra_opt"],
                           output[i, "inter_opt"])
  
  # Combine predicted and observed growth
  obs_pred_t <- focals_t %>%
    left_join(growth_t, by = c("tree_id" = "ids")) %>%
    rename(observations = annual_growth,
           predictions = pred_grow)
  obs_pred_cv <- focals_cv %>%
    left_join(growth_cv, by = c("tree_id" = "ids")) %>%
    rename(observations = annual_growth,
           predictions = pred_grow)
  
  # Calculate coefficients of determination
  output$train_r2[i] <- coef_det(obs_pred_t)
  output$cv_r2[i] <- coef_det(obs_pred_cv)
  
}

# Add results to fit data
fit_data <- bind_rows(fit_data, output)

#------------------------------------------------
# Model structure 4: Species-specific competition
#------------------------------------------------

# Add additional parameters to starting_vals
new_params <- as.data.frame(matrix(0.5, nrow = nrow(starting_vals),
                                   ncol = length(comps)))
names(new_params) <- paste("lmd", 1:length(comps), sep = "")
starting_vals <- starting_vals %>%
  bind_cols(new_params) %>%
  select(-c(intra, inter, sigma)) %>%
  mutate(sigma = rep(5, times = nrow(starting_vals)))

# Remake starting values list
start_vals <- split(starting_vals, 1:nrow(starting_vals))

# Run optimization with mclapply
optim_output <- mclapply(start_vals, nll_opt)

# Format output as data frame
optim_vals <- data.frame(matrix(unlist(optim_output),
                                nrow = nrow(starting_vals),
                                byrow = T))
optim_vals <- optim_vals[, 1:(ncol(starting_vals) + 1)]
names(optim_vals) <- c(paste(names(starting_vals), "_opt", sep = ""), "NLL")

# Combine starting and optimized values
output <- cbind(starting_vals, optim_vals)

# Add AICc column
output$AICc <- AICc_calc(ncol(starting_vals), output$NLL, nrow(focals_t))

# Add empty train_r2 and cv_r2 columns
output$train_r2 <- rep(NA, times = nrow(output))
output$cv_r2 <- rep(NA, times = nrow(output))

# Make growth predictions and calculate coefficient of determination
for(i in 1:nrow(output)){
  
  # Predict training growth
  growth_t <- growth_pred(sing_sp_t, output[i, "X0_opt"],
                          output[i, "Xb_opt"], output[i, "gmax_opt"],
                          output[i, "pet_a_opt"], output[i, "pet_b_opt"],
                          output[i, "C_opt"], output[i, "alpha_opt"],
                          output[i, "beta_opt"], output[i, "lmd1_opt"],
                          output[i, "lmd2_opt"], output[i, "lmd3_opt"],
                          output[i, "lmd4_opt"], output[i, "lmd5_opt"])
  
  # Predict test growth
  growth_cv <- growth_pred(sing_sp_cv, output[i, "X0_opt"],
                           output[i, "Xb_opt"], output[i, "gmax_opt"],
                           output[i, "pet_a_opt"], output[i, "pet_b_opt"],
                           output[i, "C_opt"], output[i, "alpha_opt"],
                           output[i, "beta_opt"], output[i, "lmd1_opt"],
                           output[i, "lmd2_opt"], output[i, "lmd3_opt"],
                           output[i, "lmd4_opt"], output[i, "lmd5_opt"])
  
  # Combine predicted and observed growth
  obs_pred_t <- focals_t %>%
    left_join(growth_t, by = c("tree_id" = "ids")) %>%
    rename(observations = annual_growth,
           predictions = pred_grow)
  obs_pred_cv <- focals_cv %>%
    left_join(growth_cv, by = c("tree_id" = "ids")) %>%
    rename(observations = annual_growth,
           predictions = pred_grow)
  
  # Calculate coefficients of determination
  output$train_r2[i] <- coef_det(obs_pred_t)
  output$cv_r2[i] <- coef_det(obs_pred_cv)
  
}

# Add results to fit data
fit_data <- bind_rows(fit_data, output)

# Write results to csv
write.csv(fit_data, paste("/gscratch/stf/sgraham3/output/cv",
                        set, "_", focal_sps, ".csv", sep = ""), row.names = F)

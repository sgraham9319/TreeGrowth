library(dplyr)
library(parallel)

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
focal_sps <- "TSME"

# Load training data
#training <- read.csv(paste("Data/Output_data/rand_training", set,
#                           ".csv", sep = ""), stringsAsFactors = F)
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

# Extract annual growth of each focal individual
focals <- sing_sp %>%
  group_by(tree_id) %>%
  summarize(dbh = dbh[1], annual_growth = annual_growth[1])

# Separate focal trees into training and cv sets
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

# Create growth prediction function
growth_pred <- function(nbhd_data, X0, Xb, gmax, pet_a, pet_b){
  
  # Get list of focal tree ids
  ids <- unique(nbhd_data$tree_id)
  
  # Create vector to store growth predictions
  pred_grow <- rep(NA, times = length(ids))
  
  # Loop through focal trees
  for(id in 1:length(ids)){
    
    # Isolate neighborhood data for one focal tree
    neighbors <- nbhd_data %>%
      filter(tree_id == ids[id])
    
    # Predict growth
    pred_grow[id] <- gmax * 
      exp((-0.5) * (log(neighbors$dbh[1] / X0) / Xb) ^ 2) *
      exp((-0.5) * ((neighbors$pet_dm[1] - pet_a) / pet_b) ^ 2)
    
  }
  
  # Join growth predictions with tree ids and output result
  results <- data.frame(ids, pred_grow, stringsAsFactors = F)
  return(results)
  
}

# Define negative log likelihood function to be minimized
neg_log_lkhd <- function(par){
  
  # Determine model structure using length of par
  if(length(par) == 6){
    mod_str <- "no_comp"
  } else if(length(par) == 9){
    mod_str <- "eq_comp"
  } else if(length(par) == 11){
    mod_str <- "int_comp"
  } else{
    mod_str <- "ss_comp"
    num_comps <- length(par) - 9
  }
  
  # Define parameters
  X0 <- par[1]
  Xb <- par[2]
  gmax <- par[3]
  pet_a <- par[4]
  pet_b <- par[5]
  if(mod_str != "no_comp"){
    C <- par[6]
    alpha <- par[7]
    beta <- par[8]
  }
  if(mod_str == "int_comp"){
    intra <- par[9]
    inter <- par[10]
  }
  if(mod_str == "ss_comp"){
    lmd1 <- par[9]
    lmd2 <- par[10]
    lmd3 <- par[11]
    lmd4 <- par[12]
    if(num_comps > 4){
      lmd5 <- par[13]
    }
    if(num_comps > 5){
      lmd6 <- par[14]
      lmd7 <- par[15]
      lmd8 <- par[16] 
    }
    if(num_comps == 9){
      lmd9 <- par[17]
    }
  }
  sigma <- par[length(par)]
  
  # Prevent parameter values from becoming nonsensical
  if(sigma < 0) {return(Inf)}
  if(X0 < 0 | X0 > 30) {return(Inf)}
  if(Xb < 0 | Xb > 3) {return(Inf)}
  if(gmax < 0) {return(Inf)}
  if(pet_a < 2 | pet_a > 6) {return(Inf)}
  if(pet_b < 0 | pet_b > 3) {return(Inf)}
  if(mod_str != "no_comp"){
    if(C < 0 | C > 10) {return(Inf)}
    if(alpha < 0 | alpha > 4) {return(Inf)}
    if(beta < 0 | beta > 4) {return(Inf)}
  }
  if(mod_str == "int_comp"){
    if(intra < 0 | intra > 1) {return(Inf)}
    if(inter < 0 | inter > 1) {return(Inf)}
  }
  if(mod_str == "ss_comp"){
    if(lmd1 < 0 | lmd1 > 1) {return(Inf)}
    if(lmd2 < 0 | lmd2 > 1) {return(Inf)}
    if(lmd3 < 0 | lmd3 > 1) {return(Inf)}
    if(lmd4 < 0 | lmd4 > 1) {return(Inf)}
    if(num_comps > 4){
      if(lmd5 < 0 | lmd5 > 1) {return(Inf)}
    }
    if(num_comps > 5){
      if(lmd6 < 0 | lmd6 > 1) {return(Inf)}
      if(lmd7 < 0 | lmd7 > 1) {return(Inf)}
      if(lmd8 < 0 | lmd8 > 1) {return(Inf)} 
    }
    if(num_comps == 9){
      if(lmd9 < 0 | lmd9 > 1) {return(Inf)}
    }
  }
  
  # Make growth predictions
  if(mod_str == "no_comp"){
    pred <- growth_pred(sing_sp_t, X0, Xb, gmax, pet_a, pet_b)
  } else if(mod_str == "eq_comp"){
    pred <- growth_pred(sing_sp_t, X0, Xb, gmax, pet_a, pet_b, C, alpha, beta)
  } else if(mod_str == "int_comp"){
    pred <- growth_pred(sing_sp_t, X0, Xb, gmax, pet_a, pet_b, C, alpha,
                        beta, intra, inter)
  } else if(mod_str == "ss_comp"){
    if(num_comps == 4){
      pred <- growth_pred(sing_sp_t, X0, Xb, gmax, pet_a, pet_b, C, alpha, beta,
                          lmd1, lmd2, lmd3, lmd4)
    } else if(num_comps == 5){
      pred <- growth_pred(sing_sp_t, X0, Xb, gmax, pet_a, pet_b, C, alpha, beta,
                          lmd1, lmd2, lmd3, lmd4, lmd5)
    } else if(num_comps == 8){
      pred <- growth_pred(sing_sp_t, X0, Xb, gmax, pet_a, pet_b, C, alpha, beta,
                          lmd1, lmd2, lmd3, lmd4, lmd5, lmd6, lmd7, lmd8)
    } else if(num_comps == 9){
      pred <- growth_pred(sing_sp_t, X0, Xb, gmax, pet_a, pet_b, C, alpha, beta,
                          lmd1, lmd2, lmd3, lmd4, lmd5, lmd6, lmd7, lmd8, lmd9)
    }
  }
  
  # Join predictions to observations by tree_id
  combined <- left_join(focals_t, pred, by = c("tree_id" = "ids"))
  
  # Calculate negative log likelihood
  NLL <- -sum(dnorm(combined$annual_growth, mean = combined$pred_grow,
                    sd = sigma, log = T))
  
  # Return value
  return(NLL)
  
}

# Create set of starting values
X0 <- 15#c(5, 15, 25)
Xb <- 2
gmax <- 1
pet_a <- 3
pet_b <- 2
sigma <- 5
starting_vals <- expand.grid(X0 = X0, Xb = Xb, gmax = gmax, pet_a = pet_a,
                             pet_b = pet_b, sigma = sigma)
starting_vals <- bind_rows(starting_vals, starting_vals)

# Convert starting values data frame to list format
start_vals <- split(starting_vals, 1:nrow(starting_vals))

# Create function for running optimization
nll_opt <- function(par_list){
  optim(par = par_list[1,], fn = neg_log_lkhd, method = "SANN")
}

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

# Write results to csv
write.csv(fit_data, paste("/gscratch/stf/sgraham3/output/cv",
                        set, "_", focal_sps, ".csv", sep = ""), row.names = F)
library(tidyverse)

# Define AICc function
AICc_calc <- function(k, NLL, n){
  (2 * (k + NLL)) + ((2 * ((k^2) + k)) / (n - k - 1))
}

# Define coefficient of determination function
coef_det <- function(x){
  1 - (sum((x$observations - x$predictions)^2) / 
         sum((x$observations - mean(x$observations))^2))
}

# Define focal species and model structures
focal_sps <- c("ABAM", "CANO", "PSME", "THPL", "TSHE", "TSME")
model_str <- c("no_comp", "eq_comp", "int_comp", "ss_comp")

# Create data frame to store model selection results
mod_sel <- data.frame(
  species = rep(focal_sps, each = length(model_str)),
  mod_str = rep(model_str, times = length(focal_sps)),
  mse = NA
)

# Loop through focal species
for(sps in 1:length(focal_sps)){
  
  # Loop through model structures
  for(mod in 1:length(model_str)){
    
    # Load optimization results
    results <- read.csv(paste("Data/Output_data/", model_str[mod], "_cv_",
                              focal_sps[sps], ".csv", sep = ""))
    
    # Subset to best model for each cross-validation set
    results <- results %>%
      arrange(mse) %>%
      group_by(cv_set) %>%
      filter(row_number() == 1)
    
    # Calculate and store mean square error across cross-validation sets
    mod_sel[(sps - 1) * 4 + mod, "mse"] <- mean(results$mse)
    
  }
}

# Find model structure with lowest mse for each species
best_mod <- mod_sel %>%
  arrange(mse) %>%
  group_by(species) %>%
  filter(row_number() == 1) %>%
  arrange(species) %>%
  mutate(test_r2 = NA)


#=====================
# Calculating test fit
#=====================

# Load training data
training <- read.csv("Data/Output_data/rand_training1.csv")

# Load test data
test <- read.csv("Data/Output_data/rand_test1.csv")

# For testing, define species
sps <- 4
for(sps in 1:length(focal_sps)){
  
  # Calculate number of focal trees in training set
  sing_sp <- training %>%
    filter(species == focal_sps[sps])
  nfocals <- length(unique(sing_sp$tree_id))
  
  # Subset test data to focal species and remove unneeded columns
  ss_test <- test %>%
    arrange(tree_id) %>%
    filter(species == focal_sps[sps]) %>%
    select(tree_id, stand_id, species, dbh, prox, sps_comp, dbh_comp,
           annual_growth, pet_mm)
  
  # Change units of variables to give them similar ranges
  ss_test <- ss_test %>%
    mutate(
      dbh = dbh / 10,                     # cm to dm
      dbh_comp = dbh_comp / 10,           # cm to dm
      annual_growth = annual_growth * 10, # cm to mm
      pet_dm = pet_mm / 100               # mm to dm
    ) %>%
    select(-pet_mm)
  
  # Extract annual growth of each focal individual
  test_focals <- ss_test %>%
    group_by(tree_id) %>%
    summarize(dbh = dbh[1], annual_growth = annual_growth[1],
              pet_dm = pet_dm[1])
  
  # Load model output
  output <- read.csv(paste("Data/Output_data/", best_mod$mod_str[sps], "_rand1_",
                           focal_sps[sps], ".csv", sep = ""))
  
  # Calculate AICc for each model
  output$AICc <- AICc_calc(length(grep("_opt", names(output))),
                           output$NLL, nfocals)
  
  # Extract best model parameters by AICc
  params <- unlist(output %>%
                     arrange(AICc) %>%
                     filter(row_number() == 1) %>%
                     select(grep("_opt", names(output))) %>%
                     select(-"sigma_opt"))
  
  # Make growth predictions (required functions depends on model structure)
  if(best_mod$mod_str[sps] == "no_comp"){
    
    # Define growth prediction function
    growth_pred <- function(nbhd_data, X0, Xb, gmax, pet_a, pet_b){
      ids <- unique(nbhd_data$tree_id)
      pred_grow <- rep(NA, times = length(ids))
      for(id in 1:length(ids)){
        neighbors <- nbhd_data %>%
          filter(tree_id == ids[id])
        pred_grow[id] <- gmax * 
          exp((-0.5) * (log(neighbors$dbh[1] / X0) / Xb) ^ 2) *
          exp((-0.5) * ((neighbors$pet_dm[1] - pet_a) / pet_b) ^ 2)
      }
      results <- data.frame(ids, pred_grow, stringsAsFactors = F)
      return(results)
    }
    
    # Predict test growth
    growth_ts <- growth_pred(ss_test, params["X0_opt"],
                             params["Xb_opt"], params["gmax_opt"],
                             params["pet_a_opt"], params["pet_b_opt"])
    
  } else if(best_mod$mod_str[sps] == "int_comp"){
    
    # Define NCI function
    nci <- function(neighbors, alpha, beta, intra, inter){
      raw <- (neighbors$dbh_comp ^ alpha) / (neighbors$prox ^ beta)
      cons <- which(neighbors$sps_comp == focal_sps[sps])
      hets <- which(neighbors$sps_comp != focal_sps[sps])
      nci_con <- raw[cons] * intra
      nci_het <- raw[hets] * inter
      return(sum(c(nci_con, nci_het)))
    }
    
    # Define growth prediction function
    growth_pred <- function(nbhd_data, X0, Xb, gmax, pet_a, pet_b,
                            C, alpha, beta, intra, inter){
      ids <- unique(nbhd_data$tree_id)
      pred_grow <- rep(NA, times = length(ids))
      for(id in 1:length(ids)){
        neighbors <- nbhd_data %>%
          filter(tree_id == ids[id])
        pred_grow[id] <- gmax * 
          exp((-0.5) * (log(neighbors$dbh[1] / X0) / Xb) ^ 2) *
          exp((-0.5) * ((neighbors$pet_dm[1] - pet_a) / pet_b) ^ 2) *
          exp(-C * nci(neighbors, alpha, beta, intra, inter))
      }
      results <- data.frame(ids, pred_grow, stringsAsFactors = F)
      return(results)
    }
    
    # Predict test growth
    growth_ts <- growth_pred(ss_test, params["X0_opt"],
                             params["Xb_opt"], params["gmax_opt"],
                             params["pet_a_opt"], params["pet_b_opt"],
                             params["C_opt"], params["alpha_opt"],
                             params["beta_opt"], params["intra_opt"],
                             params["inter_opt"])
    
  } else if(best_mod$mod_str[sps] == "ss_comp"){
    
    # Load list of common competitors
    comm_comp <- read.csv("Data/Output_data/common_comps.csv",
                          stringsAsFactors = F)
    comm_comp <- comm_comp[, focal_sps[sps]]
    
    # Convert rare competitors to OTHR in test data
    ss_test <- ss_test %>%
      mutate(sps_comp = if_else(sps_comp %in% comm_comp, sps_comp, "OTHR"))
    
    # Create vector of competitor species
    comps <- sort(unique(ss_test$sps_comp))
    
    # Define functions based on number of competitors
    if(length(comps) == 4){
      
      # Define NCI function
      nci <- function(neighbors, alpha, beta, lmd1, lmd2, lmd3, lmd4){
        raw <- (neighbors$dbh_comp ^ alpha) / (neighbors$prox ^ beta)
        nci1 <- raw[which(neighbors$sps_comp == comps[1])] * lmd1
        nci2 <- raw[which(neighbors$sps_comp == comps[2])] * lmd2
        nci3 <- raw[which(neighbors$sps_comp == comps[3])] * lmd3
        nci4 <- raw[which(neighbors$sps_comp == comps[4])] * lmd4
        return(sum(c(nci1, nci2, nci3, nci4)))
      }
      
      # Define growth prediction function
      growth_pred <- function(nbhd_data, X0, Xb, gmax, pet_a, pet_b,
                              C, alpha, beta, lmd1, lmd2, lmd3, lmd4){
        ids <- unique(nbhd_data$tree_id)
        pred_grow <- rep(NA, times = length(ids))
        for(id in 1:length(ids)){
          neighbors <- nbhd_data %>%
            filter(tree_id == ids[id])
          pred_grow[id] <- gmax *
            exp((-0.5) * ((log(neighbors$dbh[1] / X0) / Xb) ^ 2)) *
            exp((-0.5) * (((neighbors$pet_dm[1] - pet_a) / pet_b) ^ 2)) *
            exp(-C * nci(neighbors, alpha, beta, lmd1, lmd2, lmd3, lmd4))
        }
        results <- data.frame(ids, pred_grow, stringsAsFactors = F)
        return(results)
      }
      
      # Predict test growth
      growth_ts <- growth_pred(ss_test, params["X0_opt"],
                               params["Xb_opt"], params["gmax_opt"],
                               params["pet_a_opt"], params["pet_b_opt"],
                               params["C_opt"], params["alpha_opt"],
                               params["beta_opt"], params["lmd1_opt"],
                               params["lmd2_opt"], params["lmd3_opt"],
                               params["lmd4_opt"])
      
    } else if(length(comps) == 5){
      
      # Define NCI function
      nci <- function(neighbors, alpha, beta, lmd1, lmd2, lmd3, lmd4, lmd5){
        raw <- (neighbors$dbh_comp ^ alpha) / (neighbors$prox ^ beta)
        nci1 <- raw[which(neighbors$sps_comp == comps[1])] * lmd1
        nci2 <- raw[which(neighbors$sps_comp == comps[2])] * lmd2
        nci3 <- raw[which(neighbors$sps_comp == comps[3])] * lmd3
        nci4 <- raw[which(neighbors$sps_comp == comps[4])] * lmd4
        nci5 <- raw[which(neighbors$sps_comp == comps[5])] * lmd5
        return(sum(c(nci1, nci2, nci3, nci4, nci5)))
      }
      
      # Define growth prediction function
      growth_pred <- function(nbhd_data, X0, Xb, gmax, pet_a, pet_b,
                              C, alpha, beta, lmd1, lmd2, lmd3, lmd4,
                              lmd5){
        ids <- unique(nbhd_data$tree_id)
        pred_grow <- rep(NA, times = length(ids))
        for(id in 1:length(ids)){
          neighbors <- nbhd_data %>%
            filter(tree_id == ids[id])
          pred_grow[id] <- gmax *
            exp((-0.5) * ((log(neighbors$dbh[1] / X0) / Xb) ^ 2)) *
            exp((-0.5) * (((neighbors$pet_dm[1] - pet_a) / pet_b) ^ 2)) *
            exp(-C * nci(neighbors, alpha, beta, lmd1, lmd2, lmd3, lmd4,
                         lmd5))
        }
        results <- data.frame(ids, pred_grow, stringsAsFactors = F)
        return(results)
      }
      
      # Predict test growth
      growth_ts <- growth_pred(ss_test, params["X0_opt"],
                               params["Xb_opt"], params["gmax_opt"],
                               params["pet_a_opt"], params["pet_b_opt"],
                               params["C_opt"], params["alpha_opt"],
                               params["beta_opt"], params["lmd1_opt"],
                               params["lmd2_opt"], params["lmd3_opt"],
                               params["lmd4_opt"], params["lmd5_opt"])
      
    } else if(length(comps) == 9){
      
      # Define NCI function
      nci <- function(neighbors, alpha, beta, lmd1, lmd2, lmd3, lmd4, lmd5,
                      lmd6, lmd7, lmd8, lmd9){
        raw <- (neighbors$dbh_comp ^ alpha) / (neighbors$prox ^ beta)
        nci1 <- raw[which(neighbors$sps_comp == comps[1])] * lmd1
        nci2 <- raw[which(neighbors$sps_comp == comps[2])] * lmd2
        nci3 <- raw[which(neighbors$sps_comp == comps[3])] * lmd3
        nci4 <- raw[which(neighbors$sps_comp == comps[4])] * lmd4
        nci5 <- raw[which(neighbors$sps_comp == comps[5])] * lmd5
        nci6 <- raw[which(neighbors$sps_comp == comps[6])] * lmd6
        nci7 <- raw[which(neighbors$sps_comp == comps[7])] * lmd7
        nci8 <- raw[which(neighbors$sps_comp == comps[8])] * lmd8
        nci9 <- raw[which(neighbors$sps_comp == comps[9])] * lmd9
        return(sum(c(nci1, nci2, nci3, nci4, nci5, nci6, nci7, nci8, nci9)))
      }
      
      # Create growth prediction function
      growth_pred <- function(nbhd_data, X0, Xb, gmax, pet_a, pet_b,
                              C, alpha, beta, lmd1, lmd2, lmd3, lmd4,
                              lmd5, lmd6, lmd7, lmd8, lmd9){
        ids <- unique(nbhd_data$tree_id)
        pred_grow <- rep(NA, times = length(ids))
        for(id in 1:length(ids)){
          neighbors <- nbhd_data %>%
            filter(tree_id == ids[id])
          pred_grow[id] <- gmax *
            exp((-0.5) * ((log(neighbors$dbh[1] / X0) / Xb) ^ 2)) *
            exp((-0.5) * (((neighbors$pet_dm[1] - pet_a) / pet_b) ^ 2)) *
            exp(-C * nci(neighbors, alpha, beta, lmd1, lmd2, lmd3, lmd4,
                         lmd5, lmd6, lmd7, lmd8, lmd9))
        }
        results <- data.frame(ids, pred_grow, stringsAsFactors = F)
        return(results)
      }
      
      # Predict test growth
      growth_ts <- growth_pred(ss_test, params["X0_opt"],
                               params["Xb_opt"], params["gmax_opt"],
                               params["pet_a_opt"], params["pet_b_opt"],
                               params["C_opt"], params["alpha_opt"],
                               params["beta_opt"], params["lmd1_opt"],
                               params["lmd2_opt"], params["lmd3_opt"],
                               params["lmd4_opt"], params["lmd5_opt"],
                               params["lmd6_opt"], params["lmd7_opt"],
                               params["lmd8_opt"], params["lmd9_opt"])
      
    }
  }
  
  # Combine predicted and observed growth
  test_obs_pred <- test_focals %>%
    left_join(growth_ts, by = c("tree_id" = "ids")) %>%
    rename(observations = annual_growth,
           predictions = pred_grow)
  
  # Calculate and store coefficient of determination
  best_mod$test_r2[sps] <- coef_det(test_obs_pred)
  
}


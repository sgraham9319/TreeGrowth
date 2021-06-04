
quantify_lkhd_fit <- function(method = "AIC"){
  
  # Define AICc function
  AICc_calc <- function(k, NLL, n){
    (2 * (k + NLL)) + ((2 * ((k^2) + k)) / (n - k - 1))
  }
  
  # Define coefficient of determination function
  coef_det <- function(x){
    1 - (sum((x$observations - x$predictions)^2) / 
           sum((x$observations - mean(x$observations))^2))
  }
  
  # Define focal species and training sets
  sps <- c("ABAM", "CANO", "PSME", "THPL", "TSHE", "TSME")
  train_sets <- 1:4
  
  # Create output data frame
  results <- matrix(NA, nrow = length(sps) * length(train_sets),
                    ncol = 3)
  
  # Load model comparison table
  if(method == "AIC"){
    lkhd_table <- read.csv("Figures/lkhd_model_selection_rand_train.csv",
                           stringsAsFactors = F)
  } else if(method == "cv"){
    # NEED TO CREATE A CV LKHD TABLE FIRST!
  }
  
  # Create table of best model structure for each species/training set
  best_strs <- lkhd_table %>%
    mutate(str = names(lkhd_table)[apply(lkhd_table[, 3:6],
                                         1, which.min) + 2]) %>%
    select(focal_sps, training_set, str)
  
  # Loop through training sets
  for(set in train_sets){
    
    # Load training data
    training <- read.csv(paste("Data/Output_data/rand_training", set,
                                 ".csv", sep = ""), stringsAsFactors = F)
    
    # Load test data
    test <- read.csv(paste("Data/Output_data/rand_test", set,
                             ".csv", sep = ""), stringsAsFactors = F)
    
    # Loop through focal species
    for(i in 1:length(sps)){
      
      # Subset to focal species and remove unneeded columns
      sing_sp <- training %>%
        arrange(tree_id) %>%
        filter(species == sps[i]) %>%
        select(tree_id, stand_id, species, dbh, prox, sps_comp, dbh_comp,
               annual_growth, pet_mm)
      ss_test <- test %>%
        arrange(tree_id) %>%
        filter(species == sps[i]) %>%
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
      ss_test <- ss_test %>%
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
        summarize(dbh = dbh[1], annual_growth = annual_growth[1],
                  pet_dm = pet_dm[1])
      test_focals <- ss_test %>%
        group_by(tree_id) %>%
        summarize(dbh = dbh[1], annual_growth = annual_growth[1],
                  pet_dm = pet_dm[1])
      
      # Extract best model structure
      mod <- best_strs %>%
        filter(focal_sps == sps[i] & training_set == set) %>%
        pull(str)
      
      # Load model output and calculate AICc for each model
      output <- read.csv(paste("Data/Output_data/", mod, "_rand", set, "_",
                                 sps[i], ".csv", sep = ""))
      
      # Calculate AICc for each model
      output$AICc <- AICc_calc(length(grep("_opt", names(output))),
                               output$NLL, nrow(focals))
      
      # Order output by AICc
      output <- output %>% arrange(AICc)
      
      # Define growth prediction functions
      if(mod == "no_comp"){
        
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
        
        # Predict training growth
        growth_tr <- growth_pred(sing_sp, output[1, "X0_opt"],
                                 output[1, "Xb_opt"], output[1, "gmax_opt"],
                                 output[1, "pet_a_opt"], output[1, "pet_b_opt"])
        
        # Predict test growth
        growth_ts <- growth_pred(ss_test, output[1, "X0_opt"],
                                 output[1, "Xb_opt"], output[1, "gmax_opt"],
                                 output[1, "pet_a_opt"], output[1, "pet_b_opt"])
        
      } else if(mod == "eq_comp"){
        
        # Define NCI function
        nci <- function(neighbors, alpha, beta){
          sum((neighbors$dbh_comp ^ alpha) / (neighbors$prox ^ beta))
        }
        
        # Define growth prediction function
        growth_pred <- function(nbhd_data, X0, Xb, gmax, pet_a, pet_b,
                                C, alpha, beta){
          ids <- unique(nbhd_data$tree_id)
          pred_grow <- rep(NA, times = length(ids))
          for(id in 1:length(ids)){
            neighbors <- nbhd_data %>%
              filter(tree_id == ids[id])
            pred_grow[id] <- gmax * 
              exp((-0.5) * (log(neighbors$dbh[1] / X0) / Xb) ^ 2) *
              exp((-0.5) * ((neighbors$pet_dm[1] - pet_a) / pet_b) ^ 2) *
              exp(-C * nci(neighbors, alpha, beta))
          }
          results <- data.frame(ids, pred_grow, stringsAsFactors = F)
          return(results)
        }
        
        # Predict training growth
        growth_tr <- growth_pred(sing_sp, output[1, "X0_opt"],
                                 output[1, "Xb_opt"], output[1, "gmax_opt"],
                                 output[1, "pet_a_opt"], output[1, "pet_b_opt"],
                                 output[1, "C_opt"], output[1, "alpha_opt"],
                                 output[1, "beta_opt"])
        
        # Predict test growth
        growth_ts <- growth_pred(ss_test, output[1, "X0_opt"],
                                 output[1, "Xb_opt"], output[1, "gmax_opt"],
                                 output[1, "pet_a_opt"], output[1, "pet_b_opt"],
                                 output[1, "C_opt"], output[1, "alpha_opt"],
                                 output[1, "beta_opt"])
        
      } else if(mod == "int_comp"){
        
        # Define NCI function
        nci <- function(neighbors, alpha, beta, intra, inter){
          raw <- (neighbors$dbh_comp ^ alpha) / (neighbors$prox ^ beta)
          cons <- which(neighbors$sps_comp == sps[i])
          hets <- which(neighbors$sps_comp != sps[i])
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
        
        # Predict training growth
        growth_tr <- growth_pred(sing_sp, output[1, "X0_opt"],
                                 output[1, "Xb_opt"], output[1, "gmax_opt"],
                                 output[1, "pet_a_opt"], output[1, "pet_b_opt"],
                                 output[1, "C_opt"], output[1, "alpha_opt"],
                                 output[1, "beta_opt"], output[1, "intra_opt"],
                                 output[1, "inter_opt"])
        
        # Predict test growth
        growth_ts <- growth_pred(ss_test, output[1, "X0_opt"],
                                 output[1, "Xb_opt"], output[1, "gmax_opt"],
                                 output[1, "pet_a_opt"], output[1, "pet_b_opt"],
                                 output[1, "C_opt"], output[1, "alpha_opt"],
                                 output[1, "beta_opt"], output[1, "intra_opt"],
                                 output[1, "inter_opt"])
        
      } else if(mod == "ss_comp"){
        
        # Convert rare competitor species to OTHR
        comm_comp <- read.csv("Data/Output_data/common_comps.csv",
                              stringsAsFactors = F)
        comm_comp <- comm_comp[, sps[i]]
        sing_sp <- sing_sp %>%
          mutate(sps_comp = if_else(sps_comp %in% comm_comp, sps_comp, "OTHR"))
        ss_test <- ss_test %>%
          mutate(sps_comp = if_else(sps_comp %in% comm_comp, sps_comp, "OTHR"))
        
        # Create vector of competitor species
        comps <- sort(unique(sing_sp$sps_comp))
        
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
          
          # Predict training growth
          growth_tr <- growth_pred(sing_sp, output[1, "X0_opt"],
                                   output[1, "Xb_opt"], output[1, "gmax_opt"],
                                   output[1, "pet_a_opt"], output[1, "pet_b_opt"],
                                   output[1, "C_opt"], output[1, "alpha_opt"],
                                   output[1, "beta_opt"], output[1, "lmd1_opt"],
                                   output[1, "lmd2_opt"], output[1, "lmd3_opt"],
                                   output[1, "lmd4_opt"])
          
          # Predict test growth
          growth_ts <- growth_pred(ss_test, output[1, "X0_opt"],
                                   output[1, "Xb_opt"], output[1, "gmax_opt"],
                                   output[1, "pet_a_opt"], output[1, "pet_b_opt"],
                                   output[1, "C_opt"], output[1, "alpha_opt"],
                                   output[1, "beta_opt"], output[1, "lmd1_opt"],
                                   output[1, "lmd2_opt"], output[1, "lmd3_opt"],
                                   output[1, "lmd4_opt"])
          
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
          
          # Predict training growth
          growth_tr <- growth_pred(sing_sp, output[1, "X0_opt"],
                                   output[1, "Xb_opt"], output[1, "gmax_opt"],
                                   output[1, "pet_a_opt"], output[1, "pet_b_opt"],
                                   output[1, "C_opt"], output[1, "alpha_opt"],
                                   output[1, "beta_opt"], output[1, "lmd1_opt"],
                                   output[1, "lmd2_opt"], output[1, "lmd3_opt"],
                                   output[1, "lmd4_opt"], output[1, "lmd5_opt"])
          
          # Predict test growth
          growth_ts <- growth_pred(ss_test, output[1, "X0_opt"],
                                   output[1, "Xb_opt"], output[1, "gmax_opt"],
                                   output[1, "pet_a_opt"], output[1, "pet_b_opt"],
                                   output[1, "C_opt"], output[1, "alpha_opt"],
                                   output[1, "beta_opt"], output[1, "lmd1_opt"],
                                   output[1, "lmd2_opt"], output[1, "lmd3_opt"],
                                   output[1, "lmd4_opt"], output[1, "lmd5_opt"])
          
        } else if(length(comps) == 8){
          
          # Define NCI function
          nci <- function(neighbors, alpha, beta, lmd1, lmd2, lmd3, lmd4, lmd5,
                          lmd6, lmd7, lmd8){
            raw <- (neighbors$dbh_comp ^ alpha) / (neighbors$prox ^ beta)
            nci1 <- raw[which(neighbors$sps_comp == comps[1])] * lmd1
            nci2 <- raw[which(neighbors$sps_comp == comps[2])] * lmd2
            nci3 <- raw[which(neighbors$sps_comp == comps[3])] * lmd3
            nci4 <- raw[which(neighbors$sps_comp == comps[4])] * lmd4
            nci5 <- raw[which(neighbors$sps_comp == comps[5])] * lmd5
            nci6 <- raw[which(neighbors$sps_comp == comps[6])] * lmd6
            nci7 <- raw[which(neighbors$sps_comp == comps[7])] * lmd7
            nci8 <- raw[which(neighbors$sps_comp == comps[8])] * lmd8
            return(sum(c(nci1, nci2, nci3, nci4, nci5, nci6, nci7, nci8)))
          }
          
          # Define growth prediction function
          growth_pred <- function(nbhd_data, X0, Xb, gmax, pet_a, pet_b,
                                  C, alpha, beta, lmd1, lmd2, lmd3, lmd4,
                                  lmd5, lmd6, lmd7, lmd8){
            ids <- unique(nbhd_data$tree_id)
            pred_grow <- rep(NA, times = length(ids))
            for(id in 1:length(ids)){
              neighbors <- nbhd_data %>%
                filter(tree_id == ids[id])
              pred_grow[id] <- gmax *
                exp((-0.5) * ((log(neighbors$dbh[1] / X0) / Xb) ^ 2)) *
                exp((-0.5) * (((neighbors$pet_dm[1] - pet_a) / pet_b) ^ 2)) *
                exp(-C * nci(neighbors, alpha, beta, lmd1, lmd2, lmd3, lmd4,
                             lmd5, lmd6, lmd7, lmd8))
            }
            results <- data.frame(ids, pred_grow, stringsAsFactors = F)
            return(results)
          }
          
          # Predict training growth
          growth_tr <- growth_pred(sing_sp, output[1, "X0_opt"],
                                   output[1, "Xb_opt"], output[1, "gmax_opt"],
                                   output[1, "pet_a_opt"], output[1, "pet_b_opt"],
                                   output[1, "C_opt"], output[1, "alpha_opt"],
                                   output[1, "beta_opt"], output[1, "lmd1_opt"],
                                   output[1, "lmd2_opt"], output[1, "lmd3_opt"],
                                   output[1, "lmd4_opt"], output[1, "lmd5_opt"],
                                   output[1, "lmd6_opt"], output[1, "lmd7_opt"],
                                   output[1, "lmd8_opt"])
          
          # Predict test growth
          growth_ts <- growth_pred(ss_test, output[1, "X0_opt"],
                                   output[1, "Xb_opt"], output[1, "gmax_opt"],
                                   output[1, "pet_a_opt"], output[1, "pet_b_opt"],
                                   output[1, "C_opt"], output[1, "alpha_opt"],
                                   output[1, "beta_opt"], output[1, "lmd1_opt"],
                                   output[1, "lmd2_opt"], output[1, "lmd3_opt"],
                                   output[1, "lmd4_opt"], output[1, "lmd5_opt"],
                                   output[1, "lmd6_opt"], output[1, "lmd7_opt"],
                                   output[1, "lmd8_opt"])
          
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
          
          # Predict training growth
          growth_tr <- growth_pred(sing_sp, output[1, "X0_opt"],
                                   output[1, "Xb_opt"], output[1, "gmax_opt"],
                                   output[1, "pet_a_opt"], output[1, "pet_b_opt"],
                                   output[1, "C_opt"], output[1, "alpha_opt"],
                                   output[1, "beta_opt"], output[1, "lmd1_opt"],
                                   output[1, "lmd2_opt"], output[1, "lmd3_opt"],
                                   output[1, "lmd4_opt"], output[1, "lmd5_opt"],
                                   output[1, "lmd6_opt"], output[1, "lmd7_opt"],
                                   output[1, "lmd8_opt"], output[1, "lmd9_opt"])
          
          # Predict test growth
          growth_ts <- growth_pred(ss_test, output[1, "X0_opt"],
                                   output[1, "Xb_opt"], output[1, "gmax_opt"],
                                   output[1, "pet_a_opt"], output[1, "pet_b_opt"],
                                   output[1, "C_opt"], output[1, "alpha_opt"],
                                   output[1, "beta_opt"], output[1, "lmd1_opt"],
                                   output[1, "lmd2_opt"], output[1, "lmd3_opt"],
                                   output[1, "lmd4_opt"], output[1, "lmd5_opt"],
                                   output[1, "lmd6_opt"], output[1, "lmd7_opt"],
                                   output[1, "lmd8_opt"], output[1, "lmd9_opt"])
          
        }
      }
      
      # Combine predicted and observed growth
      train_obs_pred <- focals %>%
        left_join(growth_tr, by = c("tree_id" = "ids")) %>%
        rename(observations = annual_growth,
               predictions = pred_grow)
      test_obs_pred <- test_focals %>%
        left_join(growth_ts, by = c("tree_id" = "ids")) %>%
        rename(observations = annual_growth,
               predictions = pred_grow)
      
      # Add coefficients of determination to results table
      results[4 * (i - 1) + set, 2] <- coef_det(train_obs_pred)
      results[4 * (i - 1) + set, 3] <- coef_det(test_obs_pred)
      
      # Add number of focal trees in training set to results table
      results[4 * (i - 1) + set, 1] <- nrow(focals)
      
    }
  }
  
  # Format table
  results <- as.data.frame(results)
  results <- results %>%
    mutate(focal_sps = rep(sps, each = length(train_sets)),
           training_set = rep(train_sets, times = length(sps))) %>%
    rename(sample_size = "V1", train_lkhd = "V2", test_lkhd = "V3") %>%
    select(focal_sps, training_set, sample_size, train_lkhd, test_lkhd)
  
  # Return results table
  return(results)
}

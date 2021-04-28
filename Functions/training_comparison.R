

training_comparison <- function(focal_sps, model_strs, sets){
  
  # Define AICc function
  AICc_calc <- function(k, NLL, n){
    (2 * (k + NLL)) + ((2 * ((k^2) + k)) / (n - k - 1))
  }
  
  # Create empty data frame to store best models
  best_models <- data.frame(
    "gmax_opt" = double(),
    "X0_opt" = double(),
    "Xb_opt" = double(),
    "pet_a_opt" = double(),
    "pet_b_opt" = double())
  
  for(train in sets){
    
    # Load training data
    training <- read.csv(paste("Data/Output_data/training", train, ".csv",
                               sep = ""), stringsAsFactors = F)
    
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
      summarize(dbh = dbh[1], annual_growth = annual_growth[1],
                pet_dm = pet_dm[1])
    
    for(mod in model_strs){
      
      # Load model output
      output <- read.csv(paste("Data/Output_data/", mod, train, "_", focal_sps,
                               ".csv", sep = ""))
      
      # Calculate AICc for each model
      output$AICc <- AICc_calc(length(grep("_opt", names(output))),
                               output$NLL, nrow(focals))
      
      # Calculate dAICc for each model
      output$dAICc <- output$AICc - min(output$AICc)
      
      # Order output by dAICc
      output <- output %>% arrange(dAICc)
      
      # Get best model of this structure
      top_mod <- output[1, c(grep("_opt", names(output)),
                             which(names(output) == "AICc"))]
      
      # Add mean pet and max dbh for this training set
      top_mod$mean_pet <- mean(focals$pet_dm)
      top_mod$max_dbh <- max(focals$dbh)
      
      # Calculate mean NCI for each model run (if NCI included in model)
      if(length(grep("eq", mod)) == 1){
        
        # Define NCI function
        nci <- function(neighbors, alpha, beta){
          sum((neighbors$dbh_comp ^ alpha) / (neighbors$prox ^ beta))
        }
        
        # Create function to calculate nci for each individual
        nci_test <- function(nbhd_data, alpha, beta){
          ids <- unique(nbhd_data$tree_id)
          nci_vals <- rep(NA, times = length(ids))
          for(id in 1:length(ids)){
            neighbors <- nbhd_data[nbhd_data$tree_id == ids[id], ]
            nci_vals[id] <- nci(neighbors, alpha, beta)
          }
          nci_vals
        }
        
        # Calculate mean NCI using optimized alpha and beta for each model
        top_mod$mean_nci <- mean(nci_test(sing_sp, top_mod$alpha_opt,
                                          top_mod$beta_opt))
        
      } else if(length(grep("int", mod)) == 1){
        
        # Define NCI function
        nci <- function(neighbors, alpha, beta, intra, inter){
          raw <- (neighbors$dbh_comp ^ alpha) / (neighbors$prox ^ beta)
          cons <- which(neighbors$sps_comp == focal_sps)
          hets <- which(neighbors$sps_comp != focal_sps)
          nci_con <- raw[cons] * intra
          nci_het <- raw[hets] * inter
          return(sum(c(nci_con, nci_het)))
        }
        
        # Create function to calculate nci for each individual
        nci_test <- function(nbhd_data, alpha, beta, intra, inter){
          ids <- unique(nbhd_data$tree_id)
          nci_vals <- rep(NA, times = length(ids))
          for(id in 1:length(ids)){
            neighbors <- nbhd_data[nbhd_data$tree_id == ids[id], ]
            nci_vals[id] <- nci(neighbors, alpha, beta, intra, inter)
          }
          nci_vals
        }
        
        # Calculate mean NCI using optimized alpha and beta for each model
        top_mod$mean_nci <- mean(nci_test(sing_sp, top_mod$alpha_opt,
                                          top_mod$beta_opt,
                                          top_mod$intra_opt,
                                          top_mod$inter_opt))
      }
      
      # Add model to collection of best models
      best_models <- bind_rows(best_models, top_mod)
      
    }
  }
  
  # Add training set and model structure columns
  best_models$training <- rep(1:4, each = length(model_strs))
  best_models$mod_str <- rep(model_strs, times = length(sets))
  
  
  # Get best overall model per training set
  best_model_per_training <- best_models %>%
    arrange(AICc) %>%
    group_by(training) %>%
    summarize(model_str = mod_str[1],
              gmax = gmax_opt[1],
              X0 = X0_opt[1],
              Xb = Xb_opt[1],
              mean_pet = mean_pet[1],
              pet_a = pet_a_opt[1],
              pet_b = pet_b_opt[1],
              C = C_opt[1],
              mean_nci = mean_nci[1],
              max_dbh = max_dbh[1])
  
  # Call size effect function on combined data
  comp_plot <- size_effect_comp(best_model_per_training, focal_sps,
                                cols = c("red", "blue", "gold", "green"))
  
  # Combine result
  result <- list(best_models = best_models, comp_plot = comp_plot)
  
  # Return result
  return(result)
  
}
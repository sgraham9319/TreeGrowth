

lkhd_model_select <- function(training_type = "regular"){
  
  # Define AICc function
  AICc_calc <- function(k, NLL, n){
    (2 * (k + NLL)) + ((2 * ((k^2) + k)) / (n - k - 1))
  }
  
  # Create vectors of focal species and training sets
  focal_sps <- c("ABAM", "CANO", "PSME", "THPL", "TSHE", "TSME")
  train_sets <- 1:4
  
  # Create vector of model structures
  model_strs <- c("no_comp", "eq_comp", "int_comp", "ss_comp")
  
  # Create matrix to store AIC output
  aic_table <- matrix(NA, ncol = length(model_strs),
                      nrow = length(focal_sps) * length(train_sets))
  
  # Loop through training sets
  for(set in train_sets){
    
    # Load training data
    if(training_type == "regular"){
      training <- read.csv(paste("Data/Output_data/training", set, ".csv",
                                 sep = ""), stringsAsFactors = F)
    } else if(training_type == "random"){
      training <- read.csv(paste("Data/Output_data/rand_training", set, ".csv",
                                 sep = ""), stringsAsFactors = F)
    }
    
    # Loop through focal species
    for(sps in 1:length(focal_sps)){
      
      # Subset training to focal species
      sing_sp <- training %>%
        filter(species == focal_sps[sps])
      
      # Calculate number of focals
      num_focals <- length(unique(sing_sp$tree_id))
      
      # Create empty vector to store AICc values
      aic_vals <- rep(NA, times = 4)
      
      # Loop through model structures
      for(i in 1:length(model_strs)){
        
        # Load fitted models
        if(training_type == "regular"){
          output <- read.csv(paste("Data/Output_data/", model_strs[i],
                                   set, "_", focal_sps[sps], ".csv", sep = ""))
        } else if(training_type == "random"){
          output <- read.csv(paste("Data/Output_data/", model_strs[i], "_rand",
                                   set, "_", focal_sps[sps], ".csv", sep = ""))
        }
        
        # Calculate AICc for each model
        output$AICc <- AICc_calc(length(grep("_opt", names(output))),
                                 output$NLL, num_focals)
        
        # Add best model AICc to AICc vector
        aic_vals[i] <- min(output$AICc)
      }
      
      # Convert AICc vector to dAICc
      aic_vals <- round(aic_vals - min(aic_vals), 2)
      
      # Add to output table
      aic_table[set + ((sps - 1) * 4), ] <- aic_vals
      
    }
  }
  
  # Convert AIC table to data frame
  aic_table <- as.data.frame(aic_table)
  
  # Add focal species and training set columns
  aic_table$focal_sps <- rep(focal_sps, each = length(train_sets))
  aic_table$training_set <- rep(train_sets, times = length(focal_sps))
  
  # Change column names
  names(aic_table)[1:4] <- model_strs
  
  # Re-order columns
  aic_table <- aic_table[, c(5, 6, 1:4)]
  
  # Return AIC table
  return(aic_table)
  
}

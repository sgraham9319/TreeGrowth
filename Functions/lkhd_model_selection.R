

lkhd_model_select <- function(method = "AIC", num_train_sets = 4){
  
  # Define AICc function
  AICc_calc <- function(k, NLL, n){
    (2 * (k + NLL)) + ((2 * ((k^2) + k)) / (n - k - 1))
  }
  
  # Create vectors of focal species and training sets
  focal_sps <- c("ABAM", "CANO", "PSME", "THPL", "TSHE", "TSME")
  train_sets <- 1:num_train_sets
  
  # Create vector of model structures
  model_strs <- c("no_comp", "eq_comp", "int_comp", "ss_comp")
  
  # Create matrix to store AIC output
  comp_table <- matrix(NA, ncol = length(model_strs),
                      nrow = length(focal_sps) * length(train_sets))
  
  # Loop through training sets
  for(set in train_sets){
    
    # Load training data
    training <- read.csv(paste("Data/Output_data/rand_training", set, ".csv",
                                 sep = ""), stringsAsFactors = F)
    
    # Loop through focal species
    for(sps in 1:length(focal_sps)){
      
      # Subset training to focal species
      sing_sp <- training %>%
        filter(species == focal_sps[sps])
      
      # Calculate number of focals
      num_focals <- length(unique(sing_sp$tree_id))
      
      # Create empty vector to store model performance values
      perf_vals <- rep(NA, times = 4)
      
      # Loop through model structures
      for(i in 1:length(model_strs)){
        
        # Load fitted models
        if(method == "AIC"){
          output <- read.csv(paste("Data/Output_data/", model_strs[i], "_rand",
                                   set, "_", focal_sps[sps], ".csv", sep = ""))
        } else if(method == "cv"){
          output <- read.csv(paste("Data/Output_data/", model_strs[i], "_cv",
                                   set, "_", focal_sps[sps], ".csv", sep = ""))
        }
        
        # Calculate performance metric for associated model type
        if(method == "AIC"){
          
          # Calculate AICc for each model
          output$AICc <- AICc_calc(length(grep("_opt", names(output))),
                                   output$NLL, num_focals)
          
          # Add best model AICc to performance metric vector
          perf_vals[i] <- min(output$AICc)
          
        } else if(method == "cv"){
          
          # Extract minimum mse for each cv set
          min_mse <- output %>%
            arrange(mse) %>%
            group_by(cv_set) %>%
            filter(row_number() == 1)
          
          # Calculate and store mean square error across cross-validation sets
          perf_vals[i] <- mean(min_mse$mse)
          
        }
      }
      
      # If using AIC models convert AICc to dAICc
      if(method == "AIC"){
        perf_vals <- round(perf_vals - min(perf_vals), 2)
      }
      
      # Add to output table
      comp_table[set + ((sps - 1) * length(train_sets)), ] <- perf_vals
      
    }
  }
  
  # Convert comparison table to data frame
  comp_table <- as.data.frame(comp_table)
  
  # Add focal species and training set columns
  comp_table$focal_sps <- rep(focal_sps, each = length(train_sets))
  comp_table$training_set <- rep(train_sets, times = length(focal_sps))
  
  # Change column names
  names(comp_table)[1:4] <- model_strs
  
  # Re-order columns
  comp_table <- comp_table[, c(5, 6, 1:4)]
  
  # Return comparison table
  return(comp_table)
  
}

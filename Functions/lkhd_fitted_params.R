

lkhd_fitted_params <- function(){
  
  # Load model selection data
  lkhd_table <- read.csv("Figures/lkhd_model_selection_rand_train.csv",
                         stringsAsFactors = F)
  
  # Get best model structure for each species/training set
  best_strs <- lkhd_table %>%
    mutate(str = names(lkhd_table)[apply(lkhd_table[, 3:6],
                                         1, which.min) + 2]) %>%
    select(focal_sps, training_set, str)
  
  # Create empty table to store fitted parameter values
  param_table <- data.frame(X0_opt = double())
  
  # Define names of parameters
  par_names <- paste(c("X0", "Xb", "gmax", "pet_a", "pet_b", "C",
                       "alpha", "beta", "intra", "inter", "lmd1",
                       "lmd2", "lmd3", "lmd4", "lmd5", "lmd6", "lmd7",
                       "lmd8", "lmd9", "sigma"), "opt", sep = "_")
  
  # Loop through best models table extracting parameter values
  for(i in 1:nrow(best_strs)){
    
    # Load fitted model
    output <- read.csv(paste("Data/Output_data/", best_strs$str[i], "_rand",
                             best_strs$training_set[i], "_",
                             best_strs$focal_sps[i], ".csv", sep = ""))
    
    # Order by NLL
    output <- output %>%
      arrange(NLL)
    
    # Add best model params to overall table
    param_table <- param_table %>%
      bind_rows(output[1, which(names(output) %in% par_names)])
  }
  
  # Combine with best model identification information
  results <- cbind(best_strs, param_table)
  
}

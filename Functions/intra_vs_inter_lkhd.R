
intra_vs_inter_lkhd <- function(train_type){
  
  # Load AICc function
  AICc_calc <- function(k, NLL, n){
    (2 * (k + NLL)) + ((2 * ((k^2) + k)) / (n - k - 1))
  }
  
  # Create vectors of training sets and species
  train_sets <- 1:4
  sps <- c("ABAM", "CANO", "PSME", "THPL", "TSHE", "TSME")
  
  # Create empty results table
  intra <- matrix(NA, nrow = length(sps), ncol = length(train_sets))
  
  # Load table containing sample size information
  if(train_type == "regular"){
    fit <- read.csv("Figures/train_test_fit.csv")
  } else if(train_type == "random"){
    fit <- read.csv("Figures/train_test_fit_rand.csv")
  }
  
  # Extract sample size information
  samp_size <- fit %>%
    select(focal_sps, training_set, sample_size) %>%
    pivot_wider(names_from = training_set, values_from = sample_size)
  
  # Loop through training sets and species
  for(set in train_sets){
    for(i in 1:length(sps)){
      
      # Load int comp model output
      if(train_type == "regular"){
        output <- read.csv(paste("Data/Output_data/int_comp", set, "_",
                                 sps[i], ".csv", sep = ""))
      } else if(train_type == "random"){
        output <- read.csv(paste("Data/Output_data/int_comp_rand", set, "_",
                                 sps[i], ".csv", sep = ""))
      }
      
      # Calculate AICc
      output$AICc <- AICc_calc(length(grep("_opt", names(output))), output$NLL,
                               pull(samp_size[i, set + 1]))
      
      # Arrange output by AICc
      output <- output %>%
        arrange(AICc)
      
      # Determine whether intraspecific competition strongest
      if(output$intra_opt[1] > output$inter_opt[1]){
        intra[i, set] <- "-"
      } else{
        intra[i, set] <- "+"
      }
    }
  }
  
  # Format output table
  intra <- as.data.frame(intra)
  intra <- intra %>%
    rename(L1 = V1, L2 = V2, L3 = V3, L4 = V4) %>%
    mutate(species = sps) %>%
    select(species, L1, L2, L3, L4)
  
  # Return results
  return(intra)
}
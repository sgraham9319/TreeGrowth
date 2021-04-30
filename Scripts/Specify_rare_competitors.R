library(tidyverse)

# Load training sets
train1 <- read.csv("Data/Output_data/training1.csv", stringsAsFactors = F)
train2 <- read.csv("Data/Output_data/training2.csv", stringsAsFactors = F)
train3 <- read.csv("Data/Output_data/training3.csv", stringsAsFactors = F)
train4 <- read.csv("Data/Output_data/training4.csv", stringsAsFactors = F)

# Create vector of species
species <- c("ABAM", "CANO", "PSME", "THPL", "TSHE", "TSME")

# Create list to store final output
output <- vector(mode = "list", length = length(species))
names(output) <- species

# Loop through species adding rare species to output list
for(focal_sps in species){
  
  # Create output table for this species
  sing_sp <- data.frame(
    "ABAM" = integer())
  
  for(i in 1:4){
    
    # Subset to focal species
    training <- get(paste("train", i, sep = "")) %>%
      filter(species == focal_sps)
    
    # Extract competitor species numbers as data frame
    comp_freq <- as.data.frame(t(as.matrix(table(training$sps_comp))))
    
    # Add to output data frame
    sing_sp <- sing_sp %>%
      bind_rows(comp_freq)
  }
  
  # Calculate mean number interactions for each competitor species
  mean_ints <- apply(sing_sp, 2, mean, na.rm = T)
  
  # Create vector of rare competitors
  rare_sps <- names(which(mean_ints < 100))
  
  # Add to output list
  output[[focal_sps]] <- rare_sps
  
}

# Convert output list to data frame
output <- sapply(output, '[', seq(max(lengths(output))))

# Save to csv
write.csv(output, "Data/Output_data/rare_comps.csv", row.names = F)

# Load functions
source("Functions/nci.R")

# Specify training set
#set <- 1

# Load fitted parameters
fits <- read.csv("Figures/lkhd_AIC_param_table.csv")

# Create empty data frame to store average DBH, PET and NCI values
mean_vals <- data.frame(
  mean_dbh = double(),
  mean_pet = double(),
  mean_nci = double(),
  max_dbh = double(),
  min_pet = double(),
  max_pet = double()
)

# Loop through training sets
for(set in 1:4){
  
  # Extract fits for current training set
  set_fits <- fits %>%
    filter(training_set == set) %>%
    arrange(focal_sps)
  
  # Create empty columns for DBH, PET, and NCI
  set_fits <- set_fits %>%
    mutate(mean_dbh = NA,
           mean_pet = NA,
           mean_nci = NA,
           max_dbh = NA,
           min_pet = NA,
           max_pet = NA)
  
  # Load training data
  train <- read.csv(paste("Data/Output_data/rand_training", set, ".csv",
                          sep = ""))
  
  # Loop through focal species calculating average DBH, PET and NCI
  for(sps in 1:nrow(set_fits)){
    
    # Define focal species as character
    focal_sps <- set_fits$focal_sps[sps]
    
    # Subset training to focal species and remove unneeded columns
    sing_sp <- train %>%
      arrange(tree_id) %>%
      filter(species == focal_sps) %>%
      select(tree_id, stand_id, species, dbh, prox, sps_comp, dbh_comp,
             annual_growth, pet_mm)
    
    # Change units of variables to match fitting procedure
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
    
    # Calculate and store average and max DBH
    set_fits$mean_dbh[sps] <- mean(focals$dbh)
    set_fits$max_dbh[sps] <- max(focals$dbh)
    
    # Calculate and store average, min and max PET
    set_fits$mean_pet[sps] <- mean(focals$pet_dm)
    set_fits$min_pet[sps] <- min(focals$pet_dm)
    set_fits$max_pet[sps] <- max(focals$pet_dm)
    
    # Calculate and store average NCI
    if(set_fits$str[sps] == "eq_comp"){
      
      # Create function to calculate NCI for each individual
      nci_test <- function(nbhd_data, alpha, beta){
        ids <- unique(nbhd_data$tree_id)
        nci_vals <- rep(NA, times = length(ids))
        for(id in 1:length(ids)){
          neighbors <- nbhd_data[nbhd_data$tree_id == ids[id], ]
          nci_vals[id] <- nci_eq(neighbors, alpha, beta)
        }
        nci_vals
      }
      
      # Calculate and store average NCI
      set_fits$mean_nci[sps] <- mean(nci_test(sing_sp,
                                              set_fits[sps, "alpha_opt"],
                                              set_fits[sps, "beta_opt"]))
      
    } else if(set_fits$str[sps] == "int_comp"){
      
      # Create function to calculate NCI for each individual
      nci_test <- function(nbhd_data, alpha, beta, intra, inter){
        ids <- unique(nbhd_data$tree_id)
        nci_vals <- rep(NA, times = length(ids))
        for(id in 1:length(ids)){
          neighbors <- nbhd_data[nbhd_data$tree_id == ids[id], ]
          nci_vals[id] <- nci_int(neighbors, alpha, beta, intra, inter)
        }
        nci_vals
      }
      
      # Calculate and store average NCI
      set_fits$mean_nci[sps] <- mean(nci_test(sing_sp,
                                              set_fits[sps, "alpha_opt"],
                                              set_fits[sps, "beta_opt"],
                                              set_fits[sps, "intra_opt"],
                                              set_fits[sps, "inter_opt"]))
      
    } else if(set_fits$str[sps] == "ss_comp"){
      
      # Load table of common competitors
      comm_comp <- read.csv("Data/Output_data/common_comps.csv",
                            stringsAsFactors = F)
      
      # Extract common competitors for focal species
      comm_comp <- comm_comp[, focal_sps]
      
      # Convert rare competitor species to OTHR
      sing_sp <- sing_sp %>%
        mutate(sps_comp = if_else(sps_comp %in% comm_comp, sps_comp, "OTHR"))
      
      # Create vector of competitor species
      comps <- sort(unique(sing_sp$sps_comp))
      
      # Define functions based on number of competitors
      if(length(comps) == 4){
        
        # Create function to calculate NCI for each individual
        nci_test <- function(nbhd_data, alpha, beta, lmd1, lmd2, lmd3, lmd4){
          ids <- unique(nbhd_data$tree_id)
          nci_vals <- rep(NA, times = length(ids))
          for(id in 1:length(ids)){
            neighbors <- nbhd_data[nbhd_data$tree_id == ids[id], ]
            nci_vals[id] <- nci_ss4(neighbors, alpha, beta, lmd1, lmd2, lmd3,
                                    lmd4)
          }
          nci_vals
        }
        
        # Calculate and store average NCI
        set_fits$mean_nci[sps] <- mean(nci_test(sing_sp,
                                                set_fits[sps, "alpha_opt"],
                                                set_fits[sps, "beta_opt"],
                                                set_fits[sps, "lmd1_opt"],
                                                set_fits[sps, "lmd2_opt"],
                                                set_fits[sps, "lmd3_opt"],
                                                set_fits[sps, "lmd4_opt"]))
        
      } else if(length(comps) == 5){
        
        # Create function to calculate NCI for each individual
        nci_test <- function(nbhd_data, alpha, beta, lmd1, lmd2, lmd3,
                             lmd4, lmd5){
          ids <- unique(nbhd_data$tree_id)
          nci_vals <- rep(NA, times = length(ids))
          for(id in 1:length(ids)){
            neighbors <- nbhd_data[nbhd_data$tree_id == ids[id], ]
            nci_vals[id] <- nci_ss5(neighbors, alpha, beta, lmd1, lmd2, lmd3,
                                    lmd4, lmd5)
          }
          nci_vals
        }
        
        # Calculate and store average NCI
        set_fits$mean_nci[sps] <- mean(nci_test(sing_sp,
                                                set_fits[sps, "alpha_opt"],
                                                set_fits[sps, "beta_opt"],
                                                set_fits[sps, "lmd1_opt"],
                                                set_fits[sps, "lmd2_opt"],
                                                set_fits[sps, "lmd3_opt"],
                                                set_fits[sps, "lmd4_opt"],
                                                set_fits[sps, "lmd5_opt"]))
        
      } else if(length(comps) == 8){
        
        # Create function to calculate NCI for each individual
        nci_test <- function(nbhd_data, alpha, beta, lmd1, lmd2, lmd3,
                             lmd4, lmd5, lmd6, lmd7, lmd8){
          ids <- unique(nbhd_data$tree_id)
          nci_vals <- rep(NA, times = length(ids))
          for(id in 1:length(ids)){
            neighbors <- nbhd_data[nbhd_data$tree_id == ids[id], ]
            nci_vals[id] <- nci_ss8(neighbors, alpha, beta, lmd1, lmd2, lmd3,
                                    lmd4, lmd5, lmd6, lmd7, lmd8)
          }
          nci_vals
        }
        
        # Calculate and store average NCI
        set_fits$mean_nci[sps] <- mean(nci_test(sing_sp,
                                                set_fits[sps, "alpha_opt"],
                                                set_fits[sps, "beta_opt"],
                                                set_fits[sps, "lmd1_opt"],
                                                set_fits[sps, "lmd2_opt"],
                                                set_fits[sps, "lmd3_opt"],
                                                set_fits[sps, "lmd4_opt"],
                                                set_fits[sps, "lmd5_opt"],
                                                set_fits[sps, "lmd6_opt"],
                                                set_fits[sps, "lmd7_opt"],
                                                set_fits[sps, "lmd8_opt"]))
        
      } else if(length(comps) == 9){
        
        # Create function to calculate NCI for each individual
        nci_test <- function(nbhd_data, alpha, beta, lmd1, lmd2, lmd3,
                             lmd4, lmd5, lmd6, lmd7, lmd8, lmd9){
          ids <- unique(nbhd_data$tree_id)
          nci_vals <- rep(NA, times = length(ids))
          for(id in 1:length(ids)){
            neighbors <- nbhd_data[nbhd_data$tree_id == ids[id], ]
            nci_vals[id] <- nci_ss9(neighbors, alpha, beta, lmd1, lmd2, lmd3,
                                    lmd4, lmd5, lmd6, lmd7, lmd8, lmd9)
          }
          nci_vals
        }
        
        # Calculate and store average NCI
        set_fits$mean_nci[sps] <- mean(nci_test(sing_sp,
                                                set_fits[sps, "alpha_opt"],
                                                set_fits[sps, "beta_opt"],
                                                set_fits[sps, "lmd1_opt"],
                                                set_fits[sps, "lmd2_opt"],
                                                set_fits[sps, "lmd3_opt"],
                                                set_fits[sps, "lmd4_opt"],
                                                set_fits[sps, "lmd5_opt"],
                                                set_fits[sps, "lmd6_opt"],
                                                set_fits[sps, "lmd7_opt"],
                                                set_fits[sps, "lmd8_opt"],
                                                set_fits[sps, "lmd9_opt"]))
        
      }
    }
  }
  
  # Store mean DBH, PET and NCI for this training set
  mean_vals <- mean_vals %>%
    bind_rows(set_fits %>% select(mean_dbh, mean_pet, mean_nci, max_dbh,
                                  min_pet, max_pet))
  
}

# Add focal species and training set columns to mean values table
mean_vals <- mean_vals %>%
  mutate(focal_sps = rep(c("ABAM", "CANO", "PSME", "THPL", "TSHE", "TSME"), 
                         times = 4),
         training_set = rep(1:4, each = 6))

# Join mean values table to overall fits table
fits <- fits %>%
  left_join(mean_vals, by = c("focal_sps", "training_set"))

# Define size effect plotting function
size_eff <- function(x){
  gmax * exp((-0.5) * ((log(x / X0) / Xb) ^ 2)) *
    exp((-0.5) * (((pet_mean - pet_a) / pet_b) ^ 2)) *
    exp(-C * nci_mean)
}

# Define PET effect plotting function
pet_eff <- function(x){
  gmax * exp((-0.5) * ((log(dbh_mean / X0) / Xb) ^ 2)) *
    exp((-0.5) * (((x - pet_a) / pet_b) ^ 2)) *
    exp(-C * nci_mean)
}

# Define colors to use for different species
cols <- c("red", "pink", "light blue", "light green", "dark blue", "orange")

# Testing
effect <- "pet"


# Adjust plotting window parameters
par(mfcol = c(2, 2))
par(mar = c(0,0,0,0))
par(oma = c(4,4,0.1,0.1))
par(mgp = c(3, 0.3, 0))

# Create first empty plot
if(effect == "size"){
  plot(NULL, ylim = c(0, 3), xlim = c(0, max(fits$max_dbh)),
       yaxt = "n", xaxt = "n", xlab = "", ylab = "")
  text(x = 12.5, y = 2.8, label = "Training 1", font = 2)
} else if(effect == "pet"){
  plot(NULL, ylim = c(0, 3), xlim = c(min(fits$min_pet), max(fits$max_pet)),
       yaxt = "n", xaxt = "n", xlab = "", ylab = "")
  text(x = 4, y = 2.8, label = "Training 1", font = 2)
}

# Extract fits for this training set
sub_fits <- fits %>%
  filter(training_set == 1) %>%
  arrange(focal_sps)

# Add curve for each species
for(i in 1:nrow(sub_fits)){
  
  # Define parameters
  gmax <- sub_fits[i, "gmax_opt"]
  X0 <- sub_fits[i, "X0_opt"]
  Xb <- sub_fits[i, "Xb_opt"]
  pet_a <- sub_fits[i, "pet_a_opt"]
  pet_b <- sub_fits[i, "pet_b_opt"]
  C <- sub_fits[i, "C_opt"]
  pet_mean <- sub_fits[i, "mean_pet"]
  dbh_mean <- sub_fits[i, "mean_dbh"]
  nci_mean <- sub_fits[i, "mean_nci"]
  max_dbh <- sub_fits[i, "max_dbh"]
  min_pet <- sub_fits[i, "min_pet"]
  max_pet <- sub_fits[i, "max_pet"]
  
  # Add curve
  if(effect == "size"){
    curve(size_eff(x), from = 0.1, to = max_dbh, add = T, col = cols[i])
  } else if(effect == "pet"){
    curve(pet_eff(x), from = min_pet, to = max_pet,
          add = T, col = cols[i])
  }
}

# Create second empty plot
if(effect == "size"){
  plot(NULL, ylim = c(0, 3), xlim = c(0, max(fits$max_dbh)),
       yaxt = "n", xaxt = "n", xlab = "", ylab = "")
  text(x = 12.5, y = 2.8, label = "Training 2", font = 2)
} else if(effect == "pet"){
  plot(NULL, ylim = c(0, 3), xlim = c(min(fits$min_pet), max(fits$max_pet)),
       yaxt = "n", xaxt = "n", xlab = "", ylab = "")
  text(x = 4, y = 2.8, label = "Training 2", font = 2)
}

# Extract fits for this training set
sub_fits <- fits %>%
  filter(training_set == 2) %>%
  arrange(focal_sps)

# Add curve for each species
for(i in 1:nrow(sub_fits)){
  
  # Define parameters
  gmax <- sub_fits[i, "gmax_opt"]
  X0 <- sub_fits[i, "X0_opt"]
  Xb <- sub_fits[i, "Xb_opt"]
  pet_a <- sub_fits[i, "pet_a_opt"]
  pet_b <- sub_fits[i, "pet_b_opt"]
  C <- sub_fits[i, "C_opt"]
  pet_mean <- sub_fits[i, "mean_pet"]
  dbh_mean <- sub_fits[i, "mean_dbh"]
  nci_mean <- sub_fits[i, "mean_nci"]
  max_dbh <- sub_fits[i, "max_dbh"]
  min_pet <- sub_fits[i, "min_pet"]
  max_pet <- sub_fits[i, "max_pet"]
  
  # Add curve
  if(effect == "size"){
    curve(size_eff(x), from = 0.1, to = max_dbh, add = T, col = cols[i])
  } else if(effect == "pet"){
    curve(pet_eff(x), from = min_pet, to = max_pet,
          add = T, col = cols[i])
  }
}

# Create third empty plot
if(effect == "size"){
  plot(NULL, ylim = c(0, 3), xlim = c(0, max(fits$max_dbh)),
       yaxt = "n", xaxt = "n", xlab = "", ylab = "")
  text(x = 12.5, y = 2.8, label = "Training 3", font = 2)
} else if(effect == "pet"){
  plot(NULL, ylim = c(0, 3), xlim = c(min(fits$min_pet), max(fits$max_pet)),
       yaxt = "n", xaxt = "n", xlab = "", ylab = "")
  text(x = 4, y = 2.8, label = "Training 3", font = 2)
}

# Extract fits for this training set
sub_fits <- fits %>%
  filter(training_set == 3) %>%
  arrange(focal_sps)

# Add curve for each species
for(i in 1:nrow(sub_fits)){
  
  # Define parameters
  gmax <- sub_fits[i, "gmax_opt"]
  X0 <- sub_fits[i, "X0_opt"]
  Xb <- sub_fits[i, "Xb_opt"]
  pet_a <- sub_fits[i, "pet_a_opt"]
  pet_b <- sub_fits[i, "pet_b_opt"]
  C <- sub_fits[i, "C_opt"]
  pet_mean <- sub_fits[i, "mean_pet"]
  dbh_mean <- sub_fits[i, "mean_dbh"]
  nci_mean <- sub_fits[i, "mean_nci"]
  max_dbh <- sub_fits[i, "max_dbh"]
  min_pet <- sub_fits[i, "min_pet"]
  max_pet <- sub_fits[i, "max_pet"]
  
  # Add curve
  if(effect == "size"){
    curve(size_eff(x), from = 0.1, to = max_dbh, add = T, col = cols[i])
  } else if(effect == "pet"){
    curve(pet_eff(x), from = min_pet, to = max_pet,
          add = T, col = cols[i])
  }
}

# Create fourth empty plot
if(effect == "size"){
  plot(NULL, ylim = c(0, 3), xlim = c(0, max(fits$max_dbh)),
       yaxt = "n", xaxt = "n", xlab = "", ylab = "")
  text(x = 12.5, y = 2.8, label = "Training 4", font = 2)
} else if(effect == "pet"){
  plot(NULL, ylim = c(0, 3), xlim = c(min(fits$min_pet), max(fits$max_pet)),
       yaxt = "n", xaxt = "n", xlab = "", ylab = "")
  text(x = 4, y = 2.8, label = "Training 4", font = 2)
}

# Extract fits for this training set
sub_fits <- fits %>%
  filter(training_set == 4) %>%
  arrange(focal_sps)

# Add curve for each species
for(i in 1:nrow(sub_fits)){
  
  # Define parameters
  gmax <- sub_fits[i, "gmax_opt"]
  X0 <- sub_fits[i, "X0_opt"]
  Xb <- sub_fits[i, "Xb_opt"]
  pet_a <- sub_fits[i, "pet_a_opt"]
  pet_b <- sub_fits[i, "pet_b_opt"]
  C <- sub_fits[i, "C_opt"]
  pet_mean <- sub_fits[i, "mean_pet"]
  dbh_mean <- sub_fits[i, "mean_dbh"]
  nci_mean <- sub_fits[i, "mean_nci"]
  max_dbh <- sub_fits[i, "max_dbh"]
  min_pet <- sub_fits[i, "min_pet"]
  max_pet <- sub_fits[i, "max_pet"]
  
  # Add curve
  if(effect == "size"){
    curve(size_eff(x), from = 0.1, to = max_dbh, add = T, col = cols[i])
  } else if(effect == "pet"){
    curve(pet_eff(x), from = min_pet, to = max_pet,
          add = T, col = cols[i])
  }
}
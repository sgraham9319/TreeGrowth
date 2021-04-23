
# Plotting size effect
size_effect <- function(tree_dat, fits, species, mod_structure){
  
  # Define different colors for different levels of model fit
  cols <- rep("green", times = nrow(fits))
  cols[which(fits$dAICc > 2)] <- "blue"  
  cols[which(fits$dAICc > 10)] <- "red"
  
  # Create dataset with one row per focal
  focals <- tree_dat %>%
    group_by(tree_id) %>%
    summarize(dbh = dbh[1], annual_growth = annual_growth[1],
              pet_dm = pet_dm[1])
  
  # Calculate mean NCI for each model run (if NCI included in model)
  if(grep("no", mod_structure) == 0){
    
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
    model_nci_means <- rep(NA, times = nrow(fits))
    for(i in 1:length(model_nci_means)){
      model_nci_means[i] <- mean(nci_test(tree_dat, fits[i, "alpha_opt"],
                                          fits[i, "beta_opt"]))
    }
  }
  
  # Plot growth vs. dbh
  plot(focals$annual_growth ~ focals$dbh,
       main = paste(species, "-", mod_structure, sep = " "),
       ylab = "Diameter growth (mm/y)", xlab = "DBH (dm)")
  
  # Define parameter fit plotting function
  if(grep("no", mod_structure) == 1){
    param_fit <- function(x){
      gmax * exp((-0.5) * ((log(x / X0) / Xb) ^ 2)) *
        exp((-0.5) * (((mean(focals$pet_dm) - pet_a) / pet_b) ^ 2))
    }
  } else if(grep("eq", mod_structure) == 1){
    param_fit <- function(x){
      gmax * exp((-0.5) * ((log(x / X0) / Xb) ^ 2)) *
        exp((-0.5) * (((mean(focals$pet_dm) - pet_a) / pet_b) ^ 2)) *
        exp(-C * nci_mean)
    }
  }
  
  # Add model fits
  for(i in 1:nrow(fits)){
    gmax <- fits$gmax_opt[i]
    X0 <- fits$X0_opt[i]
    Xb <- fits$Xb_opt[i]
    pet_a <- fits$pet_a_opt[i]
    pet_b <- fits$pet_b_opt[i]
    if(grep("no", mod_structure) == 0){
      C <- fits$C_opt[i]
      nci_mean <- model_nci_means[i]
    }
    curve(param_fit(x), from = 0.1, to = max(focals$dbh), add = T, col = cols[i])
  }
}

# Plotting PET effect
PET_effect <- function(tree_dat, fits, species, mod_structure){
  
  # Define different colors for different levels of model fit
  cols <- rep("green", times = nrow(fits))
  cols[which(fits$dAICc > 2)] <- "blue"  
  cols[which(fits$dAICc > 10)] <- "red"
  
  # Plot growth vs. PET
  plot(tree_dat$annual_growth ~ tree_dat$pet_dm,
       main = paste(species, "-", mod_structure, sep = " "),
       ylab = "Diameter growth (mm/y)", xlab = "PET (dm)")
  
  # Define parameter fit plotting function
  param_fit <- function(x){
    gmax * exp((-0.5) * ((log(mean(tree_dat$dbh) / X0) / Xb) ^ 2)) *
      exp((-0.5) * (((x - pet_a) / pet_b) ^ 2))
  }
  
  # Add model fits
  for(i in 1:nrow(fits)){
    gmax <- fits$gmax_opt[i]
    X0 <- fits$X0_opt[i]
    Xb <- fits$Xb_opt[i]
    pet_a <- fits$pet_a_opt[i]
    pet_b <- fits$pet_b_opt[i]
    curve(param_fit(x), from = 0.1, to = max(tree_dat$pet_dm), add = T, col = cols[i])
  }
}

# Plotting PET effect
PET_effect <- function(tree_dat, fits, species, mod_structure){
  
  # Define different colors for different levels of model fit
  cols <- rep("green", times = nrow(fits))
  cols[which(fits$dAICc > 2)] <- "blue"  
  cols[which(fits$dAICc > 10)] <- "red"
  
  # Create dataset with one row per focal
  focals <- tree_dat %>%
    group_by(tree_id) %>%
    summarize(dbh = dbh[1], annual_growth = annual_growth[1],
              pet_dm = pet_dm[1])
  
  # Calculate mean NCI for each model run (if NCI included in model)
  if(grep("no", mod_structure) == 0){
    
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
    model_nci_means <- rep(NA, times = nrow(fits))
    for(i in 1:length(model_nci_means)){
      model_nci_means[i] <- mean(nci_test(tree_dat, fits[i, "alpha_opt"],
                                          fits[i, "beta_opt"]))
    }
  }
  
  # Plot growth vs. PET
  plot(focals$annual_growth ~ focals$pet_dm,
       main = paste(species, "-", mod_structure, sep = " "),
       ylab = "Diameter growth (mm/y)", xlab = "PET (dm)")
  
  # Define parameter fit plotting function
  if(grep("no", mod_structure) == 1){
    param_fit <- function(x){
      gmax * exp((-0.5) * ((log(mean(focals$dbh) / X0) / Xb) ^ 2)) *
        exp((-0.5) * (((x - pet_a) / pet_b) ^ 2))
    }
  } else if(grep("eq", mod_structure) == 1){
    param_fit <- function(x){
      gmax * exp((-0.5) * ((log(mean(focals$dbh) / X0) / Xb) ^ 2)) *
        exp((-0.5) * (((x - pet_a) / pet_b) ^ 2)) *
        exp(-C * nci_mean)
    }
  }
  
  # Add model fits
  for(i in 1:nrow(fits)){
    gmax <- fits$gmax_opt[i]
    X0 <- fits$X0_opt[i]
    Xb <- fits$Xb_opt[i]
    pet_a <- fits$pet_a_opt[i]
    pet_b <- fits$pet_b_opt[i]
    if(grep("no", mod_structure) == 0){
      C <- fits$C_opt[i]
      nci_mean <- model_nci_means[i]
    }
    curve(param_fit(x), from = 0.1, to = max(focals$pet_dm), add = T, col = cols[i])
  }
}

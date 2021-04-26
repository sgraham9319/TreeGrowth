

fitted_effect <- function(tree_dat, fits, species, mod_structure, effect){
  
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
  if(length(grep("eq", mod_structure)) == 1){
    
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
    model_nci_means <- rep(NA, times = nrow(fits))
    for(i in 1:length(model_nci_means)){
      model_nci_means[i] <- mean(nci_test(tree_dat, fits[i, "alpha_opt"],
                                          fits[i, "beta_opt"]))
    }
  } else if(length(grep("int", mod_structure)) == 1){
    
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
    model_nci_means <- rep(NA, times = nrow(fits))
    for(i in 1:length(model_nci_means)){
      model_nci_means[i] <- mean(nci_test(tree_dat, fits[i, "alpha_opt"],
                                          fits[i, "beta_opt"],
                                          fits[i, "intra_opt"],
                                          fits[i, "inter_opt"]))
    }
  }
  
  # Make plot
  if(effect == "size"){
    plot(focals$annual_growth ~ focals$dbh,
         main = paste(species, "-", mod_structure, sep = " "),
         ylab = "Diameter growth (mm/y)", xlab = "DBH (dm)")
    
    # Define parameter fit plotting function
    if(length(grep("no", mod_structure)) == 1){
      param_fit <- function(x){
        gmax * exp((-0.5) * ((log(x / X0) / Xb) ^ 2)) *
          exp((-0.5) * (((mean(focals$pet_dm) - pet_a) / pet_b) ^ 2))
      }
    } else {
      param_fit <- function(x){
        gmax * exp((-0.5) * ((log(x / X0) / Xb) ^ 2)) *
          exp((-0.5) * (((mean(focals$pet_dm) - pet_a) / pet_b) ^ 2)) *
          exp(-C * nci_mean)
      }
    }
  } else if(effect == "pet"){
    
    plot(focals$annual_growth ~ focals$pet_dm,
         main = paste(species, "-", mod_structure, sep = " "),
         ylab = "Diameter growth (mm/y)", xlab = "PET (dm)")
    
    # Define parameter fit plotting function
    if(length(grep("no", mod_structure)) == 1){
      param_fit <- function(x){
        gmax * exp((-0.5) * ((log(mean(focals$dbh) / X0) / Xb) ^ 2)) *
          exp((-0.5) * (((x - pet_a) / pet_b) ^ 2))
      }
    } else {
      param_fit <- function(x){
        gmax * exp((-0.5) * ((log(mean(focals$dbh) / X0) / Xb) ^ 2)) *
          exp((-0.5) * (((x - pet_a) / pet_b) ^ 2)) *
          exp(-C * nci_mean)
      }
    }
  }
  
  # Add model fits
  for(i in 1:nrow(fits)){
    gmax <- fits$gmax_opt[i]
    X0 <- fits$X0_opt[i]
    Xb <- fits$Xb_opt[i]
    pet_a <- fits$pet_a_opt[i]
    pet_b <- fits$pet_b_opt[i]
    if(length(grep("no", mod_structure)) == 0){
      C <- fits$C_opt[i]
      nci_mean <- model_nci_means[i]
    }
    if(effect == "size"){
      curve(param_fit(x), from = 0.1, to = max(focals$dbh),
            add = T, col = cols[i])
    } else if (effect == "pet"){
      curve(param_fit(x), from = 0.1, to = max(focals$pet_dm),
            add = T, col = cols[i])
    }
  }
}
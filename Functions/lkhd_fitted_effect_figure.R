
lkhd_fitted_effect_figure <- function(effect){
  
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
      
      # Subset training to focal species and remove unneeded columns
      sing_sp <- train %>%
        arrange(tree_id) %>%
        filter(species == set_fits$focal_sps[sps]) %>%
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
        
        # Define NCI function
        nci <- function(neighbors, alpha, beta){
          raw <- (neighbors$dbh_comp ^ alpha) / (neighbors$prox ^ beta)
          return(sum(raw))
        }
        
        # Create function to calculate NCI for each individual
        nci_test <- function(nbhd_data, alpha, beta){
          ids <- unique(nbhd_data$tree_id)
          nci_vals <- rep(NA, times = length(ids))
          for(id in 1:length(ids)){
            neighbors <- nbhd_data[nbhd_data$tree_id == ids[id], ]
            nci_vals[id] <- nci(neighbors, alpha, beta)
          }
          nci_vals
        }
        
        # Calculate and store average NCI
        set_fits$mean_nci[sps] <- mean(nci_test(sing_sp,
                                                set_fits[sps, "alpha_opt"],
                                                set_fits[sps, "beta_opt"]))
        
      } else if(set_fits$str[sps] == "int_comp"){
        
        # Define NCI function
        nci <- function(neighbors, alpha, beta, intra, inter){
          focal_sps <- neighbors$species[1]
          raw <- (neighbors$dbh_comp ^ alpha) / (neighbors$prox ^ beta)
          cons <- which(neighbors$sps_comp == focal_sps)
          hets <- which(neighbors$sps_comp != focal_sps)
          nci_con <- raw[cons] * intra
          nci_het <- raw[hets] * inter
          return(sum(c(nci_con, nci_het)))
        }
        
        # Create function to calculate NCI for each individual
        nci_test <- function(nbhd_data, alpha, beta, intra, inter){
          ids <- unique(nbhd_data$tree_id)
          nci_vals <- rep(NA, times = length(ids))
          for(id in 1:length(ids)){
            neighbors <- nbhd_data[nbhd_data$tree_id == ids[id], ]
            nci_vals[id] <- nci(neighbors, alpha, beta, intra, inter)
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
        comm_comp <- comm_comp[, set_fits$focal_sps[sps]]
        
        # Convert rare competitor species to OTHR
        sing_sp <- sing_sp %>%
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
          
          # Create function to calculate NCI for each individual
          nci_test <- function(nbhd_data, alpha, beta, lmd1, lmd2, lmd3, lmd4){
            ids <- unique(nbhd_data$tree_id)
            nci_vals <- rep(NA, times = length(ids))
            for(id in 1:length(ids)){
              neighbors <- nbhd_data[nbhd_data$tree_id == ids[id], ]
              nci_vals[id] <- nci(neighbors, alpha, beta, lmd1, lmd2, lmd3,
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
          
          # Create function to calculate NCI for each individual
          nci_test <- function(nbhd_data, alpha, beta, lmd1, lmd2, lmd3,
                               lmd4, lmd5){
            ids <- unique(nbhd_data$tree_id)
            nci_vals <- rep(NA, times = length(ids))
            for(id in 1:length(ids)){
              neighbors <- nbhd_data[nbhd_data$tree_id == ids[id], ]
              nci_vals[id] <- nci(neighbors, alpha, beta, lmd1, lmd2, lmd3,
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
          
          # Create function to calculate NCI for each individual
          nci_test <- function(nbhd_data, alpha, beta, lmd1, lmd2, lmd3,
                               lmd4, lmd5, lmd6, lmd7, lmd8){
            ids <- unique(nbhd_data$tree_id)
            nci_vals <- rep(NA, times = length(ids))
            for(id in 1:length(ids)){
              neighbors <- nbhd_data[nbhd_data$tree_id == ids[id], ]
              nci_vals[id] <- nci(neighbors, alpha, beta, lmd1, lmd2, lmd3,
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
          
          # Create function to calculate NCI for each individual
          nci_test <- function(nbhd_data, alpha, beta, lmd1, lmd2, lmd3,
                               lmd4, lmd5, lmd6, lmd7, lmd8, lmd9){
            ids <- unique(nbhd_data$tree_id)
            nci_vals <- rep(NA, times = length(ids))
            for(id in 1:length(ids)){
              neighbors <- nbhd_data[nbhd_data$tree_id == ids[id], ]
              nci_vals[id] <- nci(neighbors, alpha, beta, lmd1, lmd2, lmd3,
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
  
  # Adjust plotting window parameters
  par(mfcol = c(2, 2))
  par(mar = c(0,0,0,0))
  par(oma = c(4,4,0.1,0.1))
  par(mgp = c(3, 0.3, 0))
  
  # Create first empty plot
  if(effect == "size"){
    plot(NULL, ylim = c(0, 3), xlim = c(0, 26),
         yaxt = "n", xaxt = "n", xlab = "", ylab = "")
    text(x = 12.5, y = 2.8, label = "Training 1", font = 2)
    axis(side = 2, labels = seq(0, 2.5, 0.5), at = seq(0, 2.5, 0.5), las = 1,
         tck = -0.02)
  } else if(effect == "pet"){
    plot(NULL, ylim = c(0, 3), xlim = c(min(fits$min_pet), 5.8),
         yaxt = "n", xaxt = "n", xlab = "", ylab = "")
    text(x = 4, y = 2.8, label = "Training 1", font = 2)
    axis(side = 2, labels = seq(0, 2.5, 0.5), at = seq(0, 2.5, 0.5), las = 1,
         tck = -0.02)
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
    plot(NULL, ylim = c(0, 3), xlim = c(0, 26),
         yaxt = "n", xaxt = "n", xlab = "", ylab = "")
    text(x = 12.5, y = 2.8, label = "Training 2", font = 2)
    axis(side = 2, labels = seq(0, 2.5, 0.5), at = seq(0, 2.5, 0.5), las = 1,
         tck = -0.02)
    axis(side = 1, labels = seq(0, 250, 50), at = seq(0, 25, 5), tck = -0.02)
  } else if(effect == "pet"){
    plot(NULL, ylim = c(0, 3), xlim = c(min(fits$min_pet), 5.8),
         yaxt = "n", xaxt = "n", xlab = "", ylab = "")
    text(x = 4, y = 2.8, label = "Training 2", font = 2)
    axis(side = 2, labels = seq(0, 2.5, 0.5), at = seq(0, 2.5, 0.5), las = 1,
         tck = -0.02)
    axis(side = 1, labels = seq(2.5, 5.5, 1), at = seq(2.5, 5.5, 1), tck = -0.02)
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
    plot(NULL, ylim = c(0, 3), xlim = c(0, 26),
         yaxt = "n", xaxt = "n", xlab = "", ylab = "")
    text(x = 12.5, y = 2.8, label = "Training 3", font = 2)
  } else if(effect == "pet"){
    plot(NULL, ylim = c(0, 3), xlim = c(min(fits$min_pet), 5.8),
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
    plot(NULL, ylim = c(0, 3), xlim = c(0, 26),
         yaxt = "n", xaxt = "n", xlab = "", ylab = "")
    text(x = 12.5, y = 2.8, label = "Training 4", font = 2)
    axis(side = 1, labels = seq(0, 250, 50), at = seq(0, 25, 5), tck = -0.02)
  } else if(effect == "pet"){
    plot(NULL, ylim = c(0, 3), xlim = c(min(fits$min_pet), 5.8),
         yaxt = "n", xaxt = "n", xlab = "", ylab = "")
    text(x = 4, y = 2.8, label = "Training 4", font = 2)
    axis(side = 1, labels = seq(2.5, 5.5, 1), at = seq(2.5, 5.5, 1), tck = -0.02)
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
  
  # Add y-axis label
  mtext("Predicted annual radial growth (cm)", side = 2, line = 2, outer = T)
  
  # Add x axis label
  if(effect == "size"){
    mtext("Focal tree DBH (cm)", side = 1, line = 2, outer = T)
  } else if(effect == "pet"){
    mtext("PET (dm)", side = 1, line = 2, outer = T)
  }
  
  # Add legend
  if(effect == "size"){
    legend(x = 20, y = 1.55,
           legend = c("ABAM", "CANO", "PSME", "THPL", "TSHE", "TSME"),
           col = cols, pch = 15, bty = "n", cex = 0.7)
  } else if(effect == "pet"){
    legend(x = 4.7, y = 1.55,
           legend = c("ABAM", "CANO", "PSME", "THPL", "TSHE", "TSME"),
           col = cols, pch = 15, bty = "n", cex = 0.7)
  }
}

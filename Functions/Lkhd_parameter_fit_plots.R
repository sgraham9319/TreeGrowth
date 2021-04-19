


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

PET_effect(focals, output, focal_sps, model_str)


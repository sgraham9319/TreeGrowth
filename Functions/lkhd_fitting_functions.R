
# Create growth prediction function
growth_pred <- function(nbhd_data, X0, Xb, gmax, pet_a, pet_b,
                        C = NULL, alpha = NULL, beta = NULL,
                        lmd1 = NULL, lmd2 = NULL, lmd3 = NULL,
                        lmd4 = NULL, lmd5 = NULL, lmd6 = NULL,
                        lmd7 = NULL, lmd8 = NULL, lmd9 = NULL){
  
  # Get list of focal tree ids
  ids <- unique(nbhd_data$tree_id)
  
  # Create vector to store growth predictions
  pred_grow <- rep(NA, times = length(ids))
  
  # Loop through focal trees
  for(id in 1:length(ids)){
    
    # Isolate neighborhood data for one focal tree
    neighbors <- nbhd_data %>%
      filter(tree_id == ids[id])
    
    # Predict growth
    pred_grow[id] <- gmax * 
      exp((-0.5) * (log(neighbors$dbh[1] / X0) / Xb) ^ 2) *
      exp((-0.5) * ((neighbors$pet_dm[1] - pet_a) / pet_b) ^ 2)
    
    # Modify growth prediction if required
    if(!is.null(C)){
      
      if(is.null(lmd1)){ # Equivalent competition modifier
        pred_grow[id] <- pred_grow[id] *
          exp(-C * nci_eq(neighbors, alpha, beta))
        
      } else if(is.null(lmd3)){ # Inter vs. inter competition modifier
        pred_grow[id] <- pred_grow[id] *
          exp(-C * nci_int(neighbors, alpha, beta, lmd1, lmd2))
        
      } else if(is.null(lmd5)){ # Species-specific competition (4 sps) modifier
        pred_grow[id] <- pred_grow[id] *
          exp(-C * nci_ss4(neighbors, alpha, beta, lmd1, lmd2, lmd3, lmd4))
        
      } else if(is.null(lmd6)){ # Species-specific competition (5 sps) modifier
        pred_grow[id] <- pred_grow[id] *
          exp(-C * nci_ss5(neighbors, alpha, beta, lmd1, lmd2, lmd3, lmd4,
                           lmd5))
        
      } else if(is.null(lmd9)){ # Species-specific competition (8 sps) modifier
        pred_grow[id] <- pred_grow[id] *
          exp(-C * nci_ss8(neighbors, alpha, beta, lmd1, lmd2, lmd3, lmd4,
                           lmd5, lmd6, lmd7, lmd8))
        
      } else if(!is.null(lmd9)){ # Species-specific competition (9 sps) modifier
        pred_grow[id] <- pred_grow[id] *
          exp(-C * nci_ss9(neighbors, alpha, beta, lmd1, lmd2, lmd3, lmd4,
                           lmd5, lmd6, lmd7, lmd8, lmd9))
        
      }
    }
  }
  
  # Join growth predictions with tree ids and output result
  results <- data.frame(ids, pred_grow, stringsAsFactors = F)
  return(results)
  
}

# Define negative log likelihood function to be minimized
neg_log_lkhd <- function(par){
  
  # Determine model structure using length of par
  if(length(par) == 6){
    mod_str <- "no_comp"
  } else if(length(par) == 9){
    mod_str <- "eq_comp"
  } else if(length(par) == 11){
    mod_str <- "int_comp"
  } else{
    mod_str <- "ss_comp"
    num_comps <- length(par) - 9
  }
  
  # Define parameters
  X0 <- par[1]
  Xb <- par[2]
  gmax <- par[3]
  pet_a <- par[4]
  pet_b <- par[5]
  if(mod_str != "no_comp"){
    C <- par[6]
    alpha <- par[7]
    beta <- par[8]
  }
  if(mod_str == "int_comp"){
    intra <- par[9]
    inter <- par[10]
  }
  if(mod_str == "ss_comp"){
    lmd1 <- par[9]
    lmd2 <- par[10]
    lmd3 <- par[11]
    lmd4 <- par[12]
    if(num_comps > 4){
      lmd5 <- par[13]
    }
    if(num_comps > 5){
      lmd6 <- par[14]
      lmd7 <- par[15]
      lmd8 <- par[16] 
    }
    if(num_comps == 9){
      lmd9 <- par[17]
    }
  }
  sigma <- par[length(par)]
  
  # Prevent parameter values from becoming nonsensical
  if(sigma < 0) {return(Inf)}
  if(X0 < 0 | X0 > 30) {return(Inf)}
  if(Xb < 0 | Xb > 3) {return(Inf)}
  if(gmax < 0) {return(Inf)}
  if(pet_a < 2 | pet_a > 6) {return(Inf)}
  if(pet_b < 0 | pet_b > 3) {return(Inf)}
  if(mod_str != "no_comp"){
    if(C < 0 | C > 10) {return(Inf)}
    if(alpha < 0 | alpha > 4) {return(Inf)}
    if(beta < 0 | beta > 4) {return(Inf)}
  }
  if(mod_str == "int_comp"){
    if(intra < 0 | intra > 1) {return(Inf)}
    if(inter < 0 | inter > 1) {return(Inf)}
  }
  if(mod_str == "ss_comp"){
    if(lmd1 < 0 | lmd1 > 1) {return(Inf)}
    if(lmd2 < 0 | lmd2 > 1) {return(Inf)}
    if(lmd3 < 0 | lmd3 > 1) {return(Inf)}
    if(lmd4 < 0 | lmd4 > 1) {return(Inf)}
    if(num_comps > 4){
      if(lmd5 < 0 | lmd5 > 1) {return(Inf)}
    }
    if(num_comps > 5){
      if(lmd6 < 0 | lmd6 > 1) {return(Inf)}
      if(lmd7 < 0 | lmd7 > 1) {return(Inf)}
      if(lmd8 < 0 | lmd8 > 1) {return(Inf)} 
    }
    if(num_comps == 9){
      if(lmd9 < 0 | lmd9 > 1) {return(Inf)}
    }
  }
  
  # Make growth predictions
  if(mod_str == "no_comp"){
    pred <- growth_pred(sing_sp_t, X0, Xb, gmax, pet_a, pet_b)
  } else if(mod_str == "eq_comp"){
    pred <- growth_pred(sing_sp_t, X0, Xb, gmax, pet_a, pet_b, C, alpha, beta)
  } else if(mod_str == "int_comp"){
    pred <- growth_pred(sing_sp_t, X0, Xb, gmax, pet_a, pet_b, C, alpha,
                        beta, intra, inter)
  } else if(mod_str == "ss_comp"){
    if(num_comps == 4){
      pred <- growth_pred(sing_sp_t, X0, Xb, gmax, pet_a, pet_b, C, alpha, beta,
                          lmd1, lmd2, lmd3, lmd4)
    } else if(num_comps == 5){
      pred <- growth_pred(sing_sp_t, X0, Xb, gmax, pet_a, pet_b, C, alpha, beta,
                          lmd1, lmd2, lmd3, lmd4, lmd5)
    } else if(num_comps == 8){
      pred <- growth_pred(sing_sp_t, X0, Xb, gmax, pet_a, pet_b, C, alpha, beta,
                          lmd1, lmd2, lmd3, lmd4, lmd5, lmd6, lmd7, lmd8)
    } else if(num_comps == 9){
      pred <- growth_pred(sing_sp_t, X0, Xb, gmax, pet_a, pet_b, C, alpha, beta,
                          lmd1, lmd2, lmd3, lmd4, lmd5, lmd6, lmd7, lmd8, lmd9)
    }
  }
  
  # Join predictions to observations by tree_id
  combined <- left_join(focals_t, pred, by = c("tree_id" = "ids"))
  
  # Calculate negative log likelihood
  NLL <- -sum(dnorm(combined$annual_growth, mean = combined$pred_grow,
                    sd = sigma, log = T))
  
  # Return value
  return(NLL)
  
}

# Create function for running optimization
nll_opt <- function(par_list){
  optim(par = par_list[1,], fn = neg_log_lkhd, method = "SANN")
}

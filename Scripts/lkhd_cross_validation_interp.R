library(dplyr)

# Load required functions
source("Functions/lkhd_fitting_functions.R")
source("Functions/nci.R")

# Define coefficient of determination function
coef_det <- function(x){
  1 - (sum((x$observations - x$predictions)^2) / 
         sum((x$observations - mean(x$observations))^2))
}

# Define focal species
focal_sps <- "PSME"

# Load model output
output <- read.csv(paste("Data/Output_data/cv1_", focal_sps, ".csv", sep = ""))

# Add model structure column
output$mod_str <- rep(c("no_comp", "eq_comp", "int_comp", "ss_comp"),
                      each = nrow(output) / 4)

# Create table to store fit results
fit_comp <- as.data.frame(matrix(NA, nrow = 2, ncol = 4))
names(fit_comp) <- c("mod_sel", "mod_str", "index", "test_r2")
fit_comp$mod_sel <- c("AICc", "cv")

# Store index and structure of best models according to AICc and cv
fit_comp$index <- c(which.min(output$AICc), which.max(output$cv_r2))
fit_comp$mod_str <- output$mod_str[fit_comp$index]

# Load test data
test <- read.csv("Data/Output_data/test1.csv")

# Load common competitors data and extract for focal species
comm_comp <- read.csv("Data/Output_data/common_comps.csv", stringsAsFactors = F)
comm_comp <- comm_comp[, focal_sps]

# Subset to focal species and remove unneeded columns
ss_test <- test %>%
  arrange(tree_id) %>%
  filter(species == focal_sps) %>%
  select(tree_id, stand_id, species, dbh, prox, sps_comp, dbh_comp,
         annual_growth, pet_mm)

# Change units of variables to give them similar ranges
ss_test <- ss_test %>%
  mutate(
    dbh = dbh / 10,                     # cm to dm
    dbh_comp = dbh_comp / 10,           # cm to dm
    annual_growth = annual_growth * 10, # cm to mm
    pet_dm = pet_mm / 100               # mm to dm
  ) %>%
  select(-pet_mm)

# Change rare competitors to OTHR
ss_test <- ss_test %>%
  mutate(sps_comp = if_else(sps_comp %in% comm_comp, sps_comp, "OTHR"))

# Create vector of competitor species
comps <- sort(unique(ss_test$sps_comp))

# Extract annual growth of each focal individual
focals <- ss_test %>%
  group_by(tree_id) %>%
  summarize(dbh = dbh[1], annual_growth = annual_growth[1])

# Calculate test R2 for best AICc and cv models
for(i in 1:2){
  
  # Make growth predictions using best model
  if(fit_comp$mod_str[i] == "no_comp"){
    growth_test <- growth_pred(ss_test,
                               output[fit_comp$index[i], "X0_opt"],
                               output[fit_comp$index[i], "Xb_opt"],
                               output[fit_comp$index[i], "gmax_opt"],
                               output[fit_comp$index[i], "pet_a_opt"],
                               output[fit_comp$index[i], "pet_b_opt"])
  } else if(fit_comp$mod_str[i] == "eq_comp"){
    growth_test <- growth_pred(ss_test,
                               output[fit_comp$index[i], "X0_opt"],
                               output[fit_comp$index[i], "Xb_opt"],
                               output[fit_comp$index[i], "gmax_opt"],
                               output[fit_comp$index[i], "pet_a_opt"],
                               output[fit_comp$index[i], "pet_b_opt"],
                               output[fit_comp$index[i], "C_opt"],
                               output[fit_comp$index[i], "alpha_opt"],
                               output[fit_comp$index[i], "beta_opt"])
  } else if(fit_comp$mod_str[i] == "int_comp"){
    growth_test <- growth_pred(ss_test,
                               output[fit_comp$index[i], "X0_opt"],
                               output[fit_comp$index[i], "Xb_opt"],
                               output[fit_comp$index[i], "gmax_opt"],
                               output[fit_comp$index[i], "pet_a_opt"],
                               output[fit_comp$index[i], "pet_b_opt"],
                               output[fit_comp$index[i], "C_opt"],
                               output[fit_comp$index[i], "alpha_opt"],
                               output[fit_comp$index[i], "beta_opt"],
                               output[fit_comp$index[i], "intra_opt"],
                               output[fit_comp$index[i], "inter_opt"])
  } else if(fit_comp$mod_str[i] == "ss_comp"){
    if(length(comps) == 4){
      growth_test <- growth_pred(ss_test,
                                 output[fit_comp$index[i], "X0_opt"],
                                 output[fit_comp$index[i], "Xb_opt"],
                                 output[fit_comp$index[i], "gmax_opt"],
                                 output[fit_comp$index[i], "pet_a_opt"],
                                 output[fit_comp$index[i], "pet_b_opt"],
                                 output[fit_comp$index[i], "C_opt"],
                                 output[fit_comp$index[i], "alpha_opt"],
                                 output[fit_comp$index[i], "beta_opt"],
                                 output[fit_comp$index[i], "lmd1_opt"],
                                 output[fit_comp$index[i], "lmd2_opt"],
                                 output[fit_comp$index[i], "lmd3_opt"],
                                 output[fit_comp$index[i], "lmd4_opt"])
    } else if(length(comps) == 5){
      growth_test <- growth_pred(ss_test,
                                 output[fit_comp$index[i], "X0_opt"],
                                 output[fit_comp$index[i], "Xb_opt"],
                                 output[fit_comp$index[i], "gmax_opt"],
                                 output[fit_comp$index[i], "pet_a_opt"],
                                 output[fit_comp$index[i], "pet_b_opt"],
                                 output[fit_comp$index[i], "C_opt"],
                                 output[fit_comp$index[i], "alpha_opt"],
                                 output[fit_comp$index[i], "beta_opt"],
                                 output[fit_comp$index[i], "lmd1_opt"],
                                 output[fit_comp$index[i], "lmd2_opt"],
                                 output[fit_comp$index[i], "lmd3_opt"],
                                 output[fit_comp$index[i], "lmd4_opt"],
                                 output[fit_comp$index[i], "lmd5_opt"])
    } else if(length(comps) == 8){
      growth_test <- growth_pred(ss_test,
                                 output[fit_comp$index[i], "X0_opt"],
                                 output[fit_comp$index[i], "Xb_opt"],
                                 output[fit_comp$index[i], "gmax_opt"],
                                 output[fit_comp$index[i], "pet_a_opt"],
                                 output[fit_comp$index[i], "pet_b_opt"],
                                 output[fit_comp$index[i], "C_opt"],
                                 output[fit_comp$index[i], "alpha_opt"],
                                 output[fit_comp$index[i], "beta_opt"],
                                 output[fit_comp$index[i], "lmd1_opt"],
                                 output[fit_comp$index[i], "lmd2_opt"],
                                 output[fit_comp$index[i], "lmd3_opt"],
                                 output[fit_comp$index[i], "lmd4_opt"],
                                 output[fit_comp$index[i], "lmd5_opt"],
                                 output[fit_comp$index[i], "lmd6_opt"],
                                 output[fit_comp$index[i], "lmd7_opt"],
                                 output[fit_comp$index[i], "lmd8_opt"])
    } else if(length(comps) == 9){
      growth_test <- growth_pred(ss_test,
                                 output[fit_comp$index[i], "X0_opt"],
                                 output[fit_comp$index[i], "Xb_opt"],
                                 output[fit_comp$index[i], "gmax_opt"],
                                 output[fit_comp$index[i], "pet_a_opt"],
                                 output[fit_comp$index[i], "pet_b_opt"],
                                 output[fit_comp$index[i], "C_opt"],
                                 output[fit_comp$index[i], "alpha_opt"],
                                 output[fit_comp$index[i], "beta_opt"],
                                 output[fit_comp$index[i], "lmd1_opt"],
                                 output[fit_comp$index[i], "lmd2_opt"],
                                 output[fit_comp$index[i], "lmd3_opt"],
                                 output[fit_comp$index[i], "lmd4_opt"],
                                 output[fit_comp$index[i], "lmd5_opt"],
                                 output[fit_comp$index[i], "lmd6_opt"],
                                 output[fit_comp$index[i], "lmd7_opt"],
                                 output[fit_comp$index[i], "lmd8_opt"],
                                 output[fit_comp$index[i], "lmd9_opt"])
    }
  }
  
  # Combine predicted and observed growth
  obs_pred <- focals %>%
    left_join(growth_test, by = c("tree_id" = "ids")) %>%
    rename(observations = annual_growth,
           predictions = pred_grow)
  
  # Calculate and store coefficient of determination
  fit_comp$test_r2[i] <- coef_det(obs_pred)
}


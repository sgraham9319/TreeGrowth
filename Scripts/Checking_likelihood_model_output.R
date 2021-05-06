library(dplyr)

# Source model fit plotting functions
source("Functions/Lkhd_parameter_fit_plots.R")

# Specify model structure, training set and focal sps
model_str <- "ss_comp_rand"
set <- 2
focal_sps <- "ABAM"

# Load model output
output <- read.csv(paste("Data/Output_data/", model_str, set, "_", focal_sps,
                         ".csv", sep = ""))

# Load training data
#training <- read.csv(paste("Data/Output_data/training", set, ".csv", sep = ""),
#                   stringsAsFactors = F)
training <- read.csv(paste("Data/Output_data/rand_training", set, ".csv",
                           sep = ""), stringsAsFactors = F)

# Order training data by tree_id to avoid problems later
training <- training %>% arrange(tree_id)

# Subset to focal species and remove unneeded columns
sing_sp <- training %>%
  arrange(tree_id) %>%
  filter(species == focal_sps) %>%
  select(tree_id, stand_id, species, dbh, prox, sps_comp, dbh_comp,
         annual_growth, pet_mm)

# Change units of variables to give them similar ranges
sing_sp <- sing_sp %>%
  mutate(
    dbh = dbh / 10,                     # cm to dm
    dbh_comp = dbh_comp / 10,           # cm to dm
    annual_growth = annual_growth * 10, # cm to mm
    pet_dm = pet_mm / 100               # mm to dm
  ) %>%
  select(-pet_mm)

# If a single species comp model, change rare competitors to OTHR
if(length(grep("ss", model_str)) == 1){
  comm_comp <- read.csv("Data/Output_data/common_comps.csv",
                        stringsAsFactors = F)
  comm_comp <- comm_comp[, focal_sps]
  sing_sp <- sing_sp %>%
    mutate(sps_comp = if_else(sps_comp %in% comm_comp, sps_comp, "OTHR"))
}

# Extract annual growth of each focal individual
focals <- sing_sp %>%
  group_by(tree_id) %>%
  summarize(dbh = dbh[1], annual_growth = annual_growth[1],
            pet_dm = pet_dm[1])

# Plot starting vs. optimized values
plot(X0_opt ~ X0, data = output)
plot(Xb_opt ~ Xb, data = output)
plot(gmax_opt ~ gmax, data = output)
plot(pet_a_opt ~ pet_a, data = output)
plot(pet_b_opt ~ pet_b, data = output)
plot(C_opt ~ C, data = output)
plot(alpha_opt ~ alpha, data = output)
plot(beta_opt ~ beta, data = output)
plot(lmd1_opt ~ lmd1, data = output)
plot(lmd2_opt ~ lmd2, data = output)
plot(lmd3_opt ~ lmd3, data = output)
plot(lmd4_opt ~ lmd4, data = output)
plot(lmd5_opt ~ lmd5, data = output)
plot(lmd6_opt ~ lmd6, data = output)
plot(lmd7_opt ~ lmd7, data = output)
plot(lmd8_opt ~ lmd8, data = output)
plot(lmd9_opt ~ lmd9, data = output)

# Calculate AICc for each model
AICc_calc <- function(k, NLL, n){
  (2 * (k + NLL)) + ((2 * ((k^2) + k)) / (n - k - 1))
}
output$AICc <- AICc_calc(length(grep("_opt", names(output))), output$NLL, nrow(focals))

# Calculate dAICc for each model
output$dAICc <- output$AICc - min(output$AICc)

# Order output by dAICc
output <- output %>% arrange(dAICc)

# Plot size effect
fitted_effect(sing_sp, output, focal_sps, model_str, effect = "size")

# Plot PET effect
fitted_effect(sing_sp, output, focal_sps, model_str, effect = "pet")

# Link lmd parameters to competitor species
sort(unique(sing_sp$sps_comp))

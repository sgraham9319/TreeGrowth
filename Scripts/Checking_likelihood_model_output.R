library(dplyr)

# Source model fit plotting functions
source("Functions/Lkhd_parameter_fit_plots.R")

# Specify model structure, training set and focal sps
model_str <- "no_comp"
set <- 2
focal_sps <- "ABAM"

# Load model output
output <- read.csv(paste("Data/Output_data/", model_str, set, "_", focal_sps,
                         ".csv", sep = ""))

# Load training data
training <- read.csv(paste("Data/Output_data/training", set, ".csv", sep = ""),
                   stringsAsFactors = F)

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

# Calculate AICc for each model
AICc_calc <- function(k, NLL, n){
  (2 * (k + NLL)) + ((2 * ((k^2) + k)) / (n - k - 1))
}
output$AICc <- AICc_calc(length(grep("_opt", names(output))), output$NLL, nrow(focals))

# Calculate dAICc for each model
output$dAICc <- output$AICc - min(output$AICc)

# Order output by dAICc
output <- output %>% arrange(dAICc)

# Plot PET effect
PET_effect(focals, output, focal_sps, model_str)

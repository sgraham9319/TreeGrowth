library(dplyr)
library(ggplot2)

# Load plotting function
source("Functions/size_effect.R")

# Need an output table like this (column order matters!)
plot_data <- data.frame(
  label = paste(rep("Set", times = 4), 1:4),
  model_str = c("no_comp", "eq_comp", "eq_comp", "no_comp"),
  gmax = c(4, 3.5, 4.2, 4.5),
  X0 = c(8, 10, 20, 11),
  Xb = c(0.8, 1, 2, 0.9),
  mean_pet = c(3.9, 3.8, 3.95, 4),
  pet_a = c(4.5, 4.7, 4.3, 4.4),
  pet_b = c(2.2, 2.3, 2.5, 2.3),
  C = c(NA, 0.1, 0.05, NA),
  mean_nci = c(NA, 7, 4, NA),
  max_dbh = c(20, 25, 19, 30))

# Need a vector of colors
col_vec <- c("red", "blue", "gold", "green")





# Specify focal species and model structure
focal_sps <- "TSME"
model_str <- "no_comp"

set <- 1

# Load model fits
output <- read.csv(paste("Data/Output_data/", model_str, set, "_", focal_sps,
                         ".csv", sep = ""))

# Load training data
training <- read.csv(paste("Data/Output_data/training", set, ".csv", sep = ""),
                     stringsAsFactors = F)

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

# Calculate AICc for each model
AICc_calc <- function(k, NLL, n){
  (2 * (k + NLL)) + ((2 * ((k^2) + k)) / (n - k - 1))
}
output$AICc <- AICc_calc(length(grep("_opt", names(output))), output$NLL, nrow(focals))

# Calculate dAICc for each model
output$dAICc <- output$AICc - min(output$AICc)

# Order output by dAICc
output <- output %>% arrange(dAICc)


# Create example output table
plot_data <- data.frame(
  label = paste(rep("Set", times = 4), 1:4),
  model_str = c("no_comp", "eq_comp", "eq_comp", "no_comp"),
  gmax = c(4, 3.5, 4.2, 4.5),
  X0 = c(8, 10, 20, 11),
  Xb = c(0.8, 1, 2, 0.9),
  mean_pet = c(3.9, 3.8, 3.95, 4),
  pet_a = c(4.5, 4.7, 4.3, 4.4),
  pet_b = c(2.2, 2.3, 2.5, 2.3),
  C = c(NA, 0.1, 0.05, NA),
  mean_nci = c(NA, 7, 4, NA),
  max_dbh = c(20, 25, 19, 30))
col_vec <- c("red", "blue", "gold", "green")

size_effect(plot_data, focal_sps, col_vec)




library(tidyverse)
library(ForestPlot)

# Select focal species
focal_sps <- "PSME"

# Load training data
train <- read.csv("Data/Output_data/training1.csv", stringsAsFactors = F)

# Remove unneeded columns
train <- train %>%
  select(tree_id, species, sps_comp, dbh_comp, prox, all_density, ABAM_density,
         ABPR_density, CANO_density, PSME_density, TABR_density, THPL_density,
         TSHE_density, ABLA_density, TSME_density, ALSI_density, ALVI_density,
         PIMO_density, ALRU_density, PICO_density, PIEN_density, POBA_density,
         ABGR_density, PISI_density, pet_mm, size_corr_growth)

# Create additional explanatory variables
train <- train %>%
  mutate(
    intra = if_else(sps_comp == focal_sps, 1, 0),
    inter_dens = all_density - get(paste(focal_sps, "density", sep = "_")))

# Load common competitors data and extract for focal species
comm_comp <- read.csv("Data/Output_data/common_comps.csv", stringsAsFactors = F)
comm_comp <- comm_comp[, focal_sps]

# Convert rare competitors to OTHR
train <- train %>%
  mutate(sps_comp = if_else(sps_comp %in% comm_comp, sps_comp, "OTHR"))

# Run regularized regression model
mod_res <- growth_model(train, "size_corr_growth", focal_sps,
                        iterations = 100)


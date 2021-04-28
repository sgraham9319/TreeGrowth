library(dplyr)
library(ggplot2)

# Load plotting functions
source("Functions/size_effect_comp.R")
source("Functions/training_comparison.R")

# Calculate comparison
comparison <- training_comparison(focal_sps = "TSHE",
                    model_strs = c("no_comp", "eq_comp"),
                    sets = 1:4)

# Plot comparison
comparison$comp_plot

# View all results
View(comparison$best_models)

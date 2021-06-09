library(dplyr)
library(ggplot2)

# Load plotting functions
source("Functions/size_effect_comp.R")
source("Functions/training_comparison.R")
source("Functions/lkhd_model_selection.R")
source("Functions/quantify_lkhd_fit.R")
source("Functions/r2_figure.R")

# Calculate comparison
comparison <- training_comparison(focal_sps = "TSME",
                    model_strs = c("no_comp", "eq_comp", "int_comp", "ss_comp"),
                    sets = 1:4, train_type = "random")

# Plot comparison
comparison$comp_plot

# View all results
View(comparison$best_models)

# Create likelihood AIC model selection table
lkhd_table_rand <- lkhd_model_select()

# Create likelihood cv model selection table
lkhd_table_cv <- lkhd_model_select(method = "cv", num_train_sets = 2)

# Save model selection tables
write.csv(lkhd_table_rand, "Figures/lkhd_model_selection_rand_train.csv",
          row.names = F)
write.csv(lkhd_table_cv, "Figures/lkhd_model_selection_cv.csv",
          row.names = F)

# Calculate likelihood r2 values
lkhd_r2_rand <- quantify_lkhd_fit()
lkhd_r2_cv <- quantify_lkhd_fit(method = "cv", num_train_sets = 2)

# Load regularized regression r2
rr_r2_rand <- read.csv("Data/Figure_data/RR_r2_rand.csv")

# Combine into single r2 table
r2_table <- lkhd_r2_rand %>%
  left_join(lkhd_r2_cv %>% select(-sample_size),
            by = c("focal_sps", "training_set"),
            suffix = c("_AIC", "_cv")) %>%
  left_join(rr_r2_rand, by = c("focal_sps" = "species",
                               "training_set" = "training")) %>%
  rename(train_rr = train_r2, test_rr = test_r2) %>%
  select(focal_sps, training_set, sample_size, train_lkhd_AIC, train_lkhd_cv,
         train_rr, test_lkhd_AIC, test_lkhd_cv, test_rr) %>%
  mutate(train_lkhd_AIC = round(train_lkhd_AIC, 2),
         train_lkhd_cv = round(train_lkhd_cv, 2),
         train_rr = round(train_rr, 2),
         test_lkhd_AIC = round(test_lkhd_AIC, 2),
         test_lkhd_cv = round(test_lkhd_cv, 2),
         test_rr = round(test_rr, 2))

# Save r2 table
write.csv(r2_table, "Figures/all_train_test_fit.csv", row.names = F)

# Make R2 figures
r2_all <- read.csv("Figures/all_train_test_fit.csv", stringsAsFactors = F)
r2_figure(r2_rand, "train_test_fit_rand")
r2_figure(r2_rand, "predict_fits", cv = T)

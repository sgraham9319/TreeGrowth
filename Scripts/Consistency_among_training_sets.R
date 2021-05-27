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

# Create likelihood model selection table
lkhd_table <- lkhd_model_select(training_type = "regular")
lkhd_table_rand <- lkhd_model_select(training_type = "random")

# Save model selection tables
#write.csv(lkhd_table, "Figures/lkhd_model_selection.csv", row.names = F)
#write.csv(lkhd_table_rand, "Figures/lkhd_model_selection_rand_train.csv",
#          row.names = F)

# Calculate likelihood r2 values
lkhd_r2 <- quantify_lkhd_fit(train_type = "regular")
lkhd_r2_rand <- quantify_lkhd_fit(train_type = "random")

# Load regularized regression r2
rr_r2 <- read.csv("Data/Figure_data/RR_r2.csv")
rr_r2_rand <- read.csv("Data/Figure_data/RR_r2_rand.csv")

# Combine into single r2 table
r2_table <- lkhd_r2 %>%
  left_join(rr_r2, by = c("focal_sps" = "species",
                          "training_set" = "training")) %>%
  rename(train_rr = train_r2, test_rr = test_r2) %>%
  select(focal_sps, training_set, sample_size, train_lkhd, train_rr, test_lkhd,
         test_rr) %>%
  mutate(train_lkhd = round(train_lkhd, 2),
         train_rr = round(train_rr, 2),
         test_lkhd = round(test_lkhd, 2),
         test_rr = round(test_rr, 2))
r2_rand_table <- lkhd_r2_rand %>%
  left_join(rr_r2_rand, by = c("focal_sps" = "species",
                               "training_set" = "training")) %>%
  rename(train_rr = train_r2, test_rr = test_r2) %>%
  select(focal_sps, training_set, sample_size, train_lkhd, train_rr, test_lkhd,
         test_rr) %>%
  mutate(train_lkhd = round(train_lkhd, 2),
         train_rr = round(train_rr, 2),
         test_lkhd = round(test_lkhd, 2),
         test_rr = round(test_rr, 2))

# Save r2 tables
write.csv(r2_table, "Figures/train_test_fit.csv", row.names = F)
write.csv(r2_rand_table, "Figures/train_test_fit_rand.csv", row.names = F)

# Make R2 figures
r2 <- read.csv("Figures/train_test_fit.csv", stringsAsFactors = F)
r2_rand <- read.csv("Figures/train_test_fit_rand.csv", stringsAsFactors = F)
r2_figure(r2, "train_test_fit")
r2_figure(r2_rand, "train_test_fit_rand")
r2_figure(r2_rand, "predict_fits", cv = T)

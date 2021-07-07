library(dplyr)
library(ggplot2)

# Load required plotting functions
source("Functions/lkhd_model_select.R")
source("Functions/quantify_lkhd_fit.R")
source("Functions/r2_figure.R")
source("Functions/lkhd_fitted_params.R")
source("Functions/nci.R")
source("Functions/lkhd_fitted_effect_figure.R")

#====================================================
# Tables S2-7: AIC likelihood fitted parameter values
#====================================================

# Create and save table of fitted parameter values for AIC likelihood
lkhd_params <- lkhd_fitted_params()
write.csv(lkhd_params, "Figures/lkhd_AIC_param_table.csv", row.names = F)

#============================================
# Figures S1-2: AIC likelihood fitted effects
#============================================

# Create size effect comparison figure
lkhd_fitted_effect_figure(effect = "size")

# Create pet effect comparison figure
lkhd_fitted_effect_figure(effect = "pet")

#=========================================
# Table S1: AIC likelihood model selection
#=========================================

# Create and save likelihood AIC model selection table
lkhd_table_rand <- lkhd_model_select()
write.csv(lkhd_table_rand, "Figures/lkhd_model_selection_rand_train.csv",
          row.names = F)

#======================================================
# Table S18: Cross-validated likelihood model selection
#======================================================

# Create and save likelihood cv model selection table
lkhd_table_cv <- lkhd_model_select(method = "cv", num_train_sets = 4)
write.csv(lkhd_table_cv, "Figures/lkhd_model_selection_cv.csv",
          row.names = F)

#=====================================================
# Table S17: Training and test fit for all model types
#=====================================================

# Quantify AIC likelihood fit
lkhd_r2_rand <- quantify_lkhd_fit()

# Quantify cross-validation likelihood fit
lkhd_r2_cv <- quantify_lkhd_fit(method = "cv", num_train_sets = 4)

# Load regularized regression fit
rr_r2_rand <- read.csv("Data/Figure_data/RR_r2_rand.csv")

# Combine into single fit table
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

# Save fit table
write.csv(r2_table, "Figures/all_train_test_fit.csv", row.names = F)

#===================================================
# Figure 3: training and test fit of all model types
#===================================================

# Make training and test fit figures
r2_all <- read.csv("Figures/all_train_test_fit.csv", stringsAsFactors = F)
r2_figure(r2_all, "train_test_fit_all")

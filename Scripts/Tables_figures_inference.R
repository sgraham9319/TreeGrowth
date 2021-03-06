library(tidyverse)
library(paletteer)

source("Functions/intra_vs_inter_lkhd.R")
source("Functions/sps_int_plot.R")

#=====================================================
# 1. Does neighborhood crowding influence tree growth?
#=====================================================

# Load likelihood model comparison table
lkhd_table_rd <- read.csv("Figures/lkhd_model_selection_rand_train.csv",
                          stringsAsFactors = F)

# Create table of best model structure for each species/training set
best_strs_rd <- lkhd_table_rd %>%
  mutate(str = names(lkhd_table_rd)[apply(lkhd_table_rd[, 3:6], 1, which.min) +
                                      2]) %>%
  select(focal_sps, training_set, str)

# Create likelihood part of results table
lkhd_nbhd_rd <- best_strs_rd %>%
  mutate(nbhd = if_else(str != "no_comp", "Yes", "No")) %>%
  select(-str) %>%
  pivot_wider(names_from = training_set, values_from = nbhd)

# Load regularized regression results
rr_nbhd_rd <- read.csv("Data/Figure_data/nbhd_inf_rand.csv")

# Combine likelihood and regularized regression data
nbhd_rd <- lkhd_nbhd_rd %>%
  left_join(rr_nbhd_rd, by = c("focal_sps" = "species"))
names(nbhd_rd)[2:9] <- c(paste("L", 1:4, sep = ""), paste("RR", 1:4, sep = ""))

# Save tables
write.csv(nbhd_rd, "Figures/nbhd_influence_rand.csv", row.names = F)

#=============================================================
# 2. Does species identity of neighbors influence tree growth?
#=============================================================

# Create likelihood part of results table
lkhd_compid_rd <- best_strs_rd %>%
  mutate(compid = if_else(str %in% c("int_comp", "ss_comp"), "Yes", "No")) %>%
  select(-str) %>%
  pivot_wider(names_from = training_set, values_from = compid)

# Load regularized regression results
rr_compid_rd <- read.csv("Data/Figure_data/comp_id_inf_rand.csv")

# Combine likelihood and regularized regression data
compid_rd <- lkhd_compid_rd %>%
  left_join(rr_compid_rd, by = c("focal_sps" = "species"))
names(compid_rd)[2:9] <- c(paste("L", 1:4, sep = ""),
                           paste("RR", 1:4, sep = ""))

# Save tables
write.csv(compid_rd, "Figures/comp_id_influence_rand.csv", row.names = F)

#=======================================================================
# 3. Is growth more impacted by conspecific or heterospecific neighbors?
#=======================================================================

# Create likelihood part of table 
lkhd_intra_rd <- intra_vs_inter_lkhd("random")

# Load regularized regression results
rr_het_rd <- read.csv("Data/Figure_data/inter_str_rand.csv")
rr_con_rd <- read.csv("Data/Figure_data/intra_str_rand.csv")

# Calculate regularized regression part of table
rr_intra_rd <- cbind(rr_het_rd[, 1], rr_het_rd[, 2:5] - rr_con_rd[, 2:5])
names(rr_intra_rd) <- c("species", paste("RR", 1:4, sep = ""))

# Combine likelihood and regularized regression data
intra_rd <- lkhd_intra_rd %>%
  left_join(rr_intra_rd, by = "species")

# Save tables
write.csv(intra_rd, "Figures/intra_vs_inter_rand.csv", row.names = F)

#====================================================
# 4. Which neighbor species have the largest effects?
#====================================================

# Create vector of species
species <- c("ABAM", "CANO", "PSME", "THPL", "TSHE", "TSME")

# Create plots
for(i in species){
  sps_int_plot(i, train_type = "random")
}

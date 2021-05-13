library(tidyverse)
library(paletteer)

source("Functions/intra_vs_inter_lkhd.R")
source("Functions/sps_int_plot.R")

#============================================================
# 1. Does the competitive neighborhood influence tree growth?
#============================================================

# Load likelihood model comparison tables
lkhd_table <- read.csv("Figures/lkhd_model_selection.csv", stringsAsFactors = F)
lkhd_table_rd <- read.csv("Figures/lkhd_model_selection_rand_train.csv",
                          stringsAsFactors = F)

# Create table of best model structure for each species/training set
best_strs <- lkhd_table %>%
  mutate(str = names(lkhd_table)[apply(lkhd_table[, 3:6], 1, which.min) +
                                   2]) %>%
  select(focal_sps, training_set, str)
best_strs_rd <- lkhd_table_rd %>%
  mutate(str = names(lkhd_table_rd)[apply(lkhd_table_rd[, 3:6], 1, which.min) +
                                      2]) %>%
  select(focal_sps, training_set, str)

# Create likelihood parts of results tables
lkhd_nbhd <- best_strs %>%
  mutate(nbhd = if_else(str != "no_comp", "Yes", "No")) %>%
  select(-str) %>%
  pivot_wider(names_from = training_set, values_from = nbhd)
lkhd_nbhd_rd <- best_strs_rd %>%
  mutate(nbhd = if_else(str != "no_comp", "Yes", "No")) %>%
  select(-str) %>%
  pivot_wider(names_from = training_set, values_from = nbhd)

# Load regularized regression results
rr_nbhd <- read.csv("Data/Figure_data/nbhd_inf.csv")
rr_nbhd_rd <- read.csv("Data/Figure_data/nbhd_inf_rand.csv")

# Combine likelihood and regularized regression data
nbhd <- lkhd_nbhd %>%
  left_join(rr_nbhd, by = c("focal_sps" = "species"))
names(nbhd)[2:9] <- c(paste("L", 1:4, sep = ""), paste("RR", 1:4, sep = "")) 
nbhd_rd <- lkhd_nbhd_rd %>%
  left_join(rr_nbhd_rd, by = c("focal_sps" = "species"))
names(nbhd_rd)[2:9] <- c(paste("L", 1:4, sep = ""), paste("RR", 1:4, sep = ""))

# Save tables
write.csv(nbhd, "Figures/nbhd_influence.csv", row.names = F)
write.csv(nbhd_rd, "Figures/nbhd_influence_rand.csv", row.names = F)

#=============================================================
# 2. Does species identity of neighbors influence tree growth?
#=============================================================

# Create likelihood parts of results tables
lkhd_compid <- best_strs %>%
  mutate(compid = if_else(str %in% c("int_comp", "ss_comp"), "Yes", "No")) %>%
  select(-str) %>%
  pivot_wider(names_from = training_set, values_from = compid)
lkhd_compid_rd <- best_strs_rd %>%
  mutate(compid = if_else(str %in% c("int_comp", "ss_comp"), "Yes", "No")) %>%
  select(-str) %>%
  pivot_wider(names_from = training_set, values_from = compid)

# Load regularized regression results
rr_compid <- read.csv("Data/Figure_data/comp_id_inf.csv")
rr_compid_rd <- read.csv("Data/Figure_data/comp_id_inf_rand.csv")

# Combine likelihood and regularized regression data
compid <- lkhd_compid %>%
  left_join(rr_compid, by = c("focal_sps" = "species"))
names(compid)[2:9] <- c(paste("L", 1:4, sep = ""), paste("RR", 1:4, sep = "")) 
compid_rd <- lkhd_compid_rd %>%
  left_join(rr_compid_rd, by = c("focal_sps" = "species"))
names(compid_rd)[2:9] <- c(paste("L", 1:4, sep = ""),
                           paste("RR", 1:4, sep = ""))

# Save tables
write.csv(compid, "Figures/comp_id_influence.csv", row.names = F)
write.csv(compid_rd, "Figures/comp_id_influence_rand.csv", row.names = F)

#===========================================================
# 3. Is conspecific or heterospecific competition strongest?
#===========================================================

# Create likelihood part of table 
lkhd_intra <- intra_vs_inter_lkhd("regular")
lkhd_intra_rd <- intra_vs_inter_lkhd("random")

# Load regularized regression results
rr_het <- read.csv("Data/Figure_data/inter_str.csv")
rr_con <- read.csv("Data/Figure_data/intra_str.csv")
rr_het_rd <- read.csv("Data/Figure_data/inter_str_rand.csv")
rr_con_rd <- read.csv("Data/Figure_data/intra_str_rand.csv")

# Calculate regularized regression part of table
rr_intra <- cbind(rr_het[, 1], rr_het[, 2:5] - rr_con[, 2:5])
names(rr_intra) <- c("species", paste("RR", 1:4, sep = ""))
rr_intra_rd <- cbind(rr_het_rd[, 1], rr_het_rd[, 2:5] - rr_con_rd[, 2:5])
names(rr_intra_rd) <- c("species", paste("RR", 1:4, sep = ""))

# Combine likelihood and regularized regression data
intra <- lkhd_intra %>%
  left_join(rr_intra, by = "species")
intra_rd <- lkhd_intra_rd %>%
  left_join(rr_intra_rd, by = "species")

# Save tables
write.csv(intra, "Figures/intra_vs_inter.csv", row.names = F)
write.csv(intra_rd, "Figures/intra_vs_inter_rand.csv", row.names = F)

#================================================
# 4. Which species are the strongest competitors?
#================================================

sps_int_plot("ABAM", train_type = "regular")
sps_int_plot("TSME", train_type = "random")
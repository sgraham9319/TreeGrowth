library(tidyverse)
library(ForestPlot)

# NOTE: this script takes a while to run (about 15 minutes on my laptop), but
# for checking purposes it can be run much faster if you reduce the number
# of focal species and training sets provided on lines 16 and 17 (the number
# of times the for loop runs is the product of the lengths of the focal_sps
# and train_sets vectors)

# Define function to calculate species interaction coefficients
rescale <- function(x){
  (sum(x < 0) - sum(x > 0) + 100) / 200
}

# Create vectors of focal species and training sets
focal_sps <- c("ABAM", "CANO", "PSME", "THPL", "TSHE", "TSME")
train_sets <- 1:4

# Create matrix to hold R squared values
rsq_vals <- matrix(NA, ncol = 2,
                   nrow = length(focal_sps) * length(train_sets))

# Create matrices to hold ecological interpretation results
nbhd_inf <- matrix(NA, ncol = length(train_sets), nrow = length(focal_sps))
comp_id_inf <- nbhd_inf
inter_str <- nbhd_inf
intra_str <- nbhd_inf

# Create list to store model coefficients
coef_tables <- vector(mode = "list", length = length(focal_sps))
names(coef_tables) <- focal_sps

# Create list to store interaction coefficients
sps_int <- vector(mode = "list", length = length(focal_sps))
names(sps_int) <- focal_sps

# Loop through training sets
for(set in train_sets){
  
  # Load training data
  #train <- read.csv(paste("Data/Output_data/training", set, ".csv", sep = ""),
  #                 stringsAsFactors = F)
  train <- read.csv(paste("Data/Output_data/rand_training", set, ".csv",
                          sep = ""), stringsAsFactors = F)
  
  # Load test data
  #test <- read.csv(paste("Data/Output_data/test", set, ".csv", sep = ""),
  #                 stringsAsFactors = F)
  test <- read.csv(paste("Data/Output_data/rand_test", set, ".csv", sep = ""),
                   stringsAsFactors = F)
  
  # Remove unneeded columns from training and test
  train <- train %>%
    select(tree_id, species, sps_comp, dbh_comp, prox, all_density,
           ABAM_density, ABPR_density, CANO_density, PSME_density, TABR_density,
           THPL_density, TSHE_density, ABLA_density, TSME_density, ALSI_density,
           ALVI_density, PIMO_density, ALRU_density, PICO_density, PIEN_density,
           POBA_density, ABGR_density, PISI_density, pet_mm, size_corr_growth)
  test <- test %>%
    select(tree_id, species, sps_comp, dbh_comp, prox, all_density,
           ABAM_density, ABPR_density, CANO_density, PSME_density, TABR_density,
           THPL_density, TSHE_density, ABLA_density, TSME_density, ALSI_density,
           ALVI_density, PIMO_density, ALRU_density, PICO_density, PIEN_density,
           POBA_density, ABGR_density, PISI_density, pet_mm, size_corr_growth)
  
  # Loop through focal species
  for(i in 1:length(focal_sps)){
    
    # Subset to focal species
    sing_sp <- train %>%
      filter(species == focal_sps[i]) # REMOVE SPECIES COL
    ss_test <- test %>%
      filter(species == focal_sps[i])
    
    # Load common competitors data and extract for focal species
    comm_comp <- read.csv("Data/Output_data/common_comps.csv",
                          stringsAsFactors = F)
    comm_comp <- comm_comp[, focal_sps[i]]
    
    # Get rare competitor density columns
    rare_dens <- setdiff(names(sing_sp)[grep("density", names(sing_sp))],
                         c(paste(comm_comp, "density", sep = "_"),
                           "all_density"))
    
    # Convert rare competitors to OTHR and remove unneeded densities
    sing_sp <- sing_sp %>%
      mutate(sps_comp = if_else(sps_comp %in% comm_comp, sps_comp, "OTHR"),
             OTHR_density = apply(sing_sp %>% select(all_of(rare_dens)), 1,
                                  sum)) %>%
      select(-all_of(rare_dens))
    ss_test <- ss_test %>%
      mutate(sps_comp = if_else(sps_comp %in% comm_comp, sps_comp, "OTHR"),
             OTHR_density = apply(ss_test %>% select(all_of(rare_dens)), 1,
                                  sum)) %>%
      select(-all_of(rare_dens))
    
    # Run regularized regression model
    set.seed(100)
    mod <- growth_model(sing_sp, "size_corr_growth", focal_sps[i], # REMOVE FOCAL SPECIES ARGUMENT
                        iterations = 100, test = ss_test)
    
    # Store training R squared value
    rsq_vals[set + ((i - 1) * length(train_sets)), 1] <- mod$R_squared
    
    # Store test R squared value
    rsq_vals[set + ((i - 1) * length(train_sets)), 2] <- mod$test_R_squared
    
    # Extract best model coefficients and interaction coefficients
    mod_coef_all <- as.matrix(coef(mod$mod))[c(1, 3:nrow(coef(mod$mod))), ]
    mod_coef_all <- as.data.frame(t(mod_coef_all))
    mod_sps_int <- apply(mod$mod_coef[, grep("sps_comp", names(mod$mod_coef))],
                         2, rescale)
    mod_sps_int <- as.data.frame(t(mod_sps_int))
    
    # Store coefficients
    if(set == 1){
      coef_tables[[focal_sps[i]]] <- mod_coef_all
      sps_int[[focal_sps[i]]] <- mod_sps_int
    } else {
      coef_tables[[focal_sps[i]]] <- bind_rows(coef_tables[[focal_sps[i]]],
                                               mod_coef_all)
      sps_int[[focal_sps[i]]] <- bind_rows(sps_int[[focal_sps[i]]], mod_sps_int)
    }
    
    # Subset coefficients table to neighborhood coefficients
    coef_sub <- mod$mod_coef %>%
      select(-c("(Intercept)", "mse", "R_squared", "pet_mm"))
    
    # Record number of models where neighborhood matters
    nbhd_inf[i, set] <- sum(apply(abs(coef_sub), 1, sum) != 0)
    
    # Further subset to only species identity variables
    coef_sub <- coef_sub %>%
      select(-c("dbh_comp", "prox", "all_density"))
    
    # Record number of models where neighbor species identity matters
    comp_id_inf[i, set] <- sum(apply(abs(coef_sub), 1, sum) != 0)
    
    # Further subset to conspecific coefficients
    coef_sub <- coef_sub %>%
      select(c(paste("sps_comp", focal_sps[i], sep = ""),
               paste(focal_sps[i], "density", sep = "_")))
    
    # Record number of models where interspecific competition stronger
    inter_str[i, set] <- length(which(apply(coef_sub, 1, sum) > 0 &
                                        coef_sub[, 1] >= 0 &
                                        coef_sub[, 2] >= 0))
    
    # Record number of models where intraspecific competition stronger
    intra_str[i, set] <- length(which(apply(coef_sub, 1, sum) < 0 &
                                        coef_sub[, 1] <= 0 &
                                        coef_sub[, 2] <= 0))
    
  }
}

# Format and save R squared values table
rsq_vals <- as.data.frame(rsq_vals)
rsq_vals <- rsq_vals %>%
  mutate(species = rep(focal_sps, each = length(train_sets)),
         training = rep(train_sets, times = length(focal_sps))) %>%
  rename(train_r2 = V1, test_r2 = V2) %>%
  select(species, training, train_r2, test_r2)
#write.csv(rsq_vals, "Data/Figure_data/RR_r2.csv", row.names = F)
write.csv(rsq_vals, "Data/Figure_data/RR_r2_rand.csv", row.names = F)

# Format and save coefficient tables
for(i in 1:length(focal_sps)){
  coef_tab <- as.data.frame(t(coef_tables[[focal_sps[i]]]))
  coef_tab <- coef_tab[2:nrow(coef_tab),]
  names(coef_tab) <- train_sets
  #write.csv(coef_tab, paste("Data/Figure_data/RR_coef_", focal_sps[i],
  #                          ".csv", sep = ""))
  write.csv(coef_tab, paste("Data/Figure_data/RR_coef_rand_", focal_sps[i],
                            ".csv", sep = ""))
}

# Format and save species interaction coefficients
int_coef <- data.frame(sps_compABAM = double())
for(i in 1:length(focal_sps)){
  int_coef <- bind_rows(int_coef, as.data.frame(sps_int[[focal_sps[i]]]))
}
id_cols <- data.frame(species = rep(focal_sps, each = length(train_sets)),
                      training = rep(train_sets, times = length(focal_sps)))
int_coef <- cbind(id_cols, int_coef)
#write.csv(int_coef, "Data/Figure_data/RR_sps_ints.csv", row.names = F)
write.csv(int_coef, "Data/Figure_data/RR_sps_ints_rand.csv", row.names = F)

# Format and save other ecological interpretation tables
nbhd_inf <- as.data.frame(nbhd_inf)
nbhd_inf <- nbhd_inf %>%
  mutate(species = focal_sps) %>%
  rename(train1 = V1, train2 = V2, train3 = V3, train4 = V4) %>%
  select(species, train1, train2, train3, train4)
#write.csv(nbhd_inf, "Data/Figure_data/nbhd_inf.csv", row.names = F)
write.csv(nbhd_inf, "Data/Figure_data/nbhd_inf_rand.csv", row.names = F)

comp_id_inf <- as.data.frame(comp_id_inf)
comp_id_inf <- comp_id_inf %>%
  mutate(species = focal_sps) %>%
  rename(train1 = V1, train2 = V2, train3 = V3, train4 = V4) %>%
  select(species, train1, train2, train3, train4)
#write.csv(comp_id_inf, "Data/Figure_data/comp_id_inf.csv", row.names = F)
write.csv(comp_id_inf, "Data/Figure_data/comp_id_inf_rand.csv", row.names = F)

inter_str <- as.data.frame(inter_str)
inter_str <- inter_str %>%
  mutate(species = focal_sps) %>%
  rename(train1 = V1, train2 = V2, train3 = V3, train4 = V4) %>%
  select(species, train1, train2, train3, train4)
#write.csv(inter_str, "Data/Figure_data/inter_str.csv", row.names = F)
write.csv(inter_str, "Data/Figure_data/inter_str_rand.csv", row.names = F)

intra_str <- as.data.frame(intra_str)
intra_str <- intra_str %>%
  mutate(species = focal_sps) %>%
  rename(train1 = V1, train2 = V2, train3 = V3, train4 = V4) %>%
  select(species, train1, train2, train3, train4)
#write.csv(intra_str, "Data/Figure_data/intra_str.csv", row.names = F)
write.csv(intra_str, "Data/Figure_data/intra_str_rand.csv", row.names = F)

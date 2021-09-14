library(dplyr)

# Load neighborhoods data
nbhds <- read.csv("Data/Output_data/neighborhoods.csv")

# Load tree measurement data
tree <- read.csv("Data/Raw_data/tree_growth_2017.csv", stringsAsFactors = F)

# Extract most recent size for each tree
records <- tree %>%
  group_by(tree_id) %>%
  arrange(desc(year)) %>%
  summarize(dbh = dbh[1]) %>%
  right_join(nbhds %>%
               select(tree_id, stand_id, species) %>%
               distinct(),
             by = "tree_id")

# Calculate area at breast height in m^2 for each tree
records <- records %>%
  mutate(abh_m2 = pi * ((dbh / 100) / 2) ^ 2)

# Summarize by species and stand
records <- records %>%
  group_by(stand_id, species) %>%
  summarize(num_trees = n(),
            mean_dbh = mean(dbh),
            dens = sum(abh_m2))

# Add a total abh per stand column using self-join
records <- records %>%
  left_join(records %>%
              group_by(stand_id) %>%
              summarize(total_abh = sum(dens)),
            by = "stand_id")

# Calculate proportional densities
records <- records %>%
  mutate(prop_dens = dens / total_abh)

# Summarize for focal species
tableS2 <- records %>%
  filter(species %in% c("ABAM", "CANO", "PSME", "THPL", "TSHE", "TSME")) %>%
  group_by(species) %>%
  summarize(plots = n(),
            trees_mean = mean(num_trees),
            trees_sd = sd(num_trees),
            dbh_mean = mean(mean_dbh),
            dbh_sd = sd(mean_dbh),
            dens_mean = mean(dens),
            dens_sd = sd(dens),
            prop_dens_mean = mean(prop_dens),
            prop_dens_sd = sd(prop_dens))

# Save table to csv
write.csv(tableS2, "Figures/species_by_plot_summary.csv", row.names = F)

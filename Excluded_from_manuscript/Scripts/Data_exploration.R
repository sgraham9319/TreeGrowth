# Load required packages
library(dplyr)
library(ggplot2)

# Load abiotic data
abio <- read.csv("Data/stand_abiotic_data.csv")

# Create side of mountain variable
abio <- abio %>%
  mutate(side = if_else(stand_id %in% c("AB08", "PP17", "AO03", "AV02", "TA01"),
                        "east", "west"))

# Plot PET vs. elevation
ggplot(data = abio, aes(x = elev_m, y = pet_mm, col = side)) +
  geom_point() +
  geom_smooth(method = "lm", fill = NA)
# Elevation effect captured but not side effect

# Try AET instead
ggplot(data = abio, aes(x = elev_m, y = aet_mm, col = side)) +
  geom_point() +
  geom_smooth(method = "lm", fill = NA)

ggplot(data = abio, aes(x = temp_C, y = pet_mm, col = side)) +
  geom_point()


# Load neighborhood data
nbhd <- read.csv("Data/Output_data/neighborhoods.csv")

# Summarize by tree
nbhd <- nbhd %>%
  group_by(tree_id) %>%
  summarize(stand_id = stand_id[1],
            species = species[1],
            dbh = dbh[1],
            all_density = all_density[1],
            annual_growth = annual_growth[1],
            size_corr_growth = size_corr_growth[1],
            precip_mm = precip_mm[1],
            temp_C = temp_C[1],
            elev_m = elev_m[1],
            aet_mm = aet_mm[1],
            pet_mm = pet_mm[1])

# Specify focal species
focal_sps <- "CANO"

# Plot PET vs. growth relationship
sing_sp <- nbhd %>%
  filter(species == focal_sps)
ggplot(sing_sp, aes(x = pet_mm, y = annual_growth)) +
  geom_point()

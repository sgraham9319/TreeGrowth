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

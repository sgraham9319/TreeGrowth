library(dplyr)
library(ggplot2)

# Specify focal species and model structure
focal_sps <- "TSME"
model_str <- "no_comp"

set <- 1

# Load model fits
output <- read.csv(paste("Data/Output_data/", model_str, set, "_", focal_sps,
                         ".csv", sep = ""))

# Load training data
training <- read.csv(paste("Data/Output_data/training", set, ".csv", sep = ""),
                     stringsAsFactors = F)

# Subset to focal species and remove unneeded columns
sing_sp <- training %>%
  arrange(tree_id) %>%
  filter(species == focal_sps) %>%
  select(tree_id, stand_id, species, dbh, prox, sps_comp, dbh_comp,
         annual_growth, pet_mm)

# Change units of variables to give them similar ranges
sing_sp <- sing_sp %>%
  mutate(
    dbh = dbh / 10,                     # cm to dm
    dbh_comp = dbh_comp / 10,           # cm to dm
    annual_growth = annual_growth * 10, # cm to mm
    pet_dm = pet_mm / 100               # mm to dm
  ) %>%
  select(-pet_mm)

# Extract annual growth of each focal individual
focals <- sing_sp %>%
  group_by(tree_id) %>%
  summarize(dbh = dbh[1], annual_growth = annual_growth[1],
            pet_dm = pet_dm[1])

# Calculate AICc for each model
AICc_calc <- function(k, NLL, n){
  (2 * (k + NLL)) + ((2 * ((k^2) + k)) / (n - k - 1))
}
output$AICc <- AICc_calc(length(grep("_opt", names(output))), output$NLL, nrow(focals))

# Calculate dAICc for each model
output$dAICc <- output$AICc - min(output$AICc)

# Order output by dAICc
output <- output %>% arrange(dAICc)

ggplot(data.frame(x=c(-2, 2)), aes(x)) +
  stat_function(fun=function(x) x ^ 2, n = 10000) +
  stat_function(fun=function(x) 0.5 * x ^ 2, n = 10000)


x <- ggplot(data.frame(x=c(-2, 2)), aes(x)) +
  stat_function(fun=function(x) x ^ 2, n = 10000)
x <- x +
  stat_function(fun=function(x) 0.5 * x ^ 2, xlim = c(-1.5, 1.5), n = 10000)
x


# Create example output table
plot_data <- data.frame(
  set = 1:4,
  model_str = c("no_comp", "eq_comp", "eq_comp", "no_comp"),
  gmax = c(4, 3.5, 4.2, 4.5),
  X0 = c(8, 10, 20, 11),
  Xb = c(0.8, 1, 2, 0.9),
  pet_a = c(4.5, 4.7, 4.3, 4.4),
  pet_b = c(2.2, 2.3, 2.5, 2.3),
  C = c(NA, 0.1, 0.05, NA),
  max_dbh = c(20, 25, 19, 30),
  mean_nci = c(NA, 7, 4, NA),
  mean_pet = c(3.9, 3.8, 3.95, 4))

nc_fit <- function(x, gmax, X0, Xb, mean_pet, pet_a, pet_b){
  gmax * exp((-0.5) * ((log(x / X0) / Xb) ^ 2)) *
    exp((-0.5) * (((mean_pet - pet_a) / pet_b) ^ 2))
}
comp_fit <- function(x, gmax, X0, Xb, mean_pet, pet_a, pet_b, C, mean_nci){
  gmax * exp((-0.5) * ((log(x / X0) / Xb) ^ 2)) *
    exp((-0.5) * (((mean_pet - pet_a) / pet_b) ^ 2)) *
    exp(-C * mean_nci)
}

cols <- c("red", "blue", "gold", "green")
plot_obj <- ggplot(data.frame(x = c(0, max(plot_data$max_dbh))), aes(x)) +
  labs(title = focal_sps, x = "Diameter at breast height (dm)",
       y = "Annual growth (mm/year)")

for(i in 1:4){
  
  if(plot_data$model_str[i] == "no_comp"){
    plot_obj <- plot_obj +
      stat_function(fun = nc_fit, xlim = c(0.1, plot_data$max_dbh[i]),
                    col = cols[i],
                    args = list(gmax = plot_data$gmax[i],
                                X0 = plot_data$X0[i],
                                Xb = plot_data$Xb[i],
                                mean_pet = plot_data$mean_pet[i],
                                pet_a = plot_data$pet_a[i],
                                pet_b = plot_data$pet_b[i]))
  } else {
    plot_obj <- plot_obj +
      stat_function(fun = comp_fit, xlim = c(0.1, plot_data$max_dbh[i]),
                    col = cols[i],
                    args = list(gmax = plot_data$gmax[i],
                                X0 = plot_data$X0[i],
                                Xb = plot_data$Xb[i],
                                mean_pet = plot_data$mean_pet[i],
                                pet_a = plot_data$pet_a[i],
                                pet_b = plot_data$pet_b[i],
                                C = plot_data$C[i],
                                mean_nci = plot_data$mean_nci[i]))
  }
}

plot_obj

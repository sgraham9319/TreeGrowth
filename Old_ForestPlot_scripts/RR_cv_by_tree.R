library(dplyr)

# Define required functions
z_trans <- function(x){(x - mean(x)) / sd(x)}
coef_det <- function(x){
  1 - (sum((x$observations - x$predictions)^2) / 
         sum((x$observations - mean(x$observations))^2))
}

# Define focal species
focal_sps <- "ABAM"

# Load training data
train <- read.csv(paste("Data/Output_data/training1.csv", sep = ""),
                  stringsAsFactors = F)

# Subset to focal species
train <- train %>%
  filter(species == focal_sps)

# Create vector of common competitors
comm_comps <- names(which(table(train$sps_comp) > 100))

# Reduce to one row per tree id
train <- train %>%
  group_by(tree_id) %>%
  filter(row_number() == 1) %>%
  ungroup()

# Extract observations
obs <- train %>%
  select(tree_id, size_corr_growth)

# Create column for density of rare species
train <- train %>%
  mutate(rare_density = abs(
    all_density - apply(train %>%
                          select(paste(comm_comps, "density", sep = "_")),
                        1, sum)))

# Select columns required by model
train <- train %>%
  select(size_corr_growth, all_density,
         paste(comm_comps, "density", sep = "_"), rare_density, 
         precip_mm, temp_C, elev_m, aet_mm, pet_mm)

# Remove any columns containing only zeros
train <- train %>%
  select(-which(apply(train, 2, sum) == 0))

# Define outcome variable
outcome_var <- "size_corr_growth"

# Create model formula object
mod_form <- as.formula(paste(outcome_var, "~",
                             paste(setdiff(names(train), outcome_var),
                                   collapse = "+")))

# Create design matrix
dm <- model.matrix(mod_form, train)

# Standardize variables except for first column (intercept)
dm[, 2:ncol(dm)] <- apply(dm[, 2:ncol(dm)], 2, z_trans)

# Change any columns of NaNs (no variation) to zeros
dm[, which(is.nan(dm[1, ]))] <- 0

# Fit glmnet model
mod <- glmnet::cv.glmnet(x = dm, y = train$size_corr_growth,
                         family = "gaussian")

# Plot lambda vs. MSE
plot(mod)

# Make predictions for training data
preds <- predict(mod, newx = dm, s = "lambda.1se")

# Combine predictions with observations
obs_pred <- cbind(obs, preds)
names(obs_pred) <- c("tree_id", "observations", "predictions")

# Calculate coefficient of determination
R_squared <- coef_det(obs_pred)

# Load test data
test <- read.csv(paste("Data/Output_data/test1.csv", sep = ""),
                 stringsAsFactors = F)

# Reduce to one row per tree id and subset to focal species
test <- test %>%
  group_by(tree_id) %>%
  filter(row_number() == 1 & species == focal_sps) %>%
  ungroup()

# Extract observations
obs_test <- test %>%
  select(tree_id, size_corr_growth)

# Create column for density of rare species
test <- test %>%
  mutate(rare_density = abs(
    all_density - apply(test %>%
                          select(paste(comm_comps, "density", sep = "_")),
                        1, sum)))

# Select columns required by model
test <- test %>%
  select(size_corr_growth, all_density,
         paste(comm_comps, "density", sep = "_"), rare_density, 
         precip_mm, temp_C, elev_m, aet_mm, pet_mm)

# Remove any columns containing only zeros
test <- test %>%
  select(-which(apply(test, 2, sum) == 0))

# Create test design matrix
dm_test <- model.matrix(mod_form, test)

# Standardize variables except for first column (intercept)
dm_test[, 2:ncol(dm_test)] <- apply(dm_test[, 2:ncol(dm_test)], 2, z_trans)

# Change any columns of NaNs (no variation) to zeros
dm_test[, which(is.nan(dm_test[1, ]))] <- 0

# Make predictions for test data
preds_test <- predict(mod, newx = dm_test, s = "lambda.1se")

# Combine predictions with observations
obs_pred_test <- cbind(obs_test, preds_test)
names(obs_pred_test) <- c("tree_id", "observations", "predictions")

# Calculate coefficient of determination
R_squared_test <- coef_det(obs_pred_test)

# Check coefficients of 1se model
coef(mod)

size_effect_comp <- function(fit_data, species, cols){
  
  # Define no competition size effect function
  nc_fit <- function(x, gmax, X0, Xb, mean_pet, pet_a, pet_b){
    gmax * exp((-0.5) * ((log(x / X0) / Xb) ^ 2)) *
      exp((-0.5) * (((mean_pet - pet_a) / pet_b) ^ 2))
  }
  
  # Define competition size effect function
  comp_fit <- function(x, gmax, X0, Xb, mean_pet, pet_a, pet_b, C, mean_nci){
    gmax * exp((-0.5) * ((log(x / X0) / Xb) ^ 2)) *
      exp((-0.5) * (((mean_pet - pet_a) / pet_b) ^ 2)) *
      exp(-C * mean_nci)
  }
  
  # Generate xy values for each model fit
  plot_format <- lapply(1:nrow(fit_data), function(i){
    
    if(fit_data$model_str[i] == "no_comp"){
      args <- fit_data[i, 3:8]
      data.frame(training = fit_data$training[i],
                 x = seq(0.1, fit_data$max_dbh[i], 0.1),
                 y = do.call(nc_fit,
                             c(list(x = seq(0.1, fit_data$max_dbh[i], 0.1)),
                               args)))
    } else {
      args <- fit_data[i, 3:10]
      data.frame(training = fit_data$training[i],
                 x = seq(0.1, fit_data$max_dbh[i], 0.1),
                 y = do.call(comp_fit,
                             c(list(x = seq(0.1, fit_data$max_dbh[i], 0.1)),
                               args)))
    }
  }) %>% bind_rows
  
  # Convert training column to factor
  plot_format$training <- as.factor(plot_format$training)
  
  # Create plot
  ggplot(plot_format, aes(x = x, y = y, color = training)) +
    geom_line(lwd = 1) +
    scale_color_manual(values = cols) +
    labs(title = species, x = "Diameter at breast height (dm)",
         y = "Annual growth (mm/year)") 
}
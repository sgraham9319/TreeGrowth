
r2_figure <- function(r2_table, name){
  
  # Create empty data frame for formatted data
  plot_data <- data.frame(
    focal_sps = character(),
    mean = double(),
    max = double(),
    min = double()
  )
  
  # Format data for plotting
  for(column in c("train_lkhd_AIC", "train_lkhd_cv", "train_rr",
                  "test_lkhd_AIC", "test_lkhd_cv", "test_rr")){
    new_data <- r2_table %>%
      group_by(focal_sps) %>%
      summarize(mean = mean(get(column), na.rm = T),
                max = max(get(column), na.rm = T),
                min = min(get(column), na.rm = T))
    plot_data <- bind_rows(plot_data, new_data)
  }
  
  # Add columns to identify models
  plot_data <- plot_data %>%
    mutate(model = rep(rep(c("lkhd_AIC", "lkhd_cv", "rr"), each = 6),
                       times = 2),
           train_test = rep(c("train", "test"), each = 18)) %>%
    select(train_test, model, focal_sps, mean, max, min)
  
  # Add x coordinates and colors
  plot_data <- plot_data %>%
    mutate(x = rep(c(seq(1:6) - 0.2, seq(1:6), seq(1:6) + 0.2), times = 2),
           cols = rep(rep(c("blue", "green", "red"), each = 6), times = 2))
  
  # Separate training and test data
  training <- plot_data %>%
    filter(train_test == "train")
  test <- plot_data %>%
    filter(train_test == "test")
  
  # Initiate plot saving
  png(filename = paste("Figures/", name, ".png", sep = ""),
      type = "cairo",
      units = "in",
      width = 4,
      height = 5,
      pointsize = 12,
      res = 500)
  
  # Set plotting window parameters
  par(mfcol = c(2, 1),
      mar = c(0,0,0,0),
      oma = c(4,4,0.1,0.3),
      mgp = c(3, 0.8, 0))
  
  # Create training plot
  plot(training$x, training$mean, pch = 19, cex = 0.7,
       col = training$cols, las = 1,
       ylim = c(min(min(training$min), 0), 0.8),
       yaxt = "n", xaxt = "n")
  axis(side = 2, at = seq(0, 0.8, 0.4), labels = seq(0, 0.8, 0.4),
       tck = -0.02, las = 1, cex.axis = 0.9, mgp = c(3, 0.5, 0))
  arrows(training$x, training$min, training$x, 
         training$max, code = 3, angle = 90, length = 0, col = training$cols)
  abline(h = 0, lty = "dashed")
  text(x = 0.82, y = 0.78, label = "A")
  
  # Create test plot
  plot(test$x, test$mean, pch = 19, cex = 0.7,
       col = test$cols, las = 1,
       ylim = c(min(min(test$min), 0), 0.9),
       yaxt = "n", xaxt = "n")
  axis(side = 2, at = seq(0, 0.8, 0.4), labels = seq(0, 0.8, 0.4),
       tck = -0.02, las = 1, cex.axis = 0.9, mgp = c(3, 0.5, 0))
  arrows(test$x, test$min, test$x, 
         test$max, code = 3, angle = 90, length = 0, col = test$cols)
  abline(h = 0, lty = "dashed")
  text(x = 0.82, y = 0.86, label = "B")
  mtext(text = "Focal species", side = 1, line = 1.5)
  mtext(text = "Coefficient of determination", side = 2, line = 2, outer = T)
  axis(side = 1, at = 1:6,
       labels = c("ABAM", "CANO", "PSME", "THPL", "TSHE", "TSME"),
       cex.axis = 0.8, tck = -0.02, padj = -1)
  
  # Save plot
  dev.off()
}

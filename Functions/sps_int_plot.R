
sps_int_plot <- function(focal_sps, train_type){
  
  # Define AICc function
  AICc_calc <- function(k, NLL, n){
    (2 * (k + NLL)) + ((2 * ((k^2) + k)) / (n - k - 1))
  }
  
  #===================================
  # Extracting likelihood interactions
  #===================================
  
  # Load lists of common competitors
  comm_comps <- read.csv("Data/Output_data/common_comps.csv")
  
  # Extract common competitors for focal species
  comps <- sort(c(comm_comps[focal_sps][!is.na(comm_comps[focal_sps])], "OTHR"))
  
  # Create empty matrix to store num interactions with each competitor
  num_ints <- matrix(NA, ncol = length(comps), nrow = 4)
  
  # Loop through training sets
  for(set in 1:4){
    
    # Load training data
    if(train_type == "regular"){
      training <- read.csv(paste("Data/Output_data/training", set, ".csv",
                                 sep = ""), stringsAsFactors = F)
    } else if(train_type == "random"){
      training <- read.csv(paste("Data/Output_data/rand_training", set, ".csv",
                                 sep = ""), stringsAsFactors = F)
    }
    
    # Subset to focal species and change rare competitors to OTHR
    training <- training %>%
      filter(species == focal_sps) %>%
      mutate(sps_comp = if_else(sps_comp %in% comps, sps_comp, "OTHR"))
    
    # Calculate number of focal trees
    focals <- length(unique(training$tree_id))
    
    # Calculate number of interactions with each competitor species
    comp_n <- table(training$sps_comp)
    num_ints[set, ] <- as.vector(comp_n[comps])
    
    # Load model output
    if(train_type == "regular"){
      output <- read.csv(paste("Data/Output_data/ss_comp", set, "_", focal_sps,
                               ".csv", sep = ""))
    } else if(train_type == "random"){
      output <- read.csv(paste("Data/Output_data/ss_comp_rand", set, "_", focal_sps,
                               ".csv", sep = ""))
    }
    
    # Calculate AICc for each model
    output$AICc <- AICc_calc(length(grep("_opt", names(output))), output$NLL,
                             focals)
    
    # Arrange in order of AICc
    output <- output %>%
      arrange(AICc)
    
    # Extract fitted interaction coefficients
    if(set == 1){
      lkhd_int <- output[1, grep("lmd[0-9]_opt", names(output))]
    } else{
      lkhd_int <- rbind(lkhd_int, output[1, grep("lmd[0-9]_opt", names(output))])
    }
  }
  
  # Calculate average num interactions with each competitor species
  av_ints <- round(apply(num_ints, 2, mean), 0)
  
  # Add competitor species names to interactions table
  names(lkhd_int) <- comps
  
  #==========================================================
  # Including regularized regression interaction coefficients
  #==========================================================
  
  # Load regularized regression interactions
  if(train_type == "regular"){
    rr_int <- read.csv("Data/Figure_data/RR_sps_ints.csv")
  } else if(train_type == "random"){
    rr_int <- read.csv("Data/Figure_data/RR_sps_ints_rand.csv")
  }
  
  # Subset to focal species and select columns
  rr_int <- rr_int %>%
    filter(species == focal_sps) %>%
    select(paste("sps_comp", comps, sep = ""))
  
  # Create single vector of all interaction coefficients
  all_int <- NULL
  for(i in 1:length(comps)){
    all_int <- c(all_int, lkhd_int[, i], rr_int[, i])
  }
  
  # Flip sizes of interaction coefficients to represent focal growth
  all_int <- 1 - all_int
  
  # Reorder so "other" is last
  other_index <- which(comps == "OTHR")
  other_coef <- (other_index * 8 - 7):(other_index * 8)
  all_int <- c(all_int[-other_coef], all_int[other_coef])
  
  #============
  # Making plot
  #============
  
  # Define vertical separator between species on plot
  separator <- 0.02
  
  # Define height of each rectangle
  perfect_height <- 0.5 * (1 - (separator * nrow(comm_comps))) /
    (nrow(comm_comps) + 1)
  rect_ht <- trunc(perfect_height * 1000) / 1000
  
  # Calculate between block distance
  block <- separator + 2 * rect_ht
  
  # Create vector of rectangle base coordinates
  constant_bases <- rep(rep(c(0, rect_ht), each = 4), times = length(comps))
  modifiers <- rep(block * 0:(length(comps) - 1), each = 8)
  base_coords <- constant_bases + modifiers
  
  # Define coordinates of rectangles
  y_bottom <- rev(base_coords)
  y_top <- y_bottom + rect_ht
  x_left <- rep(seq(0, 0.75, 0.25), times = length(comps) * 2)
  x_right <- rep(seq(0.25, 1, 0.25), times = length(comps) * 2)
  
  # Create colors
  n_colors <- 20
  #colors <- paletteer_c(palette = "grDevices::heat.colors", n = n_colors)
  #colors <- paletteer_c(palette = "ggthemes::Red", n = n_colors)
  #colors <- paletteer_c(palette = "ggthemes::Temperature Diverging", n = n_colors)
  colors <- paletteer_c(palette = "pals::coolwarm", n = n_colors)
  # TO CHANGE PALETTE DIRECTION INCLUDE: direction = -1
  #colors <- paletteer_c(palette = "oompaBase::redscale", n = n_colors)
  #colors <- paletteer_c(palette = "harrypotter::ronweasley", n = n_colors)
  
  # Assign each coefficient a color
  rank <- cut(c(all_int, 0, 1), n_colors)
  
  # Set plot saving parameters
  if(train_type == "regular"){
    jpeg(filename = paste("Figures/Sps_int", focal_sps, ".jpg", sep = ""),
        type = "cairo",
        units = "px", 
        width = 1250, 
        height = 1000, 
        pointsize = 12, 
        res = 96)
  } else if(train_type == "random"){
    jpeg(filename = paste("Figures/Sps_int_rand", focal_sps, ".jpg", sep = ""),
        type = "cairo",
        units = "px", 
        width = 1250, 
        height = 1000, 
        pointsize = 12, 
        res = 96)
  }
  
  # Change plot margins
  par(mar = c(2, 5, 0, 4))
  
  # Create space for legend
  layout(matrix(1:2, nrow = 2), width = c(1, 1), height = c(5, 1))
  
  # Create empty plot
  plot.new()
  
  # Add rectangles
  rect(x_left, y_bottom, x_right, y_top, col = colors[rank[0:length(all_int)]],
       lwd = 2)
  
  # Reorder comps and av_ints so "other" is at the bottom
  new_comps <- c(comps[-other_index], "OTHER")
  new_av_ints <- c(av_ints[-other_index], av_ints[other_index])
  
  # Add commas to large numbers
  new_av_ints <- formatC(new_av_ints, format="d", big.mark=",")
  
  # Add axis labels
  mtext(text = 1:4, side = 1, at = seq(0.125, 0.875, 0.25), line = -0.2,
        cex = 2.5)
  mtext(text = "Training set", side = 1, line = 1.5, cex = 2.5)
  mtext(text = new_comps, side = 2, line = -0.8,
        at = rev(seq(rect_ht, block * length(comps) - (rect_ht + separator),
                     block)), las = 1, cex = 2)
  mtext(text = new_av_ints, side = 4, line = -1.1,
        at = rev(seq(rect_ht, block * length(comps) - (rect_ht + separator),
                     block)), las = 1, cex = 2)
  
  # Create and add legend
  par(mar = c(1, 5, 0, 4))
  plot.new()
  legend_image <- as.raster(matrix(colors, nrow = 1))
  rasterImage(legend_image, xleft = 0.25, ybottom = 0.25,
              xright = 0.75, ytop = 0.75)
  mtext(text = c(paste("Low", focal_sps, "Growth", sep = "\n"),
                 paste("High", focal_sps, "Growth", sep = "\n")),
        side = 1, at = c(0.15, 0.85), line = -1.5, cex = 2.5)
  
  # Output plot
  dev.off()
}
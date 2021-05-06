library(paletteer)
library(dplyr)

# Define coordinates of rectangles
x_left <- rep(seq(0, 0.75, 0.25), times = 8)
x_right <- rep(seq(0.25, 1, 0.25), times = 8)
y_bottom <- rep(seq(0, 0.875, 0.125), each = 4)
y_top <- rep(seq(0.125, 1, 0.125), each = 4)

# Define interaction coefficients
int_coef <- rep(seq(0, 0.9, 0.3), times = 8)

# Create colors
n_colors <- 20
#colors <- paletteer_c(palette = "grDevices::heat.colors", n = n_colors)
#colors <- paletteer_c(palette = "ggthemes::Red", n = n_colors)
#colors <- paletteer_c(palette = "ggthemes::Temperature Diverging", n = n_colors)
colors <- paletteer_c(palette = "pals::coolwarm", n = n_colors)
#colors <- paletteer_c(palette = "oompaBase::redscale", n = n_colors)
#colors <- paletteer_c(palette = "harrypotter::ronweasley", n = n_colors)

# Assign each coefficient a color
rank <- cut(c(int_coef, 0, 1), n_colors)

# Create empty plot
plot.new()

# Add rectangles
rect(x_left, y_bottom, x_right, y_top, col = colors[rank[0:length(int_coef)]])

#===================================
# Extracting likelihood interactions
#===================================

# Load lists of common competitors
comm_comps <- read.csv("Data/Output_data/common_comps.csv")

# Define AICc function
AICc_calc <- function(k, NLL, n){
  (2 * (k + NLL)) + ((2 * ((k^2) + k)) / (n - k - 1))
}

# Specify focal species
focal_sps <- "ABAM"

for(set in 1:4){
  
  # Load training data
  training <- read.csv(paste("Data/Output_data/training", set, ".csv",
                             sep = ""), stringsAsFactors = F)
  
  # Reduce to species identity of each focal
  training <- training %>%
    group_by(tree_id) %>%
    summarize(species = species[1])
  
  # Extract number of focal trees per species
  focals <- sum(training$species == focal_sps)
  
  # Load model output
  output <- read.csv(paste("Data/Output_data/ss_comp", set, "_", focal_sps,
                           ".csv", sep = ""))
  
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

# Add competitor species names to interactions table
for(sps in focal_sps){
  comps <- sort(c(comm_comps[focal_sps][!is.na(comm_comps[focal_sps])], "OTHR"))
  names(lkhd_int) <- comps
}



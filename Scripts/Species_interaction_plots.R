library(paletteer)
library(dplyr)

# Specify focal species
focal_sps <- "TSHE"

#===================================
# Extracting likelihood interactions
#===================================

# Load lists of common competitors
comm_comps <- read.csv("Data/Output_data/common_comps.csv")

# Define AICc function
AICc_calc <- function(k, NLL, n){
  (2 * (k + NLL)) + ((2 * ((k^2) + k)) / (n - k - 1))
}

# Loop through training sets
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
comps <- sort(c(comm_comps[focal_sps][!is.na(comm_comps[focal_sps])], "OTHR"))
names(lkhd_int) <- comps

#==========================================================
# Including regularized regression interaction coefficients
#==========================================================

# Load regularized regression interactions
rr_int <- read.csv("Data/Figure_data/RR_sps_ints.csv")

# Subset to focal species and select columns
rr_int <- rr_int %>%
  filter(species == focal_sps) %>%
  select(paste("sps_comp", comps, sep = ""))

# Create single vector of all interaction coefficients
all_int <- NULL
for(i in 1:length(comps)){
  all_int <- c(all_int, lkhd_int[, i], rr_int[, i])
}

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
#colors <- paletteer_c(palette = "oompaBase::redscale", n = n_colors)
#colors <- paletteer_c(palette = "harrypotter::ronweasley", n = n_colors)

# Assign each coefficient a color
rank <- cut(c(all_int, 0, 1), n_colors)

# Change plot margins
par(mar = c(2, 2, 2, 4))

# Create empty plot
plot.new()

# Add rectangles
rect(x_left, y_bottom, x_right, y_top, col = colors[rank[0:length(all_int)]])

# Add axis labels
mtext(text = 1:4, side = 1, at = seq(0.125, 0.875, 0.25), line = -0.2)
mtext(text = "Training set", side = 1, line = 1)
mtext(text = comps, side = 2,
      at = rev(seq(rect_ht, block * length(comps) - (rect_ht + separator),
                   block)), las = 1, cex = 0.8)

library(paletteer)

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

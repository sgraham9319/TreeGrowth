

# Define coordinates of rectangles
x_left <- rep(seq(0, 0.75, 0.25), times = 8)
x_right <- rep(seq(0.25, 1, 0.25), times = 8)
y_bottom <- rep(seq(0, 0.875, 0.125), each = 4)
y_top <- rep(seq(0.125, 1, 0.125), each = 4)

# Create empty plot
plot.new()

# Add rectangles
rect(x_left, y_bottom, x_right, y_top, col = c("blue", "red"))
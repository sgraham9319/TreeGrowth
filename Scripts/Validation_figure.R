library(dplyr)
library(plotfunctions)

# Define number of trees
n_trees <- 50

# Create x and y coordinates for trees
set.seed(100)
x_coords <- runif(n_trees, min = 0.15, max = 0.85)
y_coords <- runif(n_trees, min = 0.15, max = 0.85)

# Assign trees to test sets
test_set <- rep(1:4, length.out = n_trees)

# Combine coordinates and test set into data frame
trees <- data.frame(x_coords, y_coords, test_set)

# Create subsets for each test set
test1 <- trees %>% filter(test_set == 1)
test2 <- trees %>% filter(test_set == 2)
test3 <- trees %>% filter(test_set == 3)
test4 <- trees %>% filter(test_set == 4)

# Initiate plot saving
png(filename = "Figures/Validation_figure.png",
    type = "cairo",
    units = "in",
    width = 4,
    height = 4,
    pointsize = 12,
    res = 500)

# Adjust plotting margins
par(mar = c(2,0,0,0),
    oma = c(0.1, 0.1, 0.1, 0.1))

# Create empty plot
plot(NULL, ylim = c(0, 1), xlim = c(0, 1), yaxs = "i", xaxs = "i",
     yaxt = "n", ylab = "", xaxt = "n", xlab = "")

# Shade buffer zone
rect(xleft = 0, ybottom = 0, xright = 0.15, ytop = 1, border = NA, col = "grey")
rect(xleft = 0, ybottom = 0, xright = 1, ytop = 0.15, border = NA, col = "grey")
rect(xleft = 0, ybottom = 0.85, xright = 1, ytop = 1, border = NA, col = "grey")
rect(xleft = 0.85, ybottom = 0, xright = 1, ytop = 1, border = NA, col = "grey")

# Add trees from test set 1
points(x = test1$x_coords, y = test1$y_coords, pch = 21, bg = "blue", cex = 1.5)
points(x = test2$x_coords, y = test2$y_coords, pch = 21, bg = "red", cex = 1.5)
points(x = test3$x_coords, y = test3$y_coords, pch = 21, bg = "green", cex = 1.5)
points(x = test4$x_coords, y = test4$y_coords, pch = 21, bg = "yellow", cex = 1.5)

# Add covered plot border lines
lines(x = c(0, 1), y = c(1, 1))
lines(x = c(0, 0), y = c(0, 1))

# Add legend
legend_margin(x = "bottom", legend = c("Test 1", "Test 2", "Test 3", "Test 4"),
              bty = "n", horiz = T, pch = 21,
              pt.bg = c("blue", "red", "green", "yellow"), pt.cex = 1.5)

# Save plot
dev.off()
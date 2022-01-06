# Create example data frame
data <- data.frame(x1 = 1:5,
                   x2 = letters[1:5],
                   x3 = c(4, 1, 5, 3, 1),
                   x4 = c("Male", "Female", "Male", "Male", "Female"),
                   stringsAsFactors = FALSE)

data_1 <- factor(data$x4)
names(data_1) = rownames(data)

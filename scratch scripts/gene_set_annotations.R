library(ComplexHeatmap)

dev.off()

random_text = function(n) {
  sapply(1:n, function(i) {
    paste0(sample(letters, sample(4:10, 1)), collapse = "")
  })
}
text_list = list(
  text1 = random_text(4),
  text2 = random_text(4),
  text3 = random_text(4),
  text4 = random_text(4)
)

# note how we set the width of this empty annotation
ha = rowAnnotation(foo = anno_empty(border = FALSE, width = max_text_width(unlist(text_list)) + unit(4, "mm")))

Heatmap(matrix(rnorm(1000), nrow = 100), name = "mat", row_km = 4, right_annotation = ha)

for(i in 1:4) {
  decorate_annotation("foo", slice = i, {
    grid.rect(x = 0, width = unit(2, "mm"), gp = gpar(fill = i, col = NA), just = "left")
    grid.text(paste(text_list[[i]], collapse = "\n"), x = unit(4, "mm"), just = "left")
  })
}
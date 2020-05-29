
library(magick)
img1 <- image_read("figs/study site figure.tiff")
img2 <- image_read("figs/era-specific PDO and climate - vertical.tiff")
# img <- c(image_scale(img1, "80%"), image_scale(img2, "120%"))
img <- c(img1, img2)

stack <- image_append(image_scale(img, "100%"), stack = TRUE)
image_write(stack, path = "figs/combined Fig2.tiff", format = "tiff")

image_write(stack, path = "figs/combined Fig2.png", format = "png")

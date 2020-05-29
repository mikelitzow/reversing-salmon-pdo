library(png)
library(gridExtra)
library(gtable)
library(grid)
library(ggpubr)

#Loading images
img1 <- readPNG("figs/study site figure.png")
img2 <- readPNG("figs/era-specific PDO and climate updated.png")

#Convert images to Grob (graphical objects)
grob1 <- rasterGrob(img1)
grob2 <- rasterGrob(img2)

# create the plot with data
p <- ggarrange(scatter, int, ncol=2, nrow=1, labels=c("a)", "b)"),
               label.x = 0.05, label.y = 0.95)

# convert the plot to gtable
mytable <- ggplot_gtable(ggplot_build(p))

#Add a line ssame height as line 4 after line 3
# use gtable_show_layout(mytable) to see where you want to put your line
mytable <- gtable_add_rows(mytable, mytable$height[[4]], 3)
# if needed mor row can be added 
# (for example to keep consistent spacing with other graphs)

#Insert the grob in the cells of the new line
mytable <- gtable_add_grob(mytable,grob1,4,4, name = paste(runif(1)))
mytable <- gtable_add_grob(mytable,grob1,4,6, name = paste(runif(1)))

#rendering
grid.draw(mytable)

library(magick)
img1 <- image_read("figs/study site figure.tiff")
img2 <- image_read("figs/era-specific PDO and climate updated.tiff")
# img <- c(image_scale(img1, "80%"), image_scale(img2, "120%"))
img <- c(img1, img2)

stack <- image_append(image_scale(img, "100%"), stack = TRUE)
image_write(stack, path = "figs/combined Fig2.tiff", format = "tiff")

tiff("figs/combined Fig2.tiff", 8, 9.25, units='in', res=300)
image_append(image_scale(img, "100"), stack = TRUE)
dev.off()

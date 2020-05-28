library(png)
library(gridExtra)
library(ggplot2)
library(gtable)
library(RCurl) # just to load image from URL

#Loading images
img0 <- readPNG(getURLContent('http://i.stack.imgur.com/3anUH.png'))
img1 <- readPNG(getURLContent('http://i.stack.imgur.com/3anUH.png'))

#Convert images to Grob (graphical objects)
grob0 <- rasterGrob(img0)
grob1 <- rasterGrob(img1)

# create the plot with data
p <- ggplot(subset(mpg,manufacturer=='audi'|manufacturer=='toyota'),aes(x=cty,y=displ))+facet_grid(year~manufacturer)+geom_point()

# convert the plot to gtable
mytable <- ggplot_gtable(ggplot_build(p))

#Add a line ssame height as line 4 after line 3
# use gtable_show_layout(mytable) to see where you want to put your line
mytable <- gtable_add_rows(mytable, mytable$height[[4]], 3)
# if needed mor row can be added 
# (for example to keep consistent spacing with other graphs)

#Insert the grob in the cells of the new line
mytable <- gtable_add_grob(mytable,grob0,4,4, name = paste(runif(1)))
mytable <- gtable_add_grob(mytable,grob1,4,6, name = paste(runif(1)))

#rendering
grid.draw(mytable)
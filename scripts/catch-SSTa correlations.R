library(dplyr)
library(pracma)
library(ggplot2)
library(reshape2)
library(tidyr)
library(gtools)
library(nlme)
library(ggpubr)
library(ncdf4)
library(chron)
library(maps)
library(fields)
library(maptools)
library(mapdata)
library(zoo)
library(oce)

# plot correlation maps for salmon catch and gridded SST

raw.dat <- read.csv("data/salmon.and.covariate.data.csv")

# note that years in this file are already lagged to ocean entry!
# re-lag to define catch year!
raw.dat$catch.year <- ifelse(raw.dat$species=="Sockeye", raw.dat$Year+2, raw.dat$Year+1)

nc <- nc_open("data/North.Pacific.ersst")

# extract dates
ncvar_get(nc, "time")   # seconds since 1-1-1970
raw <- ncvar_get(nc, "time")
h <- raw/(24*60*60)
sst.d <- dates(h, origin = c(1,1,1970))

sst.x <- ncvar_get(nc, "longitude")
sst.y <- ncvar_get(nc, "latitude")

SST <- ncvar_get(nc,  "sst")
# Change data from a 3-D array to a matrix of monthly data by grid point:
# First, reverse order of dimensions ("transpose" array)
SST <- aperm(SST, 3:1)  

# Change to matrix with column for each grid point, rows for monthly means
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  

# Keep track of corresponding latitudes and longitudes of each column:
sst.lat <- rep(sst.y, length(sst.x))   
sst.lon <- rep(sst.x, each = length(sst.y))   
dimnames(SST) <- list(as.character(sst.d), paste("N", sst.lat, "E", sst.lon, sep=""))

# plot to check
SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=tim.colors(64))
contour(sst.x, sst.y, z, add=T) 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)

# a-ok!

# now get NDJFM means for each cell
m <- months(sst.d)
yr <- as.numeric(as.character(years(sst.d)))

win.months <- m[m %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]
win.yrs <- ifelse(m %in% c("Nov", "Dec"), yr+1, yr)
win.yrs <- win.yrs[m %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]

# get seperate matrix for winter
win.SST <- SST[m %in% win.months,]

f <- function(x) tapply(x, win.yrs, mean)
win.SST <- apply(win.SST, 2, f)

# and smooth each cell as 3-yr running mean
# I'm doing this b/c analysis shows that catch most strongly 
# responds to SST over 3-yr time scales
# also adding 2-yr rolling mean and retaining 1-yr for consistency with PDO

f3 <- function(x) rollmean(x, 3, fill=NA)
win.SST.3 <- apply(win.SST, 2, f3)

f2 <- function(x) rollmean(x, 2, align="right", fill=NA)
win.SST.2 <- apply(win.SST, 2, f2)

win.SST.1 <- win.SST # and the unsmoothed version!

# add rownames to smoothed version
rownames(win.SST.3) <- row.names(win.SST.2) <- row.names(win.SST.1)

# plot to check
SST.mean <- colMeans(win.SST.3, na.rm=T)
z <- t(matrix(SST.mean,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=tim.colors(64))
contour(sst.x, sst.y, z, add=T) 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)

# all good!
# now make correlation plots for each spp. in each era...

# first GOA sockeye
winter.sockeye.65.88 <- winter.sockeye.89.13 <- winter.sockeye.14.19 <- NA

for(j in 1:ncol(win.SST.3)){
  # note that we are using catch year as the nominal year matching the proposed eras
   # j <- 1
  winter.sockeye.65.88[j] <- 
    cor(win.SST.3[rownames(win.SST.3) %in% 1963:1986, j], 
        raw.dat$catch[raw.dat$species=="Sockeye" & 
                                      raw.dat$catch.year %in% 1965:1988]) # 2-yr lag!
  
  winter.sockeye.89.13[j] <- 
    cor(win.SST.3[rownames(win.SST.3) %in% 1987:2011, j], 
        raw.dat$catch[raw.dat$species=="Sockeye" & 
                        raw.dat$catch.year %in% 1989:2013])
  
  winter.sockeye.14.19[j] <- 
    cor(win.SST.3[rownames(win.SST.3) %in% 2012:2017, j], 
        raw.dat$catch[raw.dat$species=="Sockeye" & 
                        raw.dat$catch.year %in% 2014:2019])
}

# now coho!
winter.coho.65.88 <- winter.coho.89.13 <- winter.coho.14.19 <- NA

for(j in 1:ncol(win.SST.3)){
  # note that we are using catch year as the nominal year matching the proposed eras
  # j <- 1
  winter.coho.65.88[j] <- 
    cor(win.SST.3[rownames(win.SST.3) %in% 1964:1987, j], 
        raw.dat$catch[raw.dat$species=="Coho" & 
                        raw.dat$catch.year %in% 1965:1988]) # one-year lag
  
  winter.coho.89.13[j] <- 
    cor(win.SST.3[rownames(win.SST.3) %in% 1988:2012, j], 
        raw.dat$catch[raw.dat$species=="Coho" & 
                        raw.dat$catch.year %in% 1989:2013])
  
  winter.coho.14.19[j] <- 
    cor(win.SST.3[rownames(win.SST.3) %in% 2013:2018, j], 
        raw.dat$catch[raw.dat$species=="Coho" & 
                        raw.dat$catch.year %in% 2014:2019])
}

# now pink!
# combining odd and even for our purposes...
winter.pink.65.88 <- winter.pink.89.13 <- winter.pink.14.19 <- NA

pink.combined <- raw.dat %>%
  filter(species %in% c("Pink-odd", "Pink-even")) %>%
  select(Year, species, catch, catch.year) %>%
  arrange(Year)

pink.combined <- na.omit(pink.combined)
ggplot(pink.combined, aes(Year, catch)) +
  geom_line()

for(j in 1:ncol(win.SST.3)){
  # note that we are using catch year as the nominal year matching the proposed eras
  # j <- 1
  winter.pink.65.88[j] <- 
    cor(win.SST.3[rownames(win.SST.3) %in% 1964:1987, j], 
        pink.combined$catch[pink.combined$catch.year %in% 1965:1988]) # 1-yr lag
  
  winter.pink.89.13[j] <- 
    cor(win.SST.3[rownames(win.SST.3) %in% 1988:2012, j], 
        pink.combined$catch[pink.combined$catch.year %in% 1989:2013])
  
  winter.pink.14.19[j] <- 
    cor(win.SST.3[rownames(win.SST.3) %in% 2013:2018, j], 
        pink.combined$catch[pink.combined$catch.year %in% 2014:2019])
}


# png("figs/correlations with 3-yr smoothed winter sst sockeye pink coho.png", 6, 6, units="in", res=300)

tiff("figs/fig1.tiff", 6, 6, units="in", res=300)

par(mfrow=c(3,3), mar=c(0,0.5,3,0.5), oma=c(2.5,1.5,2,1.7), mgp=c(3, 0.2, 0))

new.col <- oceColorsPalette(64)

lim <- c(-1,1)

# first, sockeye
z <- t(matrix(winter.sockeye.65.88,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")

axis(1, labels = F, tcl=0.25, lwd=0.5)
axis(4, labels = F, tcl=0.25, lwd=0.5)

# map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)

mtext("Sockeye", outer=T, cex=1.2, side=2, adj=0.85, line=-0.3)
mtext("1965-1988", outer=T, cex=1.2, side=3, adj=0.1, line=-2.7)

z <- t(matrix(winter.sockeye.89.13,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")

axis(1, labels = F, tcl=0.25, lwd=0.5)
axis(4, labels = F, tcl=0.25, lwd=0.5)

mtext("1989-2013", outer=T, cex=1.2, side=3, adj=0.5, line=-2.7)

z <- t(matrix(winter.sockeye.14.19,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")

axis(1, labels = F, tcl=0.25, lwd=0.5)
axis(4, labels = c("60 N", "50", "40", "30"), at= seq(60,30,-10), las=1, tcl=0.25, lwd=0.5)

mtext("2014-2019", outer=T, cex=1.2, side=3, adj=0.9, line=-2.7)

# now pink
par(mar=c(1.25,0.5,1.25,0.5))

z <- t(matrix(winter.pink.65.88,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")

axis(1, labels = F, tcl=0.25, lwd=0.5)
axis(4, labels = F, tcl=0.25, lwd=0.5)

mtext("Pink", outer=T, cex=1.2, side=2, adj=0.5, line=-0.3)

z <- t(matrix(winter.pink.89.13,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))

axis(1, labels = F, tcl=0.25, lwd=0.5)
axis(4, labels = F, tcl=0.25, lwd=0.5)
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")

z <- t(matrix(winter.pink.14.19,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")

axis(1, labels = F, tcl=0.25, lwd=0.5)
axis(4, labels = c("60", "50", "40", "30"), at= seq(60,30,-10), las=1, tcl=0.25, lwd=0.5)

# coho
par(mar=c(2.5,0.5,0,0.5))
z <- t(matrix(winter.coho.65.88,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")

par(mgp=c(3, 0.05, 0))
axis(1, labels = c("140 E", "180", "140 W"), at= seq(140,220,40), las=1, tcl=0.25, lwd=0.5)

mgp=c(3, 0.2, 0)
axis(4, labels = F, tcl=0.25)

mtext("Coho", outer=T, cex=1.2, side=2, adj=0.17, line=-0.3)

z <- t(matrix(winter.coho.89.13,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))

par(mgp=c(3, 0.05, 0))
axis(1, labels = c("140 E", "180", "140 W"), at= seq(140,220,40), las=1, tcl=0.25, lwd=0.5)

axis(4, labels = F, tcl=0.25)

map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")

# add legend strip
mt.cex <- 1.1
l.mar <- 6
l.cex <- 1
l.wd <- 0.5
l.l <- 1.2
tc.l <- -0.2


z[1,1] <- -1; z[66,26] <- 1
image.plot(z, legend.only=TRUE, horizontal =TRUE,  legend.lab = "r", 
           smallplot = c(0,1,0.03,0.08), 
           legend.cex=1, col=new.col,
           legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0))) 

par(mgp=c(3, 0.2, 0))

z <- t(matrix(winter.coho.14.19,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")

axis(1, labels = c("140 E", "180", "140 W"), at= seq(140,220,40), las=1, tcl=0.25, lwd=0.5)

par(mgp=c(3, 0.2, 0))
axis(4, labels = c("60", "50", "40", "30"), at= seq(60,30,-10), las=1, tcl=0.25, lwd=0.5)


dev.off()

#############################################################################
# now the same thing for 2-yr mean SST to evaluate the effect of smoothing! #
#############################################################################

# first GOA sockeye
winter.sockeye.65.88 <- winter.sockeye.89.13 <- winter.sockeye.14.19 <- NA

for(j in 1:ncol(win.SST.2)){
  # note that we are using catch year as the nominal year matching the proposed eras
  # j <- 1
  winter.sockeye.65.88[j] <- 
    cor(win.SST.2[rownames(win.SST.2) %in% 1963:1986, j], 
        raw.dat$catch[raw.dat$species=="Sockeye" & 
                        raw.dat$catch.year %in% 1965:1988]) # 2-yr lag!
  
  winter.sockeye.89.13[j] <- 
    cor(win.SST.2[rownames(win.SST.2) %in% 1987:2011, j], 
        raw.dat$catch[raw.dat$species=="Sockeye" & 
                        raw.dat$catch.year %in% 1989:2013])
  
  winter.sockeye.14.19[j] <- 
    cor(win.SST.2[rownames(win.SST.2) %in% 2012:2017, j], 
        raw.dat$catch[raw.dat$species=="Sockeye" & 
                        raw.dat$catch.year %in% 2014:2019])
}

# now coho!
winter.coho.65.88 <- winter.coho.89.13 <- winter.coho.14.19 <- NA

for(j in 1:ncol(win.SST.2)){
  # note that we are using catch year as the nominal year matching the proposed eras
  # j <- 1
  winter.coho.65.88[j] <- 
    cor(win.SST.2[rownames(win.SST.2) %in% 1964:1987, j], 
        raw.dat$catch[raw.dat$species=="Coho" & 
                        raw.dat$catch.year %in% 1965:1988]) # one-year lag
  
  winter.coho.89.13[j] <- 
    cor(win.SST.2[rownames(win.SST.2) %in% 1988:2012, j], 
        raw.dat$catch[raw.dat$species=="Coho" & 
                        raw.dat$catch.year %in% 1989:2013])
  
  winter.coho.14.19[j] <- 
    cor(win.SST.2[rownames(win.SST.2) %in% 2013:2018, j], 
        raw.dat$catch[raw.dat$species=="Coho" & 
                        raw.dat$catch.year %in% 2014:2019])
}

# now pink!
# combining odd and even for our purposes...
winter.pink.65.88 <- winter.pink.89.13 <- winter.pink.14.19 <- NA

for(j in 1:ncol(win.SST.2)){
  # note that we are using catch year as the nominal year matching the proposed eras
  # j <- 1
  winter.pink.65.88[j] <- 
    cor(win.SST.2[rownames(win.SST.2) %in% 1964:1987, j], 
        pink.combined$catch[pink.combined$catch.year %in% 1965:1988]) # 1-yr lag
  
  winter.pink.89.13[j] <- 
    cor(win.SST.2[rownames(win.SST.2) %in% 1988:2012, j], 
        pink.combined$catch[pink.combined$catch.year %in% 1989:2013])
  
  winter.pink.14.19[j] <- 
    cor(win.SST.2[rownames(win.SST.2) %in% 2013:2018, j], 
        pink.combined$catch[pink.combined$catch.year %in% 2014:2019])
}


png("figs/correlations with 2-yr smoothed winter sst sockeye pink coho.png", 6, 6, units="in", res=300)

par(mfrow=c(3,3), mar=c(0,0.5,3,0.5), oma=c(2.5,1.5,2,1.7), mgp=c(3, 0.2, 0))

new.col <- oceColorsPalette(64)

lim <- c(-1,1)

# first, sockeye
z <- t(matrix(winter.sockeye.65.88,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")

axis(1, labels = F, tcl=0.25, lwd=0.5)
axis(4, labels = F, tcl=0.25, lwd=0.5)

# map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)

mtext("Sockeye", outer=T, cex=1.2, side=2, adj=0.85, line=-0.3)
mtext("1965-1988", outer=T, cex=1.2, side=3, adj=0.1, line=-2.7)

z <- t(matrix(winter.sockeye.89.13,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")

axis(1, labels = F, tcl=0.25, lwd=0.5)
axis(4, labels = F, tcl=0.25, lwd=0.5)

mtext("1989-2013", outer=T, cex=1.2, side=3, adj=0.5, line=-2.7)

z <- t(matrix(winter.sockeye.14.19,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")

axis(1, labels = F, tcl=0.25, lwd=0.5)
axis(4, labels = c("60 N", "50", "40", "30"), at= seq(60,30,-10), las=1, tcl=0.25, lwd=0.5)

mtext("2014-2019", outer=T, cex=1.2, side=3, adj=0.9, line=-2.7)

# now pink
par(mar=c(1.25,0.5,1.25,0.5))

z <- t(matrix(winter.pink.65.88,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")

axis(1, labels = F, tcl=0.25, lwd=0.5)
axis(4, labels = F, tcl=0.25, lwd=0.5)

mtext("Pink", outer=T, cex=1.2, side=2, adj=0.5, line=-0.3)

z <- t(matrix(winter.pink.89.13,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))

axis(1, labels = F, tcl=0.25, lwd=0.5)
axis(4, labels = F, tcl=0.25, lwd=0.5)
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")

z <- t(matrix(winter.pink.14.19,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")

axis(1, labels = F, tcl=0.25, lwd=0.5)
axis(4, labels = c("60", "50", "40", "30"), at= seq(60,30,-10), las=1, tcl=0.25, lwd=0.5)

# coho
par(mar=c(2.5,0.5,0,0.5))
z <- t(matrix(winter.coho.65.88,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")

par(mgp=c(3, 0.05, 0))
axis(1, labels = c("140 E", "180", "140 W"), at= seq(140,220,40), las=1, tcl=0.25, lwd=0.5)

mgp=c(3, 0.2, 0)
axis(4, labels = F, tcl=0.25)

mtext("Coho", outer=T, cex=1.2, side=2, adj=0.17, line=-0.3)

z <- t(matrix(winter.coho.89.13,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))

par(mgp=c(3, 0.05, 0))
axis(1, labels = c("140 E", "180", "140 W"), at= seq(140,220,40), las=1, tcl=0.25, lwd=0.5)

axis(4, labels = F, tcl=0.25)

map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")

# add legend strip
mt.cex <- 1.1
l.mar <- 6
l.cex <- 1
l.wd <- 0.5
l.l <- 1.2
tc.l <- -0.2


z[1,1] <- -1; z[66,26] <- 1
image.plot(z, legend.only=TRUE, horizontal =TRUE,  legend.lab = "r", 
           smallplot = c(0,1,0.03,0.08), 
           legend.cex=1, col=new.col,
           legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0))) 

par(mgp=c(3, 0.2, 0))

z <- t(matrix(winter.coho.14.19,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")

axis(1, labels = c("140 E", "180", "140 W"), at= seq(140,220,40), las=1, tcl=0.25, lwd=0.5)

par(mgp=c(3, 0.2, 0))
axis(4, labels = c("60", "50", "40", "30"), at= seq(60,30,-10), las=1, tcl=0.25, lwd=0.5)


dev.off()

###############################################################
# and one-year (unsmoothed) PDO #
###############################################################

# first GOA sockeye
winter.sockeye.65.88 <- winter.sockeye.89.13 <- winter.sockeye.14.19 <- NA

for(j in 1:ncol(win.SST.1)){
  # note that we are using catch year as the nominal year matching the proposed eras
  # j <- 1
  winter.sockeye.65.88[j] <- 
    cor(win.SST.1[rownames(win.SST.1) %in% 1963:1986, j], 
        raw.dat$catch[raw.dat$species=="Sockeye" & 
                        raw.dat$catch.year %in% 1965:1988]) # 2-yr lag!
  
  winter.sockeye.89.13[j] <- 
    cor(win.SST.1[rownames(win.SST.1) %in% 1987:2011, j], 
        raw.dat$catch[raw.dat$species=="Sockeye" & 
                        raw.dat$catch.year %in% 1989:2013])
  
  winter.sockeye.14.19[j] <- 
    cor(win.SST.1[rownames(win.SST.1) %in% 2012:2017, j], 
        raw.dat$catch[raw.dat$species=="Sockeye" & 
                        raw.dat$catch.year %in% 2014:2019])
}

# now coho!
winter.coho.65.88 <- winter.coho.89.13 <- winter.coho.14.19 <- NA

for(j in 1:ncol(win.SST.1)){
  # note that we are using catch year as the nominal year matching the proposed eras
  # j <- 1
  winter.coho.65.88[j] <- 
    cor(win.SST.1[rownames(win.SST.1) %in% 1964:1987, j], 
        raw.dat$catch[raw.dat$species=="Coho" & 
                        raw.dat$catch.year %in% 1965:1988]) # one-year lag
  
  winter.coho.89.13[j] <- 
    cor(win.SST.1[rownames(win.SST.1) %in% 1988:2012, j], 
        raw.dat$catch[raw.dat$species=="Coho" & 
                        raw.dat$catch.year %in% 1989:2013])
  
  winter.coho.14.19[j] <- 
    cor(win.SST.1[rownames(win.SST.1) %in% 2013:2018, j], 
        raw.dat$catch[raw.dat$species=="Coho" & 
                        raw.dat$catch.year %in% 2014:2019])
}

# now pink!
# combining odd and even for our purposes...
winter.pink.65.88 <- winter.pink.89.13 <- winter.pink.14.19 <- NA

for(j in 1:ncol(win.SST.1)){
  # note that we are using catch year as the nominal year matching the proposed eras
  # j <- 1
  winter.pink.65.88[j] <- 
    cor(win.SST.1[rownames(win.SST.1) %in% 1964:1987, j], 
        pink.combined$catch[pink.combined$catch.year %in% 1965:1988]) # 1-yr lag
  
  winter.pink.89.13[j] <- 
    cor(win.SST.1[rownames(win.SST.1) %in% 1988:2012, j], 
        pink.combined$catch[pink.combined$catch.year %in% 1989:2013])
  
  winter.pink.14.19[j] <- 
    cor(win.SST.1[rownames(win.SST.1) %in% 2013:2018, j], 
        pink.combined$catch[pink.combined$catch.year %in% 2014:2019])
}


png("figs/correlations with 1-yr smoothed winter sst sockeye pink coho.png", 6, 6, units="in", res=300)

par(mfrow=c(3,3), mar=c(0,0.5,3,0.5), oma=c(2.5,1.5,2,1.7), mgp=c(3, 0.2, 0))

new.col <- oceColorsPalette(64)

lim <- c(-1,1)

# first, sockeye
z <- t(matrix(winter.sockeye.65.88,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")

axis(1, labels = F, tcl=0.25, lwd=0.5)
axis(4, labels = F, tcl=0.25, lwd=0.5)

# map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)

mtext("Sockeye", outer=T, cex=1.2, side=2, adj=0.85, line=-0.3)
mtext("1965-1988", outer=T, cex=1.2, side=3, adj=0.1, line=-2.7)

z <- t(matrix(winter.sockeye.89.13,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")

axis(1, labels = F, tcl=0.25, lwd=0.5)
axis(4, labels = F, tcl=0.25, lwd=0.5)

mtext("1989-2013", outer=T, cex=1.2, side=3, adj=0.5, line=-2.7)

z <- t(matrix(winter.sockeye.14.19,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")

axis(1, labels = F, tcl=0.25, lwd=0.5)
axis(4, labels = c("60 N", "50", "40", "30"), at= seq(60,30,-10), las=1, tcl=0.25, lwd=0.5)

mtext("2014-2019", outer=T, cex=1.2, side=3, adj=0.9, line=-2.7)

# now pink
par(mar=c(1.25,0.5,1.25,0.5))

z <- t(matrix(winter.pink.65.88,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")

axis(1, labels = F, tcl=0.25, lwd=0.5)
axis(4, labels = F, tcl=0.25, lwd=0.5)

mtext("Pink", outer=T, cex=1.2, side=2, adj=0.5, line=-0.3)

z <- t(matrix(winter.pink.89.13,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))

axis(1, labels = F, tcl=0.25, lwd=0.5)
axis(4, labels = F, tcl=0.25, lwd=0.5)
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")

z <- t(matrix(winter.pink.14.19,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")

axis(1, labels = F, tcl=0.25, lwd=0.5)
axis(4, labels = c("60", "50", "40", "30"), at= seq(60,30,-10), las=1, tcl=0.25, lwd=0.5)

# coho
par(mar=c(2.5,0.5,0,0.5))
z <- t(matrix(winter.coho.65.88,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")

par(mgp=c(3, 0.05, 0))
axis(1, labels = c("140 E", "180", "140 W"), at= seq(140,220,40), las=1, tcl=0.25, lwd=0.5)

mgp=c(3, 0.2, 0)
axis(4, labels = F, tcl=0.25)

mtext("Coho", outer=T, cex=1.2, side=2, adj=0.17, line=-0.3)

z <- t(matrix(winter.coho.89.13,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))

par(mgp=c(3, 0.05, 0))
axis(1, labels = c("140 E", "180", "140 W"), at= seq(140,220,40), las=1, tcl=0.25, lwd=0.5)

axis(4, labels = F, tcl=0.25)

map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")

# add legend strip
mt.cex <- 1.1
l.mar <- 6
l.cex <- 1
l.wd <- 0.5
l.l <- 1.2
tc.l <- -0.2


z[1,1] <- -1; z[66,26] <- 1
image.plot(z, legend.only=TRUE, horizontal =TRUE,  legend.lab = "r", 
           smallplot = c(0,1,0.03,0.08), 
           legend.cex=1, col=new.col,
           legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0))) 

par(mgp=c(3, 0.2, 0))

z <- t(matrix(winter.coho.14.19,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")

axis(1, labels = c("140 E", "180", "140 W"), at= seq(140,220,40), las=1, tcl=0.25, lwd=0.5)

par(mgp=c(3, 0.2, 0))
axis(4, labels = c("60", "50", "40", "30"), at= seq(60,30,-10), las=1, tcl=0.25, lwd=0.5)


dev.off()

#############
# everything below is old

# 
# 
# library(oce)
# 
# col <- oceColorsPalette(64)
# 
# z <- t(matrix(winter.coho.14.19,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
# imagep(sst.x,sst.y,z, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
# 
# image(sst.x,sst.y,z, col=col, 
#       zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
# 

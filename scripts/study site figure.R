library(ncdf4)
library(zoo)
library(gplots)
library(dplyr)
library(maps)
library(mapdata)
library(chron)
library(fields)
library(tidyr)
library(nlme)
library(ggplot2)
library(oce)


# load ERSST
nc <- nc_open("/Users/MikeLitzow 1/Documents/R/climate-data/data/North.Pacific.ersst")

# process sst first
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

# now limit to GOA
# extract study area
# 54-62 deg. N, 200-226 deg. E
keep1 <- sst.lon %in% 198:226 
keep2 <- sst.lat %in% 54:62

sst.lon <- sst.lon[keep1]
sst.lat <- sst.lat[keep2]

sst.x <- sst.x[sst.x %in% 198:226]
sst.y <- sst.y[sst.y %in% 54:62]

SST <- SST[,keep1 & keep2]

# need to drop Bristol Bay cells
BB <- c("N58E200", "N58E202", "N56E200", "N56E198", "N58E198", "N60E198")
SST[,BB] <- NA

# plot to check
SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=tim.colors(64))
contour(sst.x, sst.y, z, add=T) 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)

sst.m <- months(sst.d)
sst.yr <- years(sst.d)

# restrict to winter 1965-2019 for plotting
keep <- sst.m %in% c("Nov", "Dec", "Jan", "Feb", "Mar") & sst.yr >= 1965
SST <- SST[keep,]


# total GOA FMA wind stress

# GODAS
# first v stress
nc <- nc_open("/Users/MikeLitzow 1/Documents/R/climate-data/data/North.Pacific.godas.vflxsfc")
# view dates (middle of month):
raw <- ncvar_get(nc, "time")
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))

v.stress <- ncvar_get(nc, "vflxsfc") # get all the data!
x <- ncvar_get(nc, "longitude")     # view longitudes (degrees East)
y <- ncvar_get(nc, "latitude")     # view latitudes
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes

# process
v.stress <- aperm(v.stress, 3:1)  
v.stress <- matrix(v.stress, nrow=dim(v.stress)[1], ncol=prod(dim(v.stress)[2:3])) 

nc <- nc_open("/Users/MikeLitzow 1/Documents/R/climate-data/data/North.Pacific.godas.uflxsfc")
u.stress <- ncvar_get(nc, "uflxsfc")
u.stress <- aperm(u.stress, 3:1)  
u.stress <- matrix(u.stress, nrow=dim(u.stress)[1], ncol=prod(dim(u.stress)[2:3])) 

stress <- sqrt(v.stress^2 + u.stress^2)

dimnames(stress) <-  list(as.character(d), paste("N", lat, "E", lon, sep=""))

# check
z <- t(matrix(colMeans(stress, na.rm=T),length(y)))  
image(x,y,z, col=tim.colors(64), xlim=c(180,250), ylim=c(50,64))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
# block out evertyhing but the area we want
# drop a series of points to mark out the AK Peninsula and Bristol Bay

stress[,lat<55] <- NA
stress[,lon<200] <- NA

# check
z <- t(matrix(colMeans(stress, na.rm=T),length(y)))  
image(x,y,z, col=tim.colors(64), xlim=c(180,250), ylim=c(50,64))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)

# retain FMA only
stress.m <- months(d)

stress <- stress[stress.m %in% c("Feb", "Mar", "Apr"),]

# save coordinate and means for study site Fig.
# save x and y for study site fig
x.stress <- x
y.stress <- y

plot.stress <- colMeans(stress, na.rm = T)
sum(!is.na(plot.stress)) # just to check - 213 positive cells, which looks right

# now ssh pc1
# now the same approach for GODAS
nc <- nc_open("/Users/MikeLitzow 1/Documents/R/climate-data/data/North.Pacific.godas.sshgsfc")

# view dates (middle of month):
raw <- ncvar_get(nc, "time")
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))

ssh <- ncvar_get(nc, "sshgsfc") # get all the data!
x <- ncvar_get(nc, "longitude")     # view longitudes (degrees East)
y <- ncvar_get(nc, "latitude")     # view latitudes
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes

# process!
ssh <- aperm(ssh, 3:1)  # First, reverse order of dimensions ("transpose" array)

ssh <- matrix(ssh, nrow=dim(ssh)[1], ncol=prod(dim(ssh)[2:3]))  # Change to matrix

dimnames(ssh) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# limit to GOA
ssh[,lat < 52] <- NA
ssh[,lon < 195] <- NA

ssh[,lat > 57 & lon < 201] <- NA
ssh[,lon < 195] <- NA

# and dropping the southeast section which seems to dominate EOF1 for this version!
ssh[,lat < 57 & lon > 220] <- NA


yy <- c(54, 58)
xx <- c(195, 204.5)

mod <- lm(yy ~ xx)

newdata <- data.frame(xx=lon)
predlat <- predict(mod, newdata = newdata)

drop <- NA
for(i in 1:ncol(ssh)){
  
  drop[i] <- ifelse(lon[i] > 204.5, FALSE,
                    ifelse(lon[i] < 204.5 & lat[i] > predlat[i], TRUE, FALSE))
}

ssh[,drop] <- NA

# check
SSH.mean <- colMeans(ssh, na.rm=T)
z <- t(matrix(SSH.mean,length(y)))  
image.plot(x,y,z, col=tim.colors(64), xlim=c(180,250), ylim=c(50,64))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)


# and the reduced area
SSH.mean <- colMeans(ssh, na.rm=T)
z <- t(matrix(SSH.mean,length(y)))  
image.plot(x,y,z, col=tim.colors(64), xlim=c(200,220), ylim=c(54,60))
contour(x, y, z, add=T) 
map('world2Hires',fill=F,add=T, lwd=2)
mtext("GODAS")

# now limit to FMA and fit EOF!
m <- months(d)
yr <- years(d)

ssh.godas.fma <- ssh[m %in% c("Feb", "Mar", "Apr"),]
yr.fma <- yr[m %in% c("Feb", "Mar", "Apr")]

# check
SSH.mean <- colMeans(ssh.godas.fma)
z <- t(matrix(SSH.mean,length(y)))  
image(x,y,z, col=tim.colors(64), xlim=c(180,250), ylim=c(50,64))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)

# note that there are a few NAs - but not enough to change the results for EOF!

# use SVD on covariance
# need to drop NAs!
ssh.land <- is.na(colMeans(ssh.godas.fma))
ssh.godas.fma <- ssh.godas.fma[,!land]

eof <- svd(cov(ssh.godas.fma))
plot.ssh <- -eof$u[,1] # reversing to have the same sign as mean differences 

# plot to check
z <- rep(NA, ncol(ssh.godas.fma))
z[!ssh.land] <- plot.ssh 

z <- t(matrix(z,length(y)))  
image(x,y,z, col=tim.colors(64), xlim=c(180,250), ylim=c(50,64))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)


# make a study site figure
tiff("figs/study site figure.tiff", 8,5.5, units="in", res=300)
# png("figs/study site figure.png", 8,5.5, units="in", res=300)
par(mfrow=c(2,2), tcl=-0.2, cex.lab=0.8, cex.axis=0.8, mar=c(1,2,1,0.5))
mt.cex <- 1.1
l.mar <- 3
l.cex <- 1.3
land.col <- "lightyellow3"
xlim <- c(190, 233)
ylim <- c(49,62)
new.col <- oceColorsPalette(64)

z <- t(matrix(plot.stress,length(y.stress)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
zlim <- c(0, max(z, na.rm=T))
image.plot(x.stress,y.stress,z, col=new.col,  yaxt="n", xaxt="n", axis.args = list(mgp=c(3,0.5,0)),
           zlim=zlim, xlim=xlim, ylim=ylim, legend.cex=l.cex, xlab="", ylab="")
map('world2Hires', c('Canada', 'usa'), fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col=land.col)
gak.y <- 59 + (50.7/60)
gak.x <- 360-149 + (28/60)
points(gak.x, gak.y, pch=21, bg="#CC79A7", col="black", cex=1.7)
text(gak.x+4, gak.y-0.2, "GAK1", col="#CC79A7", cex=1.2)

par(mgp=c(3,0,0))
axis(1, labels = c("160 W", "150", "140", "130"), at= seq(200, 230, 10), las=1, tcl=0.25, lwd=0.5)
par(mgp=c(3,0.2,0))
axis(2, labels = c("60 N", "58", "56", "54", "52", "50"), at= seq(60,50,-2), las=1, tcl=0.25, lwd=0.5)

mtext("a) Wind stress (Pa) and GAK1 site", adj=0, cex=0.95)

z <- rep(NA, ncol(ssh.godas.fma))
z[!ssh.land] <- plot.ssh 
z <- t(matrix(z,length(y)))  
image.plot(x,y,z, col=new.col,  yaxt="n", xaxt="n", axis.args = list(mgp=c(3,0.5,0)),
           xlim=xlim, ylim=ylim, legend.cex=l.cex, xlab="", ylab="")
contour(x, y, z, add=T) 
map('world2Hires', c('Canada', 'usa'), fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col=land.col)
points(360-145,50, pch=21, bg="#56B4E9", col="black", cex=1.5)
text(360-145, 50, "Ocean Station Papa ", pos=2, cex=1.2, col="#56B4E9")

par(mgp=c(3,0,0))
axis(1, labels = c("160 W", "150", "140", "130"), at= seq(200, 230, 10), las=1, tcl=0.25, lwd=0.5)
par(mgp=c(3,0.2,0))
axis(2, labels = c("60 N", "58", "56", "54", "52", "50"), at= seq(60,50,-2), las=1, tcl=0.25, lwd=0.5)

mtext("b) Sea surface height (EOF1) and Station Papa", adj=0, cex=0.95)

z <- t(matrix(colMeans(SST), length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image.plot(sst.x,sst.y,z, col=new.col,  yaxt="n", xaxt="n", axis.args = list(mgp=c(3,0.5,0)),
           zlim=c(4,9), xlim=xlim, ylim=ylim, legend.cex=l.cex, xlab="", ylab="")
box()
map('world2Hires', c('Canada', 'usa'), fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col=land.col)

par(mgp=c(3,0,0))
axis(1, labels = c("160 W", "150", "140", "130"), at= seq(200, 230, 10), las=1, tcl=0.25, lwd=0.5)
par(mgp=c(3,0.2,0))
axis(2, labels = c("60 N", "58", "56", "54", "52", "50"), at= seq(60,50,-2), las=1, tcl=0.25, lwd=0.5)

mtext("c) Sea surface temperature (°C)", adj=0, cex=0.95)

# location...
# study site box
x1 <- c(190,190,233,233,190)
y1 <- c(49,62,62,49,49)

# NPI box
# calculate NDJFM NPI - 
# 30º-65ºN, 160ºE-140ºW 
npi.x <- c(160-1.25,160-1.25,221.25,221.25,160-1.25)
npi.y <- c(30-1.25,66.25,66.25,30-1.25,30-1.25)
par(mar=c(1,2,1,4))
plot(10,10, type="n", xlim=c(130,250), ylim=c(20,66), xlab="", ylab="", xaxt="n", yaxt="n")
map('world2Hires', c('Canada', 'usa', 'Mexico', 'hawaii', 'ussr', 'China', 'Japan', 'South Korea', 'North Korea'), fill=T, add=T, lwd=0.5, col=land.col)

# map('world2Hires',fill=F,add=T, lwd=1)
lines(x1,y1, lwd=2, col="#D55E00")
lines(npi.x,npi.y, lwd=2, col="#0072B2")

text(130, 24, "Study site", cex=1.2, col="#D55E00", pos=4)
text(130, 20.5, "North Pacific Index", cex=1.2, col="#0072B2", pos=4)

par(mgp=c(3,0,0))
axis(1, labels = c("140 E", "160", "180", "160 W", "140", "120"), at= seq(140, 240, 20), las=1, tcl=0.25, lwd=0.5)
par(mgp=c(3,0.2,0))
axis(2, labels = c("60 N", "50", "40", "30", "20"), at= seq(60,20,-10), las=1, tcl=0.25, lwd=0.5)

mtext("d) Study site and North Pacific Index", adj=0, cex=0.95)
dev.off()

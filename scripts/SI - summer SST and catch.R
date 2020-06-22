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

theme_set(theme_bw())

raw.dat <- read.csv("data/salmon.and.covariate.data.csv")

# note that years in this file are already lagged to ocean entry!
# re-lag to define catch year!
raw.dat$catch.year <- ifelse(raw.dat$species=="Sockeye", raw.dat$Year+2, raw.dat$Year+1)

# load ERSST
# uncomment these lines to download data
# identify latest year and month needed
year <- 2019
month <- "12"

URL <- paste("https://coastwatch.pfeg.noaa.gov/erddap/griddap/nceiErsstv5.nc?sst[(1854-01-01):1:(", year, "-", month, "-01T00:00:00Z)][(0.0):1:(0.0)][(20):1:(70)][(120):1:(250)]", sep="")

download.file(URL, "data/North.Pacific.ersst")

# open netcdf file of SST data
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

# a-ok!

# now get MJJA means
m <- months(sst.d)
yr <- as.numeric(as.character(years(sst.d)))

SST <- rowMeans(SST, na.rm=T)

# limit to summer
SST <- SST[m %in% c("May", "Jun", "Jul", "Aug")]
summer.yr <- yr[m %in% c("May", "Jun", "Jul", "Aug")] 

# and get annual summer means
summer.sst <- tapply(SST, summer.yr, mean)

plot(names(summer.sst), summer.sst, type="l") # old data clearly not good at this scale!!

plot(names(summer.sst), summer.sst, type="l", xlim=c(1960,2020), ylim=c(8,12))

# and get 3-yr rolling mean 
summer.sst.3 <- rollmean(summer.sst, 3, fill=NA)

plot(names(summer.sst), summer.sst, type="l", xlim=c(1960,2020), ylim=c(8,12))
lines(names(summer.sst), summer.sst.3, col="red")

# and replace the old (winter) SST with summer
new.dat <- raw.dat %>%
  select(Year, species, catch, catch.year)

new.dat <- na.omit(new.dat)

new.dat$SST1 <- summer.sst[match(new.dat$Year, names(summer.sst))]
new.dat$SST3 <- summer.sst.3[match(new.dat$Year, names(summer.sst.3))]

new.dat <-  new.dat %>%
  pivot_longer(cols=c(SST1, SST3))

# and set era
new.dat$era <- ifelse(new.dat$catch.year <= 1988, "1965-1988", "1989-2019")

cb <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
        "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


ggplot(new.dat, aes(value, catch, color=era)) +
  geom_point() +
  geom_smooth(method="gam", se=F) +
  facet_grid(species ~ name)

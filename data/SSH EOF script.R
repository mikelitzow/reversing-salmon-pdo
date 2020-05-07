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
library(pracma)

# coastal ssh
# creating an EOF to capture variability across the GOA in response to reviewer #1 suggestion

# first SODA
nc <- nc_open("/Users/MikeLitzow 1/Documents/R/NSF-GOA/paper 1 final/submitted/hawaii_3e19_7ccd_16ff_e908_64db_63e9.nc")

# view dates (middle of month):
raw <- ncvar_get(nc, "time")
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))

ssh <- ncvar_get(nc, "ssh") # get all the data!
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

# and dropping the southeast section which seems to dominate EOF1 for godas!
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
SSH.mean <- colMeans(ssh)
z <- t(matrix(SSH.mean,length(y)))  
image(x,y,z, col=tim.colors(64), xlim=c(180,250), ylim=c(50,64))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)

# and reduced areas for plots...
SSH.mean <- colMeans(ssh)
z <- t(matrix(SSH.mean,length(y)))  
image(x,y,z, col=tim.colors(64), xlim=c(200,220), ylim=c(54,60))
contour(x, y, z, add=T)  
map('world2Hires',fill=F, xlim=c(200,220), ylim=c(54,60), add=T, lwd=2)
mtext("SODA")
lines(box1.xx, box1.yy)
lines(box2.xx, box2.yy)


# now limit to FMA and fit EOF!
m <- months(d)
yr <- years(d)

ssh.soda.fma <- ssh[m %in% c("Feb", "Mar", "Apr"),]
yr.fma <- yr[m %in% c("Feb", "Mar", "Apr")]

# check
SSH.mean <- colMeans(ssh.soda.fma)
z <- t(matrix(SSH.mean,length(y)))  
image(x,y,z, col=tim.colors(64), xlim=c(180,250), ylim=c(50,64))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)

# now get box1 and box2 time series for soda

keep <- inpolygon(lon, lat, box1.xx, box1.yy)
box1.soda <- ssh.soda.fma
box1.soda[,!keep] <- NA

keep <- inpolygon(lon, lat, box2.xx, box2.yy)
box2.soda <- ssh.soda.fma
box2.soda[,!keep] <- NA

# check
SSH.mean <- colMeans(box2.soda)
z <- t(matrix(SSH.mean,length(y)))  
image(x,y,z, col=tim.colors(64), xlim=c(200,220), ylim=c(54,60))
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
# correcto!

# and get area means, then average by year
box1.soda <- rowMeans(box1.soda, na.rm = T)
box2.soda <- rowMeans(box2.soda, na.rm = T)

box1.soda <- tapply(box1.soda, yr.fma, mean)
box2.soda <- tapply(box2.soda, yr.fma, mean)
box.diff.soda <- box1.soda - box2.soda

plot <- data.frame(year=1949:2010,
                   box1=box1.soda,
                   box2=box2.soda)

ggplot(plot, aes(box1, box2)) +
  theme_bw() +
  geom_point()

plot <- plot %>%
  pivot_longer(cols = -year)

ggplot(plot, aes(year, value, color=name)) +
  theme_bw() +
  geom_line()

plot <- data.frame(year=1949:2010,
                   diff=box1.soda-box2.soda)

ggplot(plot, aes(year, diff)) +
  theme_bw() +
  geom_line()

# use SVD on covariance
# need to drop NAs!
land <- is.na(colMeans(ssh.soda.fma))
ssh.soda.fma <- ssh.soda.fma[,!land]

eof <- svd(cov(ssh.soda.fma))

# plot to check
z <- rep(NA, ncol(ssh))
z[!land] <- eof$u[,1]

z <- t(matrix(z,length(y)))  
image(x,y,z, col=tim.colors(64), xlim=c(180,250), ylim=c(50,64))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)

# and the time series of PC1
soda.ssh.pc1 <- ssh.soda.fma %*% eof$u[,1]

# and annual FMA means for this value
soda.ssh.pc1 <- tapply(soda.ssh.pc1, yr.fma, mean)

# plot to check
plot(1949:2010, soda.ssh.pc1, type="l") # looks right re. 1976/77

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

box1.xx <- c(207.5, 210.5, 210.5, 207.5, 207.5)
box1.yy <- c(59.5, 59.5, 58.5, 58.5, 59.5)


box2.xx <- c(212, 215, 215, 212, 212)
box2.yy <- c(56.5, 56.5, 57.5, 57.5, 56.5)
# and the reuced area
SSH.mean <- colMeans(ssh, na.rm=T)
z <- t(matrix(SSH.mean,length(y)))  
image.plot(x,y,z, col=tim.colors(64), xlim=c(200,220), ylim=c(54,60))
contour(x, y, z, add=T) 
map('world2Hires',fill=F,add=T, lwd=2)
mtext("GODAS")
lines(box1.xx, box1.yy)
lines(box2.xx, box2.yy)

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
land <- is.na(colMeans(ssh.godas.fma))
ssh.godas.fma <- ssh.godas.fma[,!land]

eof <- svd(cov(ssh.godas.fma))

# plot to check
z <- rep(NA, ncol(ssh))
z[!land] <- eof$u[,1]

z <- t(matrix(z,length(y)))  
image(x,y,z, col=tim.colors(64), xlim=c(180,250), ylim=c(50,64))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)

# and the time series of PC1
godas.ssh.pc1 <- ssh.godas.fma %*% eof$u[,1]

# and annual FMA means for this value
godas.ssh.pc1 <- tapply(godas.ssh.pc1, yr.fma, mean)

# combine and plot to check
soda <- data.frame(year=1949:2010, name="soda", value=soda.ssh.pc1)
godas <- data.frame(year=1980:2019, name="godas", value=-godas.ssh.pc1)

plot.dat <- full_join(soda, godas)

ggplot(plot.dat, aes(year, value, color=name)) +
  theme_bw() +
  geom_line()
  
plot.dat <- plot.dat %>%
  filter(year %in% 1980:2010) %>%
  pivot_wider(names_from = name, values_from = value)

ggplot(plot.dat, aes(soda, godas)) +
  theme_bw() +
  geom_point()

cor(plot.dat$soda, plot.dat$godas) # 0.73 - worse - 0.69!!

# might be worth trying the stan models both ways to make sure that we're not getting different results for different 
# ways of combining!

mod1 <- lm(plot.dat$soda ~ plot.dat$godas)
summary(mod1)
mod1$coefficients

soda.estimate <- coef(mod1)[2]*godas$value + coef(mod1)[1]
names(soda.estimate) <- 1980:2019
  
plot(1980:2019, soda.estimate, type="l", ylim=c(-1.5,2.5))
lines(1949:2010, soda$value, type="l", col="red")

soda.predicted.by.godas <- c(soda$value[soda$year %in% 1965:1979], 
                         soda.estimate)

plot(1965:2019, soda.predicted.by.godas, type="l")

# estimate pre-1980 GODAS values from SODA
mod2 <- lm(plot.dat$godas ~ plot.dat$soda)
summary(mod2)
mod2$coefficients

godas.estimate <- coef(mod2)[2]*soda$value + coef(mod2)[1]

plot(1949:2010, godas.estimate, type="l", ylim=c(-4, -1))
lines(1980:2019, godas$value, type="l", col="red")

names(godas.estimate) <- 1949:2010
godas.predicted.by.soda <- c(godas.estimate[names(godas.estimate) %in% 1965:1979], godas$value)

plot(1965:2019, scale(godas.predicted.by.soda), type="l")
lines(1965:2019, scale(soda.predicted.by.godas), type="l", col="red")

# that's a reassuring degree of correspondence!

# and add to the covariates object
salmon.covariates <- read.csv("data/salmon.covariates.csv", row.names = 1)

salmon.covariates$FMA.ssh.PC1.GODAS.pred.by.SODA <-
  godas.predicted.by.soda[match(salmon.covariates$year, names(godas.predicted.by.soda))]

salmon.covariates$FMA.ssh.PC1.SODA.pred.by.GODAS <-
  soda.predicted.by.godas[match(salmon.covariates$year, names(soda.predicted.by.godas))]

# and save! 
write.csv(salmon.covariates, "data/salmon.covariates.updated.ssh.eof.csv", row.names = F)

# everything below is old...
# salmon.covariates$FMA.ssh.COMBINED <- 
#   godas.synthetic.fma[match(salmon.covariates$year, names(godas.synthetic.fma))] 
# 
# # 
# 
# # save coordinate and means for study site Fig.
# # save x and y for study site fig
# x.ssh <- x
# y.ssh <- y
# 
# m <- months(d)
# yr <- years(d)
# 
# ssh.godas.fma <- ssh.mean.godas[m %in% c("Feb", "Mar", "Apr")]
# yr.fma <- yr[m %in% c("Feb", "Mar", "Apr")]
# ssh.godas.fma <- tapply(ssh.godas.fma, yr.fma, mean)
# 
# plot(ssh.godas.fma[names(ssh.godas.fma <= 2010)], ssh.soda.fma[names(ssh.godas.fma >=1980)])
# cor(ssh.godas.fma[names(ssh.godas.fma <= 2010)], ssh.soda.fma[names(ssh.godas.fma >=1980)], use="p") # 0.778
# 
# # estimate pre-1980 GODAS values from SODA
# mod <- lm(ssh.godas.fma[names(ssh.godas.fma <= 2010)] ~ ssh.soda.fma[names(ssh.godas.fma >=1980)], 
#           na.action="na.exclude")
# godas.estimate.fma <- coef(mod)[2]*ssh.soda.fma + coef(mod)[1]
# 
# plot(1949:2010, godas.estimate.fma, type="l", col="red", xlim=c(1949,2019), ylim=c(0.02,0.22))
# lines(1980:2019, ssh.godas.fma, col="blue")
# 
# godas.synthetic.fma <- c(godas.estimate.fma[names(godas.estimate.fma) <1980], ssh.godas.fma)
# 
# plot(1949:2019, godas.synthetic.fma, type="l")
# 
# salmon.covariates$FMA.ssh.GODAS <-
#   ssh.godas.fma[match(salmon.covariates$year, names(ssh.godas.fma))]
# 
# salmon.covariates$FMA.ssh.SODA <- 
#   ssh.soda.fma[match(salmon.covariates$year, names(ssh.soda.fma))] 
# 
# salmon.covariates$FMA.ssh.COMBINED <- 
#   godas.synthetic.fma[match(salmon.covariates$year, names(godas.synthetic.fma))] 
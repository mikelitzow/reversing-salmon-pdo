# this script examines catch-PDO relationships at different lags / degrees of PDO smoothing

library(zoo)
library(tidyverse)

# set theme for ggplot
theme_set(theme_bw())


# download from SWFSC ERDDAP
# and separate time field into year / month

download.file("https://oceanview.pfeg.noaa.gov/erddap/tabledap/cciea_OC_PDO.csv?time%2CPDO&time%3E=1950-01-01&time%3C=2020-05-01T00%3A00%3A00Z", "data/latest.pdo.csv")

pdo <- read.csv("data/latest.pdo.csv", skip=1, col.names = c("time", "PDO"))

pdo <- pdo %>%
  separate(time, c("year", "month","junk"), "-") %>%
  select(-junk)
 
pdo$month <- as.numeric(pdo$month) 
pdo$year <- as.numeric(pdo$year) 


# and get a winter (NDJFM) value!
pdo$winter.year = ifelse(pdo$month %in% 11:12, pdo$year+1, pdo$year)

winter.pdo <- pdo %>%
  filter(month %in% c(11,12,1:3)) %>%
  group_by(winter.year) %>%
  summarise(winter.pdo=mean(PDO))

names(winter.pdo) <- c("Year", "PDO1")

# get smoothed values of PDO for analysis
# aligning left - i.e, year of and year after ocean entry
winter.pdo$PDO2 <- rollmean(winter.pdo$PDO1, 2, align = "left", fill=NA) 
winter.pdo$PDO3 <- rollmean(winter.pdo$PDO1, 3, align = "center", fill=NA) 

# load salmon catch data
dat <- read.csv("data/salmon.and.covariate.data.csv")

# lag back to catch year to make the era definitions match up across spp.!
dat$catch.year <- ifelse(dat$species=="Sockeye", dat$Year+2, dat$Year+1)

dat <- left_join(dat, winter.pdo)

# now get cross corrs for each species and smoothing for data in each era!
dat$era <- ifelse(dat$catch.year <= 1988, "1965-1988",
                  ifelse(dat$catch.year %in% 1989:2013, "1989-2013", NA))


# lump pinks - short time series
dat$species.grouped <- ifelse(dat$species %in% c("Pink-even", "Pink-odd"), "Pink",
                              as.character(dat$species))
spp <- as.factor(unique(dat$species.grouped))
eras <- as.factor(na.omit(unique(dat$era)))

dat$era <- as.factor(dat$era)

cross.corr <- data.frame()

for(s in spp){ # select spp.
  for(e in eras){
    
    # s <- spp[1]
    # e <- eras[1]
    
    temp <- dat %>%
      filter(species.grouped==s, era==e) %>%
      na.omit()
    
    cross.corr <- rbind(cross.corr,
                        data.frame(spp=s,
                                   era=e,
                                   lag=-4:4,
                                   PDO1=ccf(temp$PDO1, temp$catch, lag.max=4)$acf,
                                   PDO2=ccf(temp$PDO2, temp$catch, lag.max=4)$acf,
                                   PDO3=ccf(temp$PDO3, temp$catch, lag.max=4)$acf))
    
  }
}

cross.corr$lag <- as.factor(cross.corr$lag)

cross.corr <- cross.corr %>%
  pivot_longer(cols = c(-spp, -era, -lag))

# set colors
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# and reorder spp.
cross.corr$plot.order <- ifelse(cross.corr$spp=="Pink", 1, 
                                ifelse(cross.corr$spp=="Sockeye", 2, 3))

cross.corr$spp <- reorder(cross.corr$spp, cross.corr$plot.order)

ggplot(cross.corr, aes(lag, value, fill=name)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(values=cb[2:4], 
                    labels=c("One-year mean", 
                             "Two-year rolling mean", 
                             "Three-year rolling mean")) +
  facet_grid(spp~era) +
  geom_hline(yintercept = 0, col="dark grey") +
  ylab("Pearson correlation") +
  xlab("Lag (years)") +
  theme(legend.position = "top", legend.direction="horizontal", legend.title = element_blank())

ggsave("figs/PDO-catch correlations at different lags.png", width=5, height = 6, units='in')

#############
# now, save as a new data file for use in the stan models!
dat <- dat %>%
  select(-catch.year, -era, -species.grouped)
  
    
write.csv(dat, "data/salmon.and.SWFSC.PDO.csv", row.names = F)

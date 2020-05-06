library(zoo)
library(tidyverse)

# set theme for ggplot
theme_set(theme_bw())


# PDO isn't currently being updated on JISAO website - 
# download from NCDC

download.file("https://www.ncdc.noaa.gov/teleconnections/pdo/data.csv", "data/latest.pdo.csv")

temp <- read.csv("data/latest.pdo.csv")

year <- floor(as.numeric(rownames(temp))/100)
month <- as.numeric(rownames(temp)) - year*100
pdo <- as.numeric(as.character(temp[,1]))

pdo <- data.frame(year=year[2:length(year)],
                  month=month[2:length(month)],
                  pdo=pdo[2:length(pdo)])

# and get a winter (NDJFM) value!
pdo$winter.year = ifelse(pdo$month %in% 11:12, pdo$year+1, pdo$year)

winter.pdo <- pdo %>%
  filter(month %in% c(11,12,1:3)) %>%
  group_by(winter.year) %>%
  summarise(winter.pdo=mean(pdo))
names(winter.pdo) <- c("Year", "NCDC.PDO1")

# load salmon catch data
dat <- read.csv("data/salmon.and.covariate.data.csv")


# add NCDC version of PDO!
dat <- left_join(dat, winter.pdo)

# compare the new NCDC version of PDO with the version we were using before...
check.dat <- dat %>%
  select(Year, PDO1, NCDC.PDO1) %>%
  pivot_longer(cols=-Year)

ggplot(check.dat, aes(Year, value, color=name)) +
  geom_line()

# generally close - not such extreme positive values in NCDC version in recent years

ggplot(dat, aes(PDO1, NCDC.PDO1)) +
  geom_point()

cor(dat$PDO1, dat$NCDC.PDO1) # 0.932...so better than my estimates!

# these catch data are pre-lagged to ocean entry...
# get smoothed values of NCDC PDO for analysis

# aligning left - i.e, year of and year after ocean entry
winter.pdo$NCDC.PDO2 <- rollmean(winter.pdo$NCDC.PDO1, 2, align = "left", fill=NA) 
winter.pdo$NCDC.PDO3 <- rollmean(winter.pdo$NCDC.PDO1, 3, align = "center", fill=NA) 

dat <- left_join(dat, winter.pdo)

# now get cross corrs for each species and smoothing for data in each era!
dat$era <- ifelse(dat$Year <= 1988, "1965-1988",
                  ifelse(dat$Year %in% 1989:2013, "1989-2013", NA))


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
                                   PDO1=ccf(temp$NCDC.PDO1, temp$catch, lag.max=4)$acf,
                                   PDO2=ccf(temp$NCDC.PDO2, temp$catch, lag.max=4)$acf,
                                   PDO3=ccf(temp$NCDC.PDO3, temp$catch, lag.max=4)$acf))
    
  }
}

cross.corr$lag <- as.factor(cross.corr$lag)

cross.corr <- cross.corr %>%
  pivot_longer(cols = c(-spp, -era, -lag))

# set colors
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(cross.corr, aes(lag, value, fill=name)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(values=cb[2:4], 
                    labels=c("Winter PDO: year of ocean entry", 
                             "Winter PDO: year of and year after ocean entry", 
                             "Winter PDO: year before, year of, and year after")) +
  facet_grid(spp~era) +
  geom_hline(yintercept = 0, col="dark grey") +
  ylab("Pearson correlation") +
  xlab("Lag (years)") +
  theme(legend.position = "top", legend.direction="vertical", legend.title = element_blank())

ggsave("figs/PDO-catch correlations at different lags.png", width=5, height = 7, units='in')

#############
# now, save as a new data file for use in the stan models!
dat <- dat %>%
  select(Year, NCDC.PDO1, NCDC.PDO2, NCDC.PDO3, species, catch)

names(dat)[2:4] <- c("PDO1", "PDO2", "PDO3")      

write.csv(dat, "data/salmon.and.NCDC.PDO.csv", row.names = F)

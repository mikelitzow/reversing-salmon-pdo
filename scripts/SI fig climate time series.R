library(tidyverse)

dat <- read.csv("data/climate.data.csv")

# add SWFSC PDO values
dat2 <- read.csv("data/salmon.and.SWFSC.PDO.csv")

dat2 <- dat2 %>%
  group_by(Year) %>%
  summarise(NDJFM.PDO=mean(PDO1))

names(dat2)[1] <- "year"

dat <- left_join(dat, dat2)

# and replace the ssh coastal mean with ssh PC1
dat3 <- read.csv("data/salmon.covariates.updated.ssh.eof.csv")
dat$FMA.SSH.PC1 <- dat3$FMA.ssh.PC1.SODA.pred.by.GODAS[match(dat$year, dat3$year)]

dat <- dat %>%
  select(-era, -FMA.ssh.COMBINED) %>%
  filter(year %in% 1965:2019)

names(dat) <- c("Year",
                "Wind stress (Pa)",
                "Pacific Decadal Oscillation",
                "Sea surface temp. (ºC)",
                "North Pacific Index (mb)",
                "Papa advection index (ºN)",
                "GAK1 salinty (PSU)",
                "Sea surface height (PC1)")


dat <- dat %>%
  pivot_longer(cols=-Year)

dat$order <- NA

this <- grep("Pacific Decadal", dat$name)
dat$order[this] <- 1

this <- grep("North", dat$name)
dat$order[this] <- 2

this <- grep("(ºC)", dat$name)
dat$order[this] <- 3

this <- grep("(PC1)", dat$name)
dat$order[this] <- 4

this <- grep("Papa", dat$name)
dat$order[this] <- 5

this <- grep("Wind", dat$name)
dat$order[this] <- 6

this <- grep("GAK1", dat$name)
dat$order[this] <- 7

dat$name <- reorder(dat$name, dat$order)

ggplot(dat, aes(Year, value)) +
  theme_bw() +
  geom_line() +
  facet_wrap(~name, scales="free_y", ncol=2) +
  theme(axis.title = element_blank())

ggsave("figs/SI climate time series.png", width=5, height=6, units='in', dpi=300)

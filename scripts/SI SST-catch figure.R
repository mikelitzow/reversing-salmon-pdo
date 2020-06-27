library(tidyverse)
theme_set(theme_bw())

# NOTE THAT THESE YEARS ARE ALREADY LAGGED TO ENTRY YEAR
# raw.dat <- read.csv("salmon.and.covariate.data.csv")
raw.dat <- read.csv("data/salmon.and.SWFSC.PDO.csv")

# re-lagging to catch year to make this consistent across spp. -
# we will use catch year to distinguish among eras...
raw.dat$catch.year <- ifelse(raw.dat$species=="Sockeye", raw.dat$Year+2, raw.dat$Year+1)

raw.dat[["era"]] <- ifelse(raw.dat$catch.year <= 1988, "era1",
                           ifelse(raw.dat$catch.year %in% 1989:2013, "era2", "era3"))

cb <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
        "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# make a first panel of SST3 and catch for each spp. group, plotted against entry year

temp <- raw.dat %>%
  select(Year, SST3, species, catch, era) %>%
  pivot_wider(names_from = species, values_from = catch) %>%
  mutate(SST3=scale(SST3)) %>%
  pivot_longer(cols=c(-Year, -era))


temp <- na.omit(temp)

change <- temp$name=="SST3"
temp$name[change] <- "SST 3 yr."

temp$order <- ifelse(temp$name=="SST 3 yr.", 1,
                     ifelse(temp$name=="Pink-odd", 2,
                            ifelse(temp$name=="Pink-even", 3,
                                   ifelse(temp$name=="Sockeye", 4, 5))))
temp$name <- reorder(temp$name, temp$order)

panel.1<- ggplot(temp, aes(Year, value, color=name)) +
  geom_hline(yintercept = 0) +
  geom_line() +
  theme(axis.title.x = element_blank(), legend.title = element_blank()) +
  scale_color_manual(values=cb[c(1,2,7,6,4)]) +
  ylab("Anomaly") +
  scale_x_continuous(breaks=seq(1970,2020,10))

panel.1

ggsave("figs/SI - winter sst and catch time series.png", width=6, height=4, units="in")

temp2 <- temp %>%
  select(-order) %>%
  pivot_wider(names_from = name, values_from = value) 

temp2$mean.catch = apply(temp2[,4:7], 1, mean, na.rm = T)

temp2 <- temp2 %>%
  select(era, `SST 3 yr.`, mean.catch)

names(temp2)[2] <- "SST"


temp2$catch.year <- ifelse(temp2$era=="era1", "1965-1988",
                           ifelse(temp2$era=="era2", "1989-2013", "2014-2019"))

temp2$catch.year <- ifelse(temp2$era=="era1", "1965-1988", "1989-2019")

panel.2 <- ggplot(temp2, aes(SST, mean.catch, color=catch.year)) +
  geom_point() +
  geom_smooth(method="gam", se=F) +
  ylab("Mean catch anomaly") +
  xlab("SST anomaly (3-yr. mean)") +
  labs(color = "Catch year") +
  scale_color_manual(values=cb[c(2,6)])

png("figs/SI fig - two panel SST-catch.png", 6, 6, units='in', res=300)
ggpubr::ggarrange(panel.1, panel.2, ncol=1, nrow=2, labels=c("a)", "b)"))
dev.off()  

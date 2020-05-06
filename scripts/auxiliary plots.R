# make some plots for talk(s)!
library(tidyverse)

# NOTE THAT THESE YEARS ARE ALREADY LAGGED TO ENTRY YEAR
# raw.dat <- read.csv("data/salmon.and.covariate.data.csv")
raw.dat <- read.csv("data/salmon.and.NCDC.PDO.csv")
# update - using NCDC PDO values!

raw.dat[["era"]] <- ifelse(raw.dat$Year <= 1986, "era1",
                           ifelse(raw.dat$Year %in% 1987:2011, "era2", "era3"))

cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
        "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# make a plot for AMSS talk

plot.dat <- raw.dat %>%
  select(Year, catch)  %>%
  filter(Year <= 2013) %>%
  group_by(Year) %>%
  summarise(mean.catch=mean(catch, na.rm=T))

plot.dat$PDO <- raw.dat$PDO3[match(plot.dat$Year, raw.dat$Year)]

plot.dat$era <- ifelse(plot.dat$Year < 1988, "1963-1988", "1989-2013")

ggplot(plot.dat, aes(PDO, mean.catch, color=era)) +
  theme_bw() +
  geom_point() +
  geom_smooth(method="lm", se=F) +
  theme(legend.title = element_blank(), legend.position = c(0.78, 0.2),
        axis.text.y = element_blank()) +
  scale_color_manual(values=cb[c(2,3)]) +
  xlab("PDO") +
  ylab("Catch")

ggsave("AMSS two-era catch vs PDO.png", width=3, height=2.5, units='in')

# now an SST version!
plot.dat$SST <- (9/5)*raw.dat$SST3[match(plot.dat$Year, raw.dat$Year)]+32


ggplot(filter(plot.dat, Year %in% 1963:1988), aes(SST, mean.catch)) +
  theme_bw() +
  geom_point(color=cb[2]) +
  geom_smooth(method="lm", se=F, color=cb[2]) +
  theme(legend.title = element_blank(), legend.position = 'none') +
  xlab("Winter temperature (ºF)") +
  ylab("Catch anomaly")

ggsave("one-era catch vs SST.png", width=3.5, height=2.5, units='in')


ggplot(plot.dat, aes(SST, mean.catch, color=era)) +
  theme_bw() +
  geom_point() +
  geom_smooth(method="lm", se=F) +
  theme(legend.title = element_blank(), legend.position = c(0.78, 0.2)) +
  scale_color_manual(values=cb[c(2,3)]) +
  xlab("Winter temperature (ºF)") +
  ylab("Catch anomaly")

ggsave("two-era catch vs SST.png", width=3.5, height=2.5, units='in')

# and dummy tw0-distribution pdf for slope!
dummy.density <- data.frame(era=rep(c("1963-1988", "1989-2013"), each=5000),
                            y=c(rnorm(5000, 0.8, 0.1), rnorm(5000, 0.05, 0.15)))

ggplot(dummy.density, aes(y, fill=era)) +
  theme_bw() +
  geom_density(alpha=0.6, color="dark grey") +
  scale_fill_manual(values=cb[c(2,3)]) +
  xlim(-0.7,2) +
  geom_vline(xintercept = 0, lty=2) +
  ylab("Probability density") +
  xlab("Slope") +
  theme(legend.position = c(0.8, 0.8), legend.title = element_blank())


ggsave("slope dummy two era plot.png", width=3, height=2.5, units="in")


names(plot.dat)[2:3] <- c("GOA salmon catch", "Winter PDO (3-year running mean)")

plot.dat <- plot.dat %>%
  gather(key, value, -era, -Year)

plot.dat$color <- as.factor(ifelse(plot.dat$value < 0, 1, 2))
plot.dat$order <- ifelse(plot.dat$key=="Winter PDO (3-year running mean)", 1, 2)
plot.dat$key <- reorder(plot.dat$key, plot.dat$order)

ggplot(plot.dat, aes(Year, value, fill=color)) +
  theme_bw() +
  geom_bar(position="dodge", stat="identity") +
  facet_wrap(~key, scales="free", nrow=2) +
  theme(legend.position = 'none', axis.title.x = element_blank()) +
  scale_fill_manual(values=c("blue", "red")) +
  ylab("Anomaly") +
  geom_hline(yintercept = 0, color="dark grey") +
  geom_vline(xintercept = 1988.5, lty=2)

ggsave("pdo and goa catch barplot.png", width=3.5, height=4, units="in")

# now a dummy example plot
# load pdo
download.file("http://jisao.washington.edu/pdo/PDO.latest", "~pdo")
names <- read.table("~pdo", skip=30, nrows=1, as.is = T)
pdo <- read.table("~pdo", skip=31, nrows=119, fill=T, col.names = names)
pdo$YEAR <- 1900:(1899+nrow(pdo)) # drop asterisks!
pdo <- pdo %>%
  gather(month, value, -YEAR) %>%
  arrange(YEAR)

pdo <- pdo %>%
  filter(YEAR %in% 1930:1995)

pdo3 <- zoo::rollmean(tapply(pdo$value, pdo$YEAR, mean), 3, fill=NA)
error <- rnorm(length(pdo3), 0, 1)
y <- pdo3 + error

dummy.plot <- data.frame(year=1930:1995,
                         pdo=pdo3,
                         y <- scale(y))

names(dummy.plot)[2:3] <- c("PDO (3-year running mean)", "Generic biological response")

dummy.plot <- dummy.plot %>%
  gather(key, value, -year)

dummy.plot$color <- as.factor(ifelse(dummy.plot$value < 0, 1, 2))
dummy.plot$order <- ifelse(dummy.plot$key=="PDO (3-year running mean)", 1, 2)
dummy.plot$key <- reorder(dummy.plot$key, dummy.plot$order)

ggplot(dummy.plot, aes(year, value, fill=color)) +
  theme_bw() +
  geom_bar(position="dodge", stat="identity") +
  facet_wrap(~key, scales="free", nrow=2) +
  theme(legend.position = 'none', axis.title.x = element_blank()) +
  scale_fill_manual(values=c("blue", "red")) +
  ylab("Anomaly") +
  geom_hline(yintercept = 0, color="dark grey")

ggsave("pdo and biology dummy barplot.png", width=3.5, height=4, units="in")

dummy.plot2 <- data.frame(pdo=pdo3, y=y)

ggplot(dummy.plot2, aes(pdo, y)) +
  theme_bw() +
  geom_point(color=cb[2]) +
  geom_smooth(method="lm", se=F, color=cb[2]) +
  ylab("Biological response") +
  xlab("PDO")

ggsave("pdo and biology dummy scatter plot.png", width=3, height=2.5, units="in")


# and dummy pdf for slope!
dummy.density <- data.frame(y=rnorm(5000, 1, 0.2))

ggplot(dummy.density, aes(y)) +
  theme_bw() +
  geom_density(fill=cb[2], alpha=0.6, color="dark grey") +
  xlim(-0.5,1.8) +
  geom_vline(xintercept = 0, lty=2) +
  ylab("Probability density") +
  xlab("Slope")

ggsave("slope dummy plot.png", width=3, height=2.5, units="in")
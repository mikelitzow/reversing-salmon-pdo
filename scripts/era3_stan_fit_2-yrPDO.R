## Bayesian Catch:PDO reversal analysis

# this script examines catch-PDO relationships with PDO data smoothed with 2-yr running means
# and relationships with unsmoothed PDO data

library(plyr)
library(tidyverse)
library(rstan)
library(ggplot2)
library(rstanarm)
library(bayesplot)
library(overlapping)

# NOTE THAT THESE YEARS ARE ALREADY LAGGED TO ENTRY YEAR
# raw.dat <- read.csv("salmon.and.covariate.data.csv")
raw.dat <- read.csv("data/salmon.and.SWFSC.PDO.csv")

# re-lagging to catch year to make this consistent across spp. -
# we will use catch year to distinguish among eras...
raw.dat$catch.year <- ifelse(raw.dat$species=="Sockeye", raw.dat$Year+2, raw.dat$Year+1)

raw.dat[["era"]] <- ifelse(raw.dat$catch.year <= 1988, "era1",
  ifelse(raw.dat$catch.year %in% 1989:2013, "era2", "era3"))
                           

cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
        "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


## Prep data -----------------------------------------------
dat3 <- na.omit(raw.dat)
dat3 <- plyr::ddply(dat3, .(species), transform, pdo = scale(PDO2))

# dat3 <- plyr::ddply(dat3, .(species), transform, sst = (9/5)*SST3+32) # changing to raw ºF!

dat3$era <- as.factor(dat3$era)

## Find start and end indices for each species
n <- as.numeric(dat3$species)
int.end <- which(diff(n) != 0)
int.start <- int.end + 1
end <- c(int.end, length(n))
start <- c(1, int.start)

## Data for Stan models
dat3_stan <- list(y = dat3$catch,
                  x1 = dat3$pdo,
                  y_start = start,
                  y_end = end,
                  n_species = length(unique(dat3$species)),
                  N = nrow(dat3),
                  era1 = ifelse(dat3$catch.year <= 1988, 1, 0),
                  era2 = ifelse(dat3$catch.year %in% 1989:2013, 1, 0),
                  era3 = ifelse(dat3$catch.year >= 2014, 1, 0))

## 3 era: hierarchical --> model for manuscript ----------------
era3_hier_arm <- stan_glmer(catch ~ era + pdo + pdo:era + (era + pdo + pdo:era | species),
                            data = dat3,
                            chains = 4, cores = 4, thin = 1,
                            warmup = 1000, iter = 4000, refresh = 0,
                            adapt_delta = 0.99,
                            prior = normal(location = 0, scale = 5, autoscale = FALSE),
                            prior_intercept = normal(location = 0, scale = 5, autoscale = FALSE),
                            prior_aux = student_t(df = 3, location = 0, scale = 5, autoscale = FALSE),
                            prior_covariance = decov(regularization = 1,
                                                     concentration = 1,
                                                     shape = 1, scale = 1))

# do predictions for recent era
newdata = expand.grid("era"=unique(dat3$era),
                      pdo=c(mean(dat3$pdo[which(dat3$era=="era3")]),
                            mean(dat3$pdo[which(dat3$era=="era3")])+1), species=unique(dat3$species)) %>%
  dplyr::filter(era=="era3")
pred = posterior_predict(era3_hier_arm, newdata=newdata)
newdata$mean_pred = round(apply(pred,2,mean), 3)

fixef(era3_hier_arm)
ranef(era3_hier_arm)
coef(era3_hier_arm)
era3_hier_arm$covmat
print(era3_hier_arm)

mu_beta  <- as.matrix(era3_hier_arm, pars = c("pdo", "eraera2:pdo", "eraera3:pdo"))
coef_beta <- data.frame(coef = "Slope",
                        era1 = mu_beta[ , 1],
                        era2 = mu_beta[ , 1] + mu_beta[ , 2],
                        era3 = mu_beta[ , 1] + mu_beta[ , 3])
mu_alpha <- as.matrix(era3_hier_arm, pars = c("(Intercept)", "eraera2", "eraera3"))
coef_alpha <- data.frame(coef = "Intercept",
                         era1 = mu_alpha[ , 1],
                         era2 = mu_alpha[ , 1] + mu_alpha[ , 2],
                         era3 = mu_alpha[ , 1] + mu_alpha[ , 3])
mbeta  <- reshape2::melt(coef_beta, id.vars = "coef")
malpha <- reshape2::melt(coef_alpha, id.vars = "coef")
mdf_hier <- rbind(mbeta, malpha)


era_mean_slopes <- mbeta %>%
  group_by(variable) %>%
  summarize(mean=mean(value))

era_mean_slopes[3,2] - era_mean_slopes[2,2]
era_mean_slopes[3,2] - era_mean_slopes[1,2]

# and % positive / negative slopes by era!
sum(mbeta$value[mbeta$variable=="era1"]>0)/ length(mbeta$value[mbeta$variable=="era1"])
sum(mbeta$value[mbeta$variable=="era2"]>0)/ length(mbeta$value[mbeta$variable=="era2"])
sum(mbeta$value[mbeta$variable=="era3"]>0)/ length(mbeta$value[mbeta$variable=="era3"])

# calculate pairwise overlaps in slopes and intercepts
slope_overlap = overlapping::overlap(x = list(slope1 = mu_beta[,1],slope2=mu_beta[,2],slope3=mu_beta[,3]))
int_overlap = overlapping::overlap(x = list(int1 = mu_alpha[ , 1],int2=mu_alpha[ , 2],int3=mu_alpha[ , 3]))


slopes2 <- ggplot(mbeta, aes(x = value, fill = variable)) +
  theme_bw() +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(cb[2], cb[3], cb[4]),
                    labels=c("1965-1988", "1989-2013", "2014-2019")) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope (scaled anomaly)",
       y = "Posterior density") +
  theme(legend.title = element_blank(), legend.position = 'bottom',
        legend.direction = "horizontal", legend.key.size = unit(4, 'mm'))
print(slopes2)

# reset era names for plotting
dat3$era.labs <- factor(dat3$era, labels = c("1965-1988", "1989-2013", "2014-2019"))

# and order spp
dat3$plot.order <- ifelse(dat3$species=="Pink-odd", 1,
                          ifelse(dat3$species=="Pink-even", 2,
                                 ifelse(dat3$species=="Sockeye", 3, 4)))

dat3$species <- reorder(dat3$species, dat3$plot.order)

scatter2 <- ggplot(dat3) +
  aes(x = pdo, y = catch, color = species) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap( ~era.labs) +
  scale_color_manual(values=c(cb[2], cb[7], cb[6], cb[4])) +
  theme_bw() + ylab("Catch anomaly") + xlab("PDO (Nov-Mar, 2-yr running mean)") +
  theme(legend.title = element_blank(), legend.position = 'bottom')

print(scatter2)


#############
# and one-year!
## Prep data -----------------------------------------------
dat3 <- na.omit(raw.dat)
dat3 <- plyr::ddply(dat3, .(species), transform, pdo = scale(PDO1))

# dat3 <- plyr::ddply(dat3, .(species), transform, sst = (9/5)*SST3+32) # changing to raw ºF!

dat3$era <- as.factor(dat3$era)

## Find start and end indices for each species
n <- as.numeric(dat3$species)
int.end <- which(diff(n) != 0)
int.start <- int.end + 1
end <- c(int.end, length(n))
start <- c(1, int.start)

## Data for Stan models
dat3_stan <- list(y = dat3$catch,
                  x1 = dat3$pdo,
                  y_start = start,
                  y_end = end,
                  n_species = length(unique(dat3$species)),
                  N = nrow(dat3),
                  era1 = ifelse(dat3$catch.year <= 1988, 1, 0),
                  era2 = ifelse(dat3$catch.year %in% 1989:2013, 1, 0),
                  era3 = ifelse(dat3$catch.year >= 2014, 1, 0))

## 3 era: hierarchical --> model for manuscript ----------------
era3_hier_arm <- stan_glmer(catch ~ era + pdo + pdo:era + (era + pdo + pdo:era | species),
                            data = dat3,
                            chains = 4, cores = 4, thin = 1,
                            warmup = 1000, iter = 4000, refresh = 0,
                            adapt_delta = 0.99,
                            prior = normal(location = 0, scale = 5, autoscale = FALSE),
                            prior_intercept = normal(location = 0, scale = 5, autoscale = FALSE),
                            prior_aux = student_t(df = 3, location = 0, scale = 5, autoscale = FALSE),
                            prior_covariance = decov(regularization = 1,
                                                     concentration = 1,
                                                     shape = 1, scale = 1))

# do predictions for recent era
newdata = expand.grid("era"=unique(dat3$era),
                      pdo=c(mean(dat3$pdo[which(dat3$era=="era3")]),
                            mean(dat3$pdo[which(dat3$era=="era3")])+1), species=unique(dat3$species)) %>%
  dplyr::filter(era=="era3")
pred = posterior_predict(era3_hier_arm, newdata=newdata)
newdata$mean_pred = round(apply(pred,2,mean), 3)

fixef(era3_hier_arm)
ranef(era3_hier_arm)
coef(era3_hier_arm)
era3_hier_arm$covmat
print(era3_hier_arm)

mu_beta  <- as.matrix(era3_hier_arm, pars = c("pdo", "eraera2:pdo", "eraera3:pdo"))
coef_beta <- data.frame(coef = "Slope",
                        era1 = mu_beta[ , 1],
                        era2 = mu_beta[ , 1] + mu_beta[ , 2],
                        era3 = mu_beta[ , 1] + mu_beta[ , 3])
mu_alpha <- as.matrix(era3_hier_arm, pars = c("(Intercept)", "eraera2", "eraera3"))
coef_alpha <- data.frame(coef = "Intercept",
                         era1 = mu_alpha[ , 1],
                         era2 = mu_alpha[ , 1] + mu_alpha[ , 2],
                         era3 = mu_alpha[ , 1] + mu_alpha[ , 3])
mbeta  <- reshape2::melt(coef_beta, id.vars = "coef")
malpha <- reshape2::melt(coef_alpha, id.vars = "coef")
mdf_hier <- rbind(mbeta, malpha)


era_mean_slopes <- mbeta %>%
  group_by(variable) %>%
  summarize(mean=mean(value))

era_mean_slopes[3,2] - era_mean_slopes[2,2]
era_mean_slopes[3,2] - era_mean_slopes[1,2]

# and % positive / negative slopes by era!
sum(mbeta$value[mbeta$variable=="era1"]>0)/ length(mbeta$value[mbeta$variable=="era1"])
sum(mbeta$value[mbeta$variable=="era2"]>0)/ length(mbeta$value[mbeta$variable=="era2"])
sum(mbeta$value[mbeta$variable=="era3"]>0)/ length(mbeta$value[mbeta$variable=="era3"])

# calculate pairwise overlaps in slopes and intercepts
slope_overlap = overlapping::overlap(x = list(slope1 = mu_beta[,1],slope2=mu_beta[,2],slope3=mu_beta[,3]))
int_overlap = overlapping::overlap(x = list(int1 = mu_alpha[ , 1],int2=mu_alpha[ , 2],int3=mu_alpha[ , 3]))


slopes1 <- ggplot(mbeta, aes(x = value, fill = variable)) +
  theme_bw() +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(cb[2], cb[3], cb[4]),
                    labels=c("1965-1988", "1989-2013", "2014-2019")) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope (scaled anomaly)",
       y = "Posterior density") +
  theme(legend.title = element_blank(), legend.position = 'bottom',
        legend.direction = "horizontal", legend.key.size = unit(4, 'mm'))
print(slopes1)

# reset era names for plotting
dat3$era.labs <- factor(dat3$era, labels = c("1965-1988", "1989-2013", "2014-2019"))

# and order spp
dat3$plot.order <- ifelse(dat3$species=="Pink-odd", 1,
                          ifelse(dat3$species=="Pink-even", 2,
                                 ifelse(dat3$species=="Sockeye", 3, 4)))
dat3$species <- reorder(dat3$species, dat3$plot.order)

scatter1 <- ggplot(dat3) +
  aes(x = pdo, y = catch, color = species) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap( ~era.labs) +
  scale_color_manual(values=c(cb[2], cb[7], cb[6], cb[4])) +
  theme_bw() + ylab("Catch anomaly") + xlab("PDO (Nov-Mar)") +
  theme(legend.title = element_blank(), legend.position = 'bottom')

print(scatter1)

fill <- ggplot() + theme_void()

png("figs/era-specific catches and 1-y and 2-yrPDO.png", 8, 6.5, units='in', res=300)
ggpubr::ggarrange(ggpubr::ggarrange(fill, fill, nrow=1, widths = c(1, 0.7)),
                  ggpubr::ggarrange(scatter2, slopes2, nrow=1, widths = c(1, 0.7)),
                  ggpubr::ggarrange(fill, fill, nrow=1, widths = c(1, 0.7)),
                  ggpubr::ggarrange(scatter1, slopes1, nrow=1, widths = c(1, 0.7)),
                  nrow=4, heights = c(0.15, 1, 0.1, 1),
                  labels=c(NA, "a) Two-year PDO mean", NA, "b) One-year PDO mean"), vjust=0.3, hjust=0,
                  font.label = list(size=12))
dev.off()

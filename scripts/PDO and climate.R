library(gplots)
library(tidyr)
library(overlapping)
library(tidyverse)

dat <- read.csv("data/climate.data.csv")

cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# excluding variables showing negative results (these are plotted separately)
# and limit to 1965-onwards
dat <- dat %>%
  select(-Papa, -GAK1.sal, -NDJFM.PDO) %>%
  filter(year > 1964)

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

# now plot relative to PDO!
colnames(dat)[c(8,4,5,6)] <- c("Sea surface height (PC1)",
                               "Wind stress (Pa)",
                               "Sea surface temp. (ºC)",
                             "North Pacific Index (mb)")

dat <- dat %>%
  select(-year, -FMA.ssh.COMBINED) %>%
  pivot_longer(cols=c(-NDJFM.PDO, -era), names_to = 'key')

dat$plot.era <- ifelse(dat$era==1, "1965-1988",
                       ifelse(dat$era==2, "1989-2013", "2014-2019"))

dat$key.order <- ifelse(dat$key=="North Pacific Index (mb)", 1,
                        ifelse(dat$key=="Sea surface temp. (ºC)", 2,
                        ifelse(dat$key=="Sea surface height (PC1)", 3, 4)))

dat$key <- reorder(dat$key, dat$key.order)

scatter <- ggplot(dat, aes(NDJFM.PDO, value, color=plot.era)) +
  theme_bw() +
  geom_point() +
  facet_wrap(~key, scales="free_y") +
  scale_color_manual(values=cb[2:4]) +
  geom_smooth(method="lm", se=F) +
  xlab("PDO Index (Nov-Mar)") +
  theme(legend.title = element_blank(), axis.title.y = element_blank(), legend.position = 'top')

# and analyze era relationships with stan models
library(rstan)
library(ggplot2)
library(plyr)
library(rstanarm)
library(bayesplot)

# rename pdo
names(dat)[2] <- "pdo"

# remove units for posterior plots

#era as function
dat$era <- as.factor(dat$era)

# remove units from labels
dat$key <- as.character(dat$key)

change <- grep("North P", dat$key)
dat$key[change] <- "North Pacific Index"

change <- grep("Wind", dat$key)
dat$key[change] <- "Wind stress"

change <- grep("height", dat$key)
dat$key[change] <- "Sea surface height"

change <- grep("temp", dat$key)
dat$key[change] <- "Sea surface temp."

dat$key <- factor(dat$key)
dat$key <- reorder(dat$key, dat$key.order)

## fit a model with era-specific intercepts and slopes

era_NPI_2 <- stan_glm(scale(value) ~ era + pdo + pdo:era,
                      data = dat[dat$key == "North Pacific Index", ],
                      chains = 4, cores = 4, thin = 1, seed=421,
                      warmup = 1000, iter = 4000, refresh = 0,
                      prior = normal(location = 0, scale = 5, autoscale = FALSE),
                      prior_intercept = normal(location = 0, scale = 5, autoscale = FALSE),
                      prior_aux = student_t(df = 3, location = 0, scale = 5, autoscale = FALSE))

era_stress_2 <- stan_glm(scale(value) ~ era + pdo + pdo:era,
                         data = dat[dat$key == "Wind stress", ],
                         chains = 4, cores = 4, thin = 1, seed=421,
                         warmup = 1000, iter = 4000, refresh = 0,
                         prior = normal(location = 0, scale = 5, autoscale = FALSE),
                         prior_intercept = normal(location = 0, scale = 5, autoscale = FALSE),
                         prior_aux = student_t(df = 3, location = 0, scale = 5, autoscale = FALSE))

era_SSH_2 <- stan_glm(scale(value) ~ era + pdo + pdo:era,
                      data = dat[dat$key == "Sea surface height", ],
                      chains = 4, cores = 4, thin = 1, seed=421,
                      warmup = 1000, iter = 4000, refresh = 0,
                      prior = normal(location = 0, scale = 5, autoscale = FALSE),
                      prior_intercept = normal(location = 0, scale = 5, autoscale = FALSE),
                      prior_aux = student_t(df = 3, location = 0, scale = 5, autoscale = FALSE))

era_SST_2 <- stan_glm(scale(value) ~ era + pdo + pdo:era,
                      data = dat[dat$key == "Sea surface temp.", ],
                      chains = 4, cores = 4, thin = 1, seed=421,
                      warmup = 1000, iter = 4000, refresh = 0,
                      prior = normal(location = 0, scale = 5, autoscale = FALSE),
                      prior_intercept = normal(location = 0, scale = 5, autoscale = FALSE),
                      prior_aux = student_t(df = 3, location = 0, scale = 5, autoscale = FALSE))

lst <- list(era_NPI_2, era_stress_2, era_SSH_2, era_SST_2)

lst <- lapply(lst, function(x) {
  beta <- as.matrix(x, pars = c("(Intercept)", "era2", "era3"))
  data.frame(key = unique(x$data$key),
             era1 = beta[ , 1],
             era2 = beta[ , 1] + beta[ , 2],
             era3 = beta[ , 1] + beta[ , 3])
})
coef_indv_arm <- plyr::rbind.fill(lst)
mdf_indv_arm <- reshape2::melt(coef_indv_arm, id.vars = "key")


for(i in 1:length(unique(coef_indv_arm$key))) {

  sub = dplyr::filter(coef_indv_arm, key == unique(coef_indv_arm$key)[i])
  # calculate pairwise overlaps in slopes and intercepts
  int_overlap = overlapping::overlap(x = list(int1 = sub$era1,int2=sub$era2,int3=sub$era3))
  saveRDS(int_overlap$OV,file=paste0("output/", sub$key[1], "_climate_int_overlap.rds"))

}

## extract slopes
lst <- list(era_NPI_2, era_stress_2, era_SSH_2, era_SST_2)
lst.slope <- lapply(lst, function(x) {
  beta <- as.matrix(x, pars = c("pdo", "era2:pdo", "era3:pdo"))
  data.frame(key = unique(x$data$key),
             era1 = beta[ , 1],
             era2 = beta[ , 1] + beta[ , 2],
             era3 = beta[ , 1] + beta[ , 3])
})
coef_slope <- plyr::rbind.fill(lst.slope)
mdf_slope <- reshape2::melt(coef_slope, id.vars = "key")

# calculate overlap in slopes

for(i in 1:length(unique(coef_slope$key))) {
 
  sub = dplyr::filter(coef_slope, key == unique(coef_slope$key)[i])
  # calculate pairwise overlaps in slopes and intercepts
  int_overlap = overlapping::overlap(x = list(slope1 = sub$era1,slope2=sub$era2,slope3=sub$era3))
  saveRDS(int_overlap$OV,file=paste0("output/", sub$key[1], "_climate_slope_overlap.rds"))
  

}

int_tab <- plyr::ddply(mdf_indv_arm, .(key, variable), summarize,
                       mean = mean(value),
                       lower95 = quantile(value, probs = 0.025),
                       upper95 = quantile(value, probs = 0.975))

slope_tab <- plyr::ddply(mdf_slope, .(key, variable), summarize,
                         mean = mean(value),
                         lower95 = quantile(value, probs = 0.025),
                         upper95 = quantile(value, probs = 0.975))


# get intercept / slope differences in units of SD
slope_tab[3,3]-slope_tab[1,3] # NPI era3 - era1
slope_tab[3,3]-slope_tab[2,3] # NPI era3 - era2

slope_tab[6,3]-slope_tab[4,3] # sst era3 - era1
slope_tab[6,3]-slope_tab[5,3] # sst era3 - era2

slope_tab[9,3]-slope_tab[7,3] # ssh era3 - era1
slope_tab[9,3]-slope_tab[8,3] # ssh era3 - era2

slope_tab[12,3]-slope_tab[10,3] # stress era3 - era1
slope_tab[12,3]-slope_tab[11,3] # stress era3 - era2


int_tab[3,3]-int_tab[1,3] # NPI era3 - era1
int_tab[3,3]-int_tab[2,3] # NPI era3 - era2

int_tab[6,3]-int_tab[4,3] # sst era3 - era1
int_tab[6,3]-int_tab[5,3] # sst era3 - era2

int_tab[9,3]-int_tab[7,3] # ssh era3 - era1
int_tab[9,3]-int_tab[8,3] # ssh era3 - era2

int_tab[12,3]-int_tab[10,3] # stress era3 - era1
int_tab[12,3]-int_tab[11,3] # stress era3 - era2


write.csv(int_tab, "output/climate intercept table.csv")

write.csv(slope_tab, "output/climate slope table.csv")


# make a combined plot
library(ggpubr)


int <- ggplot(mdf_indv_arm, aes(x = value, fill = variable)) +
  theme_bw() +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(cb[2], cb[3], cb[4]), labels=c("1964-1988", "1989-2013", "2014-2019")) +
  theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Intercept (scaled anomaly)",
       y = "Posterior density") +
  facet_wrap( ~ key, scales="free")


# make a combined plot

tiff("figs/era-specific PDO and climate - vertical.tiff", 4, 7, units='in', res=300)
# png("era-specific PDO and climate - vertical.png", 8, 7, units='in', res=300)
ggarrange(scatter, int, ncol=1, nrow=2, labels=c("e)", "f)"),
          label.x = 0.05, label.y = 0.95, font.label = "plain")
dev.off()



# plot slopes for SI
slope <- ggplot(mdf_slope, aes(x = value, fill = variable)) +
  theme_bw() +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(cb[2], cb[3], cb[4]), labels=c("1964-1988", "1989-2013", "2014-2019")) +
  theme(legend.title = element_blank()) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope (scaled anomaly)",
       y = "Posterior density") +
  facet_wrap( ~ key, scales="free")
print(slope)

ggsave("figs/SI - era-specific climate slopes on PDO.png", width=6, height=4, units='in')

## Bayesian model diagnostics
post_era_NPI_2    <- as.array(era_NPI_2)
post_era_stress_2 <- as.array(era_stress_2)
post_era_SSH_2    <- as.array(era_SSH_2)
post_era_SST_2    <- as.array(era_SST_2)

mcmc_trace(post_era_NPI_2)
mcmc_trace(post_era_stress_2)
mcmc_trace(post_era_SSH_2)
mcmc_trace(post_era_SST_2)

mcmc_areas(post_era_NPI_2)
mcmc_areas(post_era_stress_2)
mcmc_areas(post_era_SSH_2)
mcmc_areas(post_era_SST_2)

range(summary(era_NPI_2[["stanfit"]])[["summary"]][ , "n_eff"])
range(summary(era_NPI_2[["stanfit"]])[["summary"]][ , "Rhat"])
range(summary(era_stress_2[["stanfit"]])[["summary"]][ , "n_eff"])
range(summary(era_stress_2[["stanfit"]])[["summary"]][ , "Rhat"])
range(summary(era_SSH_2[["stanfit"]])[["summary"]][ , "n_eff"])
range(summary(era_SSH_2[["stanfit"]])[["summary"]][ , "Rhat"])
range(summary(era_SST_2[["stanfit"]])[["summary"]][ , "n_eff"])
range(summary(era_SST_2[["stanfit"]])[["summary"]][ , "Rhat"])

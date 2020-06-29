library(tidyverse)
library(overlapping)

# run/plot the negative results for changing climate - PDO intercepts
# GAK1 salinity and Papa advection

dat <- read.csv("data/climate.data.csv", row.names = 1)

dat <- dat %>%
  select(Papa, GAK1.sal, era)

dat$year <- as.numeric(rownames(dat))

# limit to 1965-2019
dat <- dat %>%
  filter(year > 1964)

# add SWFSC PDO values
dat2 <- read.csv("data/salmon.and.SWFSC.PDO.csv")

dat2 <- dat2 %>%
  group_by(Year) %>%
  summarise(NDJFM.PDO=mean(PDO1))

names(dat2)[1] <- "year"

dat <- left_join(dat, dat2)

# now plot relative to PDO!
colnames(dat)[c(1,2)] <- c("Papa index (ºN, Dec-Feb)",
                           "GAK1 salinity (psu, Feb-Apr)")
dat <- dat %>%
  gather(key, value, -NDJFM.PDO, -era)

dat$plot.era <- ifelse(dat$era==1, "1964-1988",
                       ifelse(dat$era==2, "1989-2013", "2014-2019"))

# load colorblind palette
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

scatter <- ggplot(dat, aes(NDJFM.PDO, value, color=plot.era)) +
  theme_bw() +
  geom_point() +
  facet_wrap(~key, scales="free_y") +
  scale_color_manual(values=cb[2:4]) +
  geom_smooth(method="lm", se=F) +
  xlab("PDO Index (Nov-Mar)") +
  theme(legend.title = element_blank(), axis.title.y = element_blank())

# rename pdo
names(dat)[2] <- "pdo"

# remove units for posterior plots

#era as function
dat$era <- as.factor(dat$era)

# remove units from labels
dat$key <- as.character(dat$key)

change <- grep("GAK", dat$key)
dat$key[change] <- "GAK1 salinity (Feb-Apr)"

change <- grep("Papa", dat$key)
dat$key[change] <- "Papa index (Dec-Feb)"

dat$key <- factor(dat$key)

## fit a model with era-specific intercepts and slopes
library(rstan)
library(ggplot2)
library(plyr)
library(rstanarm)
library(bayesplot)

era_GAK_2 <- stan_glm(scale(value) ~ era + pdo + pdo:era,
                      data = dat[dat$key == "GAK1 salinity (Feb-Apr)", ],
                      chains = 4, cores = 4, thin = 1, seed=421,
                      warmup = 1000, iter = 4000, refresh = 0,
                      prior = normal(location = 0, scale = 5, autoscale = FALSE),
                      prior_intercept = normal(location = 0, scale = 5, autoscale = FALSE),
                      prior_aux = student_t(df = 3, location = 0, scale = 5, autoscale = FALSE))

era_Papa_2 <- stan_glm(scale(value) ~ era + pdo + pdo:era,
                       data = dat[dat$key == "Papa index (Dec-Feb)", ],
                       chains = 4, cores = 4, thin = 1, seed=421,
                       warmup = 1000, iter = 4000, refresh = 0,
                       prior = normal(location = 0, scale = 5, autoscale = FALSE),
                       prior_intercept = normal(location = 0, scale = 5, autoscale = FALSE),
                       prior_aux = student_t(df = 3, location = 0, scale = 5, autoscale = FALSE))


lst <- list(era_GAK_2, era_Papa_2)

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
  # i <- 1
  sub = dplyr::filter(coef_indv_arm, key == unique(coef_indv_arm$key)[i])
  # calculate pairwise overlaps in slopes and intercepts
  int_overlap = overlapping::overlap(x = list(int1 = sub$era1,int2=sub$era2,int3=sub$era3))
  saveRDS(int_overlap$OV,file=paste0("output/", sub$key[1], "_climate_int_overlap.rds"))
}

## extract slopes
lst <- list(era_GAK_2, era_Papa_2)
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
  # i <- 2
  sub = dplyr::filter(coef_slope, key == unique(coef_slope$key)[i])
  # calculate pairwise overlaps in slopes and intercepts
  int_overlap = overlapping::overlap(x = list(slope1 = sub$era1,slope2=sub$era2,slope3=sub$era3))
  saveRDS(int_overlap$OV,file=paste0("output/", sub$key[1], "_climate_slope_overlap.rds"))
}

int_tab_neg <- plyr::ddply(mdf_indv_arm, .(key, variable), summarize,
                           mean = mean(value),
                           lower95 = quantile(value, probs = 0.025),
                           upper95 = quantile(value, probs = 0.975))

slope_tab_neg <- plyr::ddply(mdf_slope, .(key, variable), summarize,
                             mean = mean(value),
                             lower95 = quantile(value, probs = 0.025),
                             upper95 = quantile(value, probs = 0.975))

int_tab_neg[3,3]-int_tab_neg[1,3] # gak1 era3 - era1
int_tab_neg[3,3]-int_tab_neg[2,3] # gak1 era3 - era2

int_tab_neg[6,3]-int_tab_neg[4,3] # Papa era3 - era1
int_tab_neg[6,3]-int_tab_neg[5,3] # Papa era3 - era2

slope_tab_neg[3,3]-slope_tab_neg[1,3] # gak1 era3 - era1
slope_tab_neg[3,3]-slope_tab_neg[2,3] # gak1 era3 - era2

slope_tab_neg[6,3]-slope_tab_neg[4,3] # Papa era3 - era1
slope_tab_neg[6,3]-slope_tab_neg[5,3] # Papa era3 - era2

# plot intercepts
int <- ggplot(mdf_indv_arm, aes(x = value, fill = variable)) +
  theme_bw() +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(cb[2], cb[3], cb[4]), labels=c("1964-1988", "1989-2013", "2014-2019")) +
  theme(legend.title = element_blank()) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Intercept (scaled anomaly)",
       y = "Posterior density") +
  facet_wrap( ~ key, scales="free")
print(int)

# plot slopes
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

# make a combined plot
png("figs/SI - era-specific PDO and GAK1-Papa.png", 6.5, 7, units='in', res=300)
ggpubr::ggarrange(scatter, int, slope, ncol=1, nrow=3, labels=c("a)", "b)", "c)"))
dev.off()



# make a combined plot
library(ggpubr)

png("figs/era-specific PDO - Papa and salinity.png", 6.5, 5, units='in', res=300)
ggarrange(scatter, int, ncol=1, nrow=2, labels=c("a)", "b)"))
dev.off()


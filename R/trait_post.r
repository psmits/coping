library(plyr)
library(coda)
library(arm)
library(stringr)
library(reshape2)
library(ggplot2)
library(rstan)
library(grid)
source('../R/coping_foo_post.r')
source('../data/data_dump/trait_info.data.R')
#
theme_set(theme_bw())
cbp <- c('#E69F00', '#56B4E9', '#009E73', '#F0E442', 
         '#0072B2', '#D55E00', '#CC79A7')
theme_update(axis.text = element_text(size = 10),
             axis.title = element_text(size = 20),
             legend.text = element_text(size = 15),
             legend.title = element_text(size = 16),
             legend.key.size = unit(1, 'cm'),
             strip.text = element_text(size = 7))

nsim <- 1000

post <- list.files('../data/mcmc_out', 
                   pattern = 'trait_categorical',
                   full.names = TRUE)
fit <- read_stan_csv(post[1:4])
ext <- extract(fit, permuted = TRUE)
# sample, time, response class, covariate
# TO DO i have already simulated from the posterior predictive

# got to do some plots

# change in ratio over time


# to deal with making run over all samples, split by cohort
softmax <- function(beta) {
  tot <- sum(exp(beta))
  exp(beta) / tot
}
oo <- list()
for(ii in seq(C)) {
  oo[[ii]] <- ext$hold[1, cohort == ii, ]
}
hold.prob <- llply(oo, function(x) t(apply(x, 1, softmax)))
hold.prob <- laply(hold.prob, colMeans)
colnames(hold.prob) <- c('arb', 'grnd', 'scn')
hold.prob <- melt(hold.prob)
prob.seq <- ggplot(hold.prob, aes(x = Var1, y = value, fill = Var2)) 
prob.seq <- prob.seq + geom_bar(stat = 'identity', position = 'stack')
prob.seq <- prob.seq + labs(x = 'Time', y = 'Probability of occurrence')
ggsave(plot = prob.seq, filename = '../doc/figure/occurrence_prob.png')


# effect of trait size on probability of occurrence over time
eff.mean.cohort <- eff.sd.cohort <- list()
for(ii in seq(C)) {
  eff.mean.cohort[[ii]] <- colMeans(ext$beta[, ii, , 2])
  eff.sd.cohort[[ii]] <- apply(ext$beta[, ii, , 2], 2, sd)
}
eff.mean.cohort <- Reduce(rbind, eff.mean.cohort)
eff.sd.cohort <- Reduce(rbind, eff.sd.cohort)
eff.cohort <- cbind(melt(eff.mean.cohort), melt(eff.sd.cohort)[, 3])
names(eff.cohort) <- c('time', 'response', 'mean', 'sd')
eff.cohort$time <- rep(seq(1:C), K)
eff.cohort <- data.frame(eff.cohort)
trait.eff <- ggplot(eff.cohort, aes(x = time, y = mean, 
                                    colour = factor(response)))
trait.eff <- trait.eff + geom_pointrange(aes(ymax = mean + sd, ymin = mean - sd))
trait.eff <- trait.eff + labs(x = 'Time', y = 'Effect of body size')
trait.eff <- trait.eff + scale_fill_discrete(name = 'Locomotor\ncategory')


mean.eff <- cbind(mean = colMeans(cbind(ext$beta_mu[, , 2], 0)), 
                  sd = apply(cbind(ext$beta_mu[, , 2], 0), 2, sd))
mean.eff <- cbind(mean.eff, response = seq(nrow(mean.eff)))
mean.eff <- cbind(rbind(mean.eff, mean.eff), 
                  time = c(rep(1, K), rep(C, K)))
mean.eff <- data.frame(mean.eff)
mean.eff$ymin <- mean.eff$mean - mean.eff$sd
mean.eff$ymax <- mean.eff$mean + mean.eff$sd
trait.eff <- trait.eff + geom_ribbon(data = mean.eff,
                                     mapping = aes(x = time, y = mean, 
                                                   ymax = ymax,
                                                   ymin = ymin,
                                                   fill = factor(response),
                                                   colour = NULL), 
                                     alpha = 0.3)
trait.eff <- trait.eff + scale_colour_discrete(name = 'Locomotor\ncategory')
ggsave(plot = trait.eff, filename = '../doc/figure/trait_eff.png')

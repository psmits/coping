library(plyr)
library(coda)
library(arm)
library(stringr)
library(reshape2)
library(ggplot2)
library(rstan)
library(grid)
source('../R/coping_foo_post.r')
source('../R/multiclass_roc.r')
source('../data/data_dump/trait_info.data.R')
#
theme_set(theme_bw())
cbp <- c('#E69F00', '#56B4E9', '#009E73', '#F0E442', 
         '#0072B2', '#D55E00', '#CC79A7')
theme_update(axis.text = element_text(size = 15),
             axis.title = element_text(size = 20),
             legend.text = element_text(size = 17),
             legend.title = element_text(size = 20),
             legend.key.size = unit(1, 'cm'),
             strip.text = element_text(size = 20))

nsim <- 1000


# bring in the mcmc results
post <- list.files('../data/mcmc_out', 
                   pattern = 'trait_',
                   full.names = TRUE)
fit <- read_stan_csv(post)
ext <- extract(fit, permuted = TRUE)
# matrices are structured: sample, time, response class, covariate

roc.posterior <- replicate(nsim, roc.dist(ext, y, cohort), simplify = FALSE)
roc.posterior <- data.frame(Reduce(rbind, roc.posterior))
roc.melt <- melt(roc.posterior)
roc.melt$variable <- factor(roc.melt$variable, 
                            levels = rev(paste0('X', seq(C - 1))))
roc.plot <- ggplot(roc.melt, aes(x = variable, y = value))
roc.plot <- roc.plot + geom_violin() + geom_boxplot(width = 0.1)
roc.plot <- roc.plot + labs(x = 'Cohort', y = 'AUC')
ggsave(plot = roc.plot, filename = '../doc/figure/roc_dist.png',
       width = 8, height = 6, dpi = 600)



# effect of trait size on probability of occurrence over time
eff.mean.cohort <- eff.q1.cohort <- eff.q9.cohort <- list()
for(ii in seq(K)) {
  eff.mean.cohort[[ii]] <- mean(ext$beta[, ii, 2])
  eff.q1.cohort[[ii]] <- quantile(ext$beta[, ii, 2], 0.1)
  eff.q9.cohort[[ii]] <- quantile(ext$beta[, ii, 2], 0.9)
}
eff.mean.cohort <- Reduce(rbind, eff.mean.cohort)
eff.q1.cohort <- Reduce(rbind, eff.q1.cohort)
eff.q9.cohort <- Reduce(rbind, eff.q9.cohort)

eff.cohort <- cbind(melt(eff.mean.cohort), 
                    low = melt(eff.q1.cohort)[, 3], 
                    high = melt(eff.q9.cohort)[, 3])
names(eff.cohort)[1:3] <- c('time', 'response', 'mean')
eff.cohort <- data.frame(eff.cohort)
eff.cohort$response <- as.factor(seq(K))
trait.eff <- ggplot(eff.cohort, aes(x = response, y = mean))
trait.eff <- trait.eff + geom_pointrange(aes(ymax = high, ymin = low))
trait.eff <- trait.eff + labs(x = 'Ecotype', y = 'Effect of body size')
trait.eff <- trait.eff + scale_colour_manual(values = cbp,
                                             name = 'Locomotor\ncategory')
ggsave(plot = trait.eff, filename = '../doc/figure/trait_eff.png',
       width = 8, height = 6, dpi = 600)



## change in probability as cohort
#colMeans(ext$intercept[, , 1])
#colMeans(ext$intercept[, , 2])
#colMeans(ext$intercept[, , 3])
#colMeans(ext$intercept[, , 4])
#colMeans(ext$intercept[, , 5])
#colMeans(ext$intercept[, , 6])

# effect of climate
clim.eff <- data.frame(res = rep(seq(K - 1), 2),
                       typ = rep(c('mean', 'range'), each = K - 1),
                       mid = c(apply(ext$alpha[, , 1], 2, mean), 
                               apply(ext$alpha[, , 2], 2, mean)),
                       low = c(apply(ext$alpha[, , 1], 2, 
                                     function(x) quantile(x, 0.1)), 
                               apply(ext$alpha[, , 2], 2, 
                                     function(x) quantile(x, 0.1))),
                       hgh = c(apply(ext$alpha[, , 1], 2, 
                                     function(x) quantile(x, 0.9)), 
                               apply(ext$alpha[, , 2], 2, 
                                     function(x) quantile(x, 0.9))))
clim.plot <- ggplot(clim.eff, aes(x = res, y = mid))
clim.plot <- clim.plot + geom_pointrange(mapping = aes(ymax = hgh, ymin = low))
clim.plot <- clim.plot + facet_grid(. ~ typ, switch = 'x')
ggsave(plot = clim.plot, filename = '../doc/figure/clim_eff.png',
       width = 8, height = 6, dpi = 600)

library(plyr)
library(coda)
library(arm)
library(stringr)
library(reshape2)
library(ggplot2)
library(rstan)
library(grid)
library(pROC)
library(Metrics)
source('../R/multiclass_roc.r')
source('../data/data_dump/trait_info.data.R')
#
#theme_set(theme_bw())
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
                   full.names = TRUE)
fit <- read_stan_csv(post[1:4])

ext <- extract(fit, permuted = TRUE)
# matrices are structured: sample, time, response class, covariate

# now roc is much much easier because binary instead of multi




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

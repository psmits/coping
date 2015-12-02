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
theme_update(axis.text = element_text(size = 15),
             axis.title = element_text(size = 20),
             legend.text = element_text(size = 17),
             legend.title = element_text(size = 20),
             legend.key.size = unit(1, 'cm'),
             strip.text = element_text(size = 20))

nsim <- 1000

post <- list.files('../data/mcmc_out', 
                   pattern = 'trait_cat_',
                   full.names = TRUE)
fit <- read_stan_csv(post[1:4])
ext <- extract(fit, permuted = TRUE)

#rh <- ncol(summary(fit)[[1]])
#(summary(fit)[[1]][, rh] < 1.1)

# sample, time, response class, covariate
# i have already simulated from the posterior predictive


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
oo <- oo[!(laply(llply(oo, dim), is.null))]
hold.prob <- llply(oo, function(x) t(apply(x, 1, softmax)))
hold.prob <- laply(hold.prob, colMeans)
colnames(hold.prob) <- c('arb', 'fos', 'grnd', 'scn')
hold.prob <- melt(hold.prob)
prob.seq <- ggplot(hold.prob, aes(x = Var1, y = value, fill = factor(Var2)))
prob.seq <- prob.seq + geom_bar(stat = 'identity', position = 'stack')
prob.seq <- prob.seq + scale_x_reverse()
prob.seq <- prob.seq + labs(x = 'Time', y = 'Probability of occurrence')
prob.seq <- prob.seq + scale_fill_manual(values = cbp, 
                                         name = 'Locomotor\ncategory')
ggsave(plot = prob.seq, filename = '../doc/figure/occurrence_prob.png',
       width = 8, height = 6, dpi = 600)


# effect of trait size on probability of occurrence over time
eff.mean.cohort <- eff.q1.cohort <- eff.q9.cohort <- list()
for(ii in seq(C)) {
  eff.mean.cohort[[ii]] <- colMeans(ext$beta[, ii, , 2])
  eff.q1.cohort[[ii]] <- apply(ext$beta[, ii, , 2], 2, 
                               function(x) quantile(x, 0.1))
  eff.q9.cohort[[ii]] <- apply(ext$beta[, ii, , 2], 2, 
                               function(x) quantile(x, 0.9))
}
eff.mean.cohort <- Reduce(rbind, eff.mean.cohort)
eff.q1.cohort <- Reduce(rbind, eff.q1.cohort)
eff.q9.cohort <- Reduce(rbind, eff.q9.cohort)
colnames(eff.mean.cohort) <- c('arb', 'fos', 'grnd', 'scn')

eff.cohort <- cbind(melt(eff.mean.cohort), 
                    low = melt(eff.q1.cohort)[, 3], 
                    high = melt(eff.q9.cohort)[, 3])
names(eff.cohort)[1:3] <- c('time', 'response', 'mean')
eff.cohort$time <- rep(seq(1:C), K)
eff.cohort <- data.frame(eff.cohort)
eff.cohort$response <- as.factor(eff.cohort$response)
trait.eff <- ggplot(eff.cohort, aes(x = time, y = mean, 
                                    colour = response))
trait.eff <- trait.eff + geom_pointrange(aes(ymax = high, ymin = low))
trait.eff <- trait.eff + labs(x = 'Time', y = 'Effect of body size')
trait.eff <- trait.eff + scale_colour_manual(values = cbp,
                                             name = 'Locomotor\ncategory')

mean.eff <- cbind(mean = colMeans(cbind(ext$beta_mu[, , 2], 0)))
mean.eff <- data.frame(mean.eff)
mean.eff <- cbind(mean.eff, response = c('arb', 'fos', 'grnd', 'scn'))
mean.eff$ymin <- c(apply(ext$beta_mu[, , 2], 2, 
                         function(x) quantile(x, 0.1)), 0)
mean.eff$ymax <- c(apply(ext$beta_mu[, , 2], 2, 
                         function(x) quantile(x, 0.9)), 0)
mean.eff <- cbind(rbind(mean.eff, mean.eff), 
                  time = c(rep(1, K), rep(C, K)))
trait.eff <- trait.eff + geom_ribbon(data = mean.eff,
                                     mapping = aes(x = time, y = mean, 
                                                   ymax = ymax,
                                                   ymin = ymin,
                                                   fill = response,
                                                   colour = NULL), 
                                     alpha = 0.2)
trait.eff <- trait.eff + geom_line(data = mean.eff,
                                   mapping = aes(x = time, y = mean, 
                                                 colour = response))
trait.eff <- trait.eff + scale_fill_manual(values = cbp,
                                           name = 'Locomotor\ncategory')
trait.eff <- trait.eff + scale_x_reverse()
ggsave(plot = trait.eff, filename = '../doc/figure/trait_eff.png',
       width = 8, height = 6, dpi = 600)


# effect of diet on occurrence of trait
#   order: arb, foss, grnd, scn
cb <- cbind(ext$beta_mu[, , 1], 0)
hb <- cbind(ext$beta_mu[, , 3], 0)
ib <- cbind(ext$beta_mu[, , 4], 0)
ob <- cbind(ext$beta_mu[, , 5], 0)
diet.eff <- list(data.frame(cb, cov = 'carn'), 
                 data.frame(cb + hb, cov = 'herb'),
                 data.frame(cb + ib, cov = 'insect'),
                 data.frame(cb + ob, cov = 'omni'))
diet.eff <- llply(diet.eff, function(x) {
                  names(x) <- c('arb', 'foss', 'grnd', 'scn')
                  x})
diet.diff <- llply(diet.eff, function(x) pairwise.diffs(x[, 1:4]))
diet.diff <- Map(function(x, y) cbind(x, y),
                 x = diet.diff, 
                 y = c('carn', 'herb', 'insect', 'omni'))
diet.diff <- Reduce(rbind, diet.diff)
diet.diff <- melt(diet.diff)

relab.x <- scale_x_discrete(labels = 
                            c('beta[arb] - beta[foss]' = 
                              expression(beta[arb]-beta[foss]), 
                              'beta[arb] - beta[grnd]' = 
                              expression(beta[arb]-beta[grnd]), 
                              'beta[arb] - beta[scn]' = 
                              expression(beta[arb]-beta[scn]), 
                              'beta[foss] - beta[grnd]' = 
                              expression(beta[foss]-beta[grnd]),
                              'beta[foss] - beta[scn]' = 
                              expression(beta[foss]-beta[scn]),
                              'beta[grnd] - beta[scn]' = 
                              expression(beta[grd]-beta[scn])))

dietgg <- ggplot(diet.diff, aes(x = variable, y = value))
dietgg <- dietgg + geom_hline(yintercept = 0, colour = 'grey', size = 1.5)
dietgg <- dietgg + geom_violin() + geom_boxplot(width = 0.1)
dietgg <- dietgg + facet_grid(y ~ .) + relab.x
dietgg <- dietgg + labs(x = 'Pairwise comparison', y = 'Estimate')
ggsave(plot = dietgg, filename = '../doc/figure/diet_eff.png',
       width = 8, height = 8, dpi = 600)


# effect of climate on (base) occurrence of trait
#   there are 4 different estimates (k = 1,2,3,4)
#   group-level trait

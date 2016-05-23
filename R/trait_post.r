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
source('../R/read_one_stan.r')
source('../R/sim_from_model.r')
source('../data/data_dump/trait_info.data.R')
sight.implied <- sight
source('../data/data_dump/trait_w_gaps.data.R')
sight.obs <- sight
#
nsim <- 1000

###########
# full Bayes
post <- list.files('../data/mcmc_out', pattern = '[0-9]', full.names = TRUE)

# horseshoe priors
fit <- read_stan_csv(post[5:8])
ext <- extract(fit, permuted = TRUE)
# sample, group level, individual level

# now for beta
bytime <- list()
for(ii in seq(T)) {
  byindiv <- list()
  for(jj in seq(D)) {
    byindiv[[jj]] <- quantile(ext$beta[, ii, jj], 
                              c(0.10, 0.25, 0.5, 0.75, 0.90))
    names(byindiv[[jj]]) <- c('low', 'lowmid', 'med', 'highmed', 'high')
  }
  byindiv <- Reduce(rbind, byindiv)
  bytime[[ii]] <- data.frame(byindiv, coef = seq(D), time = ii)
}
for(jj in seq(length(T))) {
  for(ii in seq(from = 2, to = D)) {
    bytime[[jj]][ii, 1:5] <- bytime[[jj]][1, 1:5] + bytime[[jj]][ii, 1:5]
  }
}

melted <- Reduce(rbind, bytime)
melted$coef <- mapvalues(melted$coef, unique(melted$coef),
                         c('arboreal/carnivore', 'mass', 'herb', 'insect', 
                           'omni', 'digit', 'foss', 'planti', 'scan', 'unguli'))
melted$coef <- factor(melted$coef, 
                      levels = c('arboreal/carnivore', 'mass', 'herb', 
                                 'insect', 'omni', 'digit', 'foss', 'planti', 
                                 'scan', 'unguli'))
melted[, 1:5] <- apply(melted[, 1:5], 2, invlogit)  # probability scale
# iter, indiv, group, value
beta.plot <- ggplot(melted, aes(x = time, y = med))
beta.plot <- beta.plot + geom_hline(yintercept = 0.5, 
                                    colour = 'grey', size = 0.5)
beta.plot <- beta.plot + geom_line()
beta.plot <- beta.plot + geom_linerange(aes(ymax = high, ymin = low))
beta.plot <- beta.plot + facet_wrap(~ coef, ncol = 2)
beta.plot <- beta.plot + labs(x = 'Time (My)', y = 'Probability')



# now for gamma
byindiv <- list()
for(ii in seq(D)) {
  bygroup <- list()
  for(jj in seq(U)) {
    bygroup[[jj]] <- quantile(ext$gamma[, jj, ii], 
                              c(0.10, 0.25, 0.5, 0.75, 0.90))
    names(bygroup[[jj]]) <- c('low', 'lowmid', 'med', 'highmed', 'high')
  }
  bygroup <- Reduce(rbind, bygroup)
  byindiv[[ii]] <- data.frame(bygroup, group = seq(U), indiv = ii)
}
melted <- Reduce(rbind, byindiv)
melted$group <- mapvalues(melted$group, unique(melted$group), 
                          c('intercept', 'mean temp', 'range temp', 
                            'phase 2', 'phase3'))
melted$group <- factor(melted$group,
                       levels = c('intercept', 'mean temp', 'range temp', 
                                  'phase 2', 'phase3'))
melted$indiv <- mapvalues(melted$indiv, unique(melted$indiv),
                          c('arboreal/carnivore', 'mass', 'herb', 
                            'insect', 'omni', 'digit', 'foss', 'planti', 
                            'scan', 'unguli'))
melted$indiv <- factor(melted$indiv,
                       levels = c('arboreal/carnivore', 'mass', 'herb', 
                                  'insect', 'omni', 'digit', 'foss', 'planti',
                                  'scan', 'unguli'))
# iter, indiv, group, value
gamma.plot <- ggplot(melted, aes(x = factor(indiv), y = med))
gamma.plot <- gamma.plot + geom_hline(yintercept = 0, colour = 'grey')
gamma.plot <- gamma.plot + geom_linerange(aes(ymax = highmed, ymin = lowmid), 
                                          size = 2)
gamma.plot <- gamma.plot + geom_pointrange(aes(ymax = high, ymin = low))
gamma.plot <- gamma.plot + facet_grid(group ~ .)
gamma.plot <- gamma.plot + coord_flip()
gamma.plot <- gamma.plot + labs(x = 'Coefficient estimate (log odds scale)',
                                y = 'Individual-level effect')

# sample, D, D
ext$Omega
#matrix(ext$Omega[1, , ], ncol = D)





############
## advi
#post <- list.files('../data/mcmc_out', pattern = 'advi',
#                   full.names = TRUE)
#
## horseshoe priors
#fit <- read_one_stan_csv(post[1])
#
## prior predictive simulation
#preds <- fit[which(str_detect(names(fit), 'pred*'))]
#pp <- alply(preds, 1, function(x) matrix(x, nrow = N, ncol = T))
#pp <- llply(pp, function(x) {x <- x[, T:1]
#            x})
#pp <- llply(pp, function(x) apply(x, 2, unlist))
#
#simout.horse <- model.simulation(N, T, pp[[1]])
#
## gammas
#gams <- fit[which(str_detect(names(fit), 'gamma*'))]
#
#st <- seq(from = 1, to = ncol(gams), by = 5)
#byU <- list()
#for(ii in seq(length(st))) {
#  byU[[ii]] <- gams[, seq(st[ii], length.out = 5)]
#}
#byU.q <- llply(byU, function(x) 
#               t(apply(x, 2, function(y) 
#                       quantile(y, c(.1, .25, .5, .75, .9)))))
#byU.q <- Reduce(rbind, Map(function(x, y) data.frame(x, d = y), 
#                           byU.q, seq(length(byU.q))))
#byU.q$coef <- rownames(byU.q)
#colnames(byU.q) <- c('low', 'lowmed', 'med', 'highmed', 'high', 'd', 'coef')
#
##byU.q$coef <- factor(rep(1:5, each = D))
### length is 10 --> individual-level traits
### each element is length 5 --> group-level traits
##gamma.plot <- ggplot(byU.q, aes(x = coef, y = med))
##gamma.plot <- gamma.plot + geom_pointrange(aes(ymax = high, ymin = low))
#
#
## model with sampling
#fit <- read_one_stan_csv(post[2])  # horseshoe
#
## prior predictive simulations
#ps <- rev(fit[which(str_detect(names(fit), 'p\\.[0-9]'))])
#
#preds <- fit[which(str_detect(names(fit), 'pred*'))]
#pp <- alply(preds, 1, function(x) matrix(x, nrow = N, ncol = T))
#pp <- llply(pp, function(x) {x <- x[, T:1]
#            x})
#pp <- llply(pp, function(x) apply(x, 2, unlist))
#
#simout.md <- model.simulation(N, T, pp[[1]], ps[1, ])
#
## gammas
#gams <- fit[which(str_detect(names(fit), 'gamma*'))]
#
#st <- seq(from = 1, to = ncol(gams), by = 5)
#byU <- list()
#for(ii in seq(length(st))) {
#  byU[[ii]] <- gams[, seq(st[ii], length.out = 5)]
#}
#byU.q <- llply(byU, function(x) 
#               apply(x, 2, function(y) 
#                     quantile(y, c(.25, .5, .75))))
## length is 10 --> individual-level traits
## each element is length 5 --> group-level traits

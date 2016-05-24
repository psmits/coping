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
samp <- sample(4000, nsim)

###########
# full Bayes
post <- list.files('../data/mcmc_out', pattern = '[0-9]', full.names = TRUE)

# horseshoe priors
fit <- read_stan_csv(post[5:8])
ext <- extract(fit, permuted = TRUE)
# sample, group level, individual level
# send ext$pred through the simulator
#   ntax <- N
#   ntime <- T
#   pred <- ext$pred[1, , ] #?
#   model.simulation(ntax, ntime, pred)
sim <- list()
for(ii in seq(nsim)) {
  sim[[ii]] <- model.simulation(N, T, ext$pred[samp[ii], , ])
}

# break the binary factors up
nn <- oo <- list()
diet <- move <- list()
for(ii in seq(T)) {
  cept <- ext$beta[, ii, 1]
  inv.cept <- invlogit(cept)

  wcept <- apply(ext$beta[ , ii, 3:10], 2, function(x) x + cept)
  inv.wcept <- apply(wcept, 2, invlogit)

  rel.wcept <- apply(inv.wcept, 2, function(x) x - inv.cept)

  nn[[ii]] <- cbind(inv.cept, inv.wcept)
  oo[[ii]] <- cbind(inv.cept, rel.wcept)
  
  diet[[ii]] <- nn[[ii]][, c(1:4)]

  move[[ii]] <- nn[[ii]][, c(1, 5:9)]
}
diet <- diet[-c(1, T)]
move <- move[-c(1, T)]


diet <- llply(diet, function(x) 
              apply(x, 2, function(y) 
                    quantile(y, c(0.1, 0.25, 0.5, 0.75, 0.9))))
diet <- llply(diet, function(x) {
              colnames(x) <- c('carnivore', 'herb', 'insect', 'omni')
              x})
diet <- llply(diet, t)
diet <- Map(function(x, y) data.frame(x, type = rownames(x), time = y), 
            diet, seq(length(diet)))
diet <- llply(diet, function(x) {
              names(x) <- c('low', 'lowmed', 'med', 'highmed', 'high', 
                            'type', 'time')
              x})
diet <- data.frame(Reduce(rbind, diet))
dietdupe <- diet
dietdupe$group <- dietdupe$type
dietdupe$type <- NULL

# plot of relative probability of occurrence based on diet
dietprob <- ggplot(diet, aes(x = time, y = med, fill = type))
dietprob <- dietprob + geom_ribbon(data = dietdupe, 
                                 aes(ymax = highmed, ymin = lowmed,
                                     fill = NULL, group = group), 
                                 alpha = 0.2, colour = 'grey')
dietprob <- dietprob + geom_ribbon(aes(ymax = high, ymin = low), 
                                   alpha = 0.25)
dietprob <- dietprob + geom_ribbon(aes(ymax = highmed, ymin = lowmed), 
                                   alpha = 0.5)
dietprob <- dietprob + facet_wrap(~ type)



move <- llply(move, function(x) 
              apply(x, 2, function(y) 
                    quantile(y, c(0.1, 0.25, 0.5, 0.75, 0.9))))
move <- llply(move, function(x) {
              colnames(x) <- c('arboreal', 'digit', 'foss', 
                               'planti', 'scan', 'unguli')
              x})
move <- llply(move, t)
move <- Map(function(x, y) data.frame(x, type = rownames(x), time = y), 
            move, seq(length(move)))
move <- llply(move, function(x) {
              names(x) <- c('low', 'lowmed', 'med', 'highmed', 'high', 
                            'type', 'time')
              x})
move <- data.frame(Reduce(rbind, move))
movedupe <- move
movedupe$group <- movedupe$type
movedupe$type <- NULL

# plot of relative probability of occurrence based on move
moveprob <- ggplot(move, aes(x = time, y = med, fill = type))
moveprob <- moveprob + geom_ribbon(data = movedupe, 
                                 aes(ymax = highmed, ymin = lowmed,
                                     fill = NULL, group = group), 
                                 alpha = 0.15)
moveprob <- moveprob + geom_ribbon(aes(ymax = high, ymin = low), 
                                   alpha = 0.25)
moveprob <- moveprob + geom_ribbon(aes(ymax = highmed, ymin = lowmed), 
                                   alpha = 0.5)
moveprob <- moveprob + facet_wrap(~ type)




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

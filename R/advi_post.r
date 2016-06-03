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



###########
# advi
post <- list.files('../data/mcmc_out', pattern = 'advi',
                   full.names = TRUE)

# horseshoe priors
fit <- read_one_stan_csv(post[1])

# prior predictive simulation
preds <- fit[which(str_detect(names(fit), 'pred*'))]
pp <- alply(preds, 1, function(x) matrix(x, nrow = N, ncol = T))
pp <- llply(pp, function(x) {x <- x[, T:1]
            x})
pp <- llply(pp, function(x) apply(x, 2, unlist))

simout.horse <- model.simulation(N, T, pp[[1]])

# gammas
gams <- fit[which(str_detect(names(fit), 'gamma*'))]

st <- seq(from = 1, to = ncol(gams), by = 5)
byU <- list()
for(ii in seq(length(st))) {
  byU[[ii]] <- gams[, seq(st[ii], length.out = 5)]
}
byU.q <- llply(byU, function(x) 
               t(apply(x, 2, function(y) 
                       quantile(y, c(.1, .25, .5, .75, .9)))))
byU.q <- Reduce(rbind, Map(function(x, y) data.frame(x, d = y), 
                           byU.q, seq(length(byU.q))))
byU.q$coef <- rownames(byU.q)
colnames(byU.q) <- c('low', 'lowmed', 'med', 'highmed', 'high', 'd', 'coef')

#byU.q$coef <- factor(rep(1:5, each = D))
## length is 10 --> individual-level traits
## each element is length 5 --> group-level traits
#gamma.plot <- ggplot(byU.q, aes(x = coef, y = med))
#gamma.plot <- gamma.plot + geom_pointrange(aes(ymax = high, ymin = low))


# model with sampling
fit <- read_one_stan_csv(post[2])  # horseshoe

# prior predictive simulations
ps <- rev(fit[which(str_detect(names(fit), 'p\\.[0-9]'))])

preds <- fit[which(str_detect(names(fit), 'pred*'))]
pp <- alply(preds, 1, function(x) matrix(x, nrow = N, ncol = T))
pp <- llply(pp, function(x) {x <- x[, T:1]
            x})
pp <- llply(pp, function(x) apply(x, 2, unlist))

simout.md <- model.simulation(N, T, pp[[1]], ps[1, ])

# gammas
gams <- fit[which(str_detect(names(fit), 'gamma*'))]

st <- seq(from = 1, to = ncol(gams), by = 5)
byU <- list()
for(ii in seq(length(st))) {
  byU[[ii]] <- gams[, seq(st[ii], length.out = 5)]
}
byU.q <- llply(byU, function(x) 
               apply(x, 2, function(y) 
                     quantile(y, c(.25, .5, .75))))
# length is 10 --> individual-level traits
# each element is length 5 --> group-level traits


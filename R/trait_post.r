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

#post <- list.files('../data/mcmc_out', pattern = '[0-9]', full.names = TRUE)
#post <- post[2:4]
#fit <- read_stan_csv(post)

# bring in the advi results
post <- list.files('../data/mcmc_out', pattern = 'advi',
                   full.names = TRUE)

# horseshoe priors
fit <- read_one_stan_csv(post[1])
# need to get gammas and betas
# also want Omega

preds <- fit[which(str_detect(names(fit), 'pred*'))]
pp <- alply(preds, 1, function(x) matrix(x, nrow = N, ncol = T))
pp <- llply(pp, function(x) {x <- x[, T:1]
            x})
pp <- llply(pp, function(x) apply(x, 2, unlist))

simout.horse <- model.simulation(N, T, pp[[1]])


# with sampling
fit <- read_one_stan_csv(post[2])  # horseshoe

ps <- rev(fit[which(str_detect(names(fit), 'p\\.[0-9]'))])

preds <- fit[which(str_detect(names(fit), 'pred*'))]
pp <- alply(preds, 1, function(x) matrix(x, nrow = N, ncol = T))
pp <- llply(pp, function(x) {x <- x[, T:1]
            x})
pp <- llply(pp, function(x) apply(x, 2, unlist))

simout.md <- model.simulation(N, T, pp[[1]], ps[1, ])

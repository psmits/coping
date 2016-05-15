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
source('../R/sim_from_model.r')
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
fit <- read_stan_csv(post)
#traceplot(fit, pars = 'lp__')

ext <- extract(fit, permuted = TRUE)
# matrices are structured: sample, time, response class, covariate


library(plyr)
library(coda)
library(arm)
library(stringr)
library(reshape2)
library(ggplot2)
library(rstan)
library(grid)
source('../R/coping_foo_post.r')
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
fit <- read_stan_csv(post)
ext <- extract(fit, permuted = TRUE)



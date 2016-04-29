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

# now roc is much much easier because binary instead of multi
d <- dim(ext$pred)
ntax <- d[2]
ntime <- d[3]
tt <- sample(length(ext$lp__), 1000)
test <- list()
for(ii in seq(length(tt))) {
  test[[ii]] <- model.simulation(ntax = ntax, ntime = ntime,
                                 p.mu = ext$p_mu[tt[ii]], 
                                 p.sigma = ext$p_sigma[tt[ii]], 
                                 p.eff = ext$p_norm[tt[ii], ],
                                 inter.mu = ext$intercept_mu[tt[ii]],
                                 inter.sigma = ext$sigma[tt[ii]], 
                                 inter.eff = ext$intercept[tt[ii], ])
}

# plot of posterior of estimated diversity plus mean line
div.est <- laply(test, function(x) colSums(x$z))
div.mean <- data.frame(x = seq(ntime), y = colMeans(div.est))
div.melt <- melt(div.est)
div.gg <- ggplot(div.melt, aes(x = Var2, y = value, group = Var1)) 
div.gg <- div.gg + geom_line(alpha = 0.1, size = 0.1)
div.gg <- div.gg + geom_line(data = div.mean, 
                             mapping = aes(x = x, y = y, group = NULL), 
                             colour = 'blue')

# plot of posterior of observed diversity plus mean line and empirical line
obs.est <- laply(test, function(x) colSums(x$z))
obs.mean <- data.frame(x = seq(ntime), y = colMeans(obs.est))
obs.emp <- data.frame(x = seq(ntime), y = colSums(sight))
obs.melt <- melt(obs.est)
obs.gg <- ggplot(obs.melt, aes(x = Var2, y = value, group = Var1)) 
obs.gg <- obs.gg + geom_line(alpha = 0.1, size = 0.1)
obs.gg <- obs.gg + geom_line(data = obs.mean, 
                             mapping = aes(x = x, y = y, group = NULL), 
                             colour = 'blue')
obs.gg <- obs.gg + geom_line(data = obs.emp, 
                             mapping = aes(x = x, y = y, group = NULL), 
                             colour = 'goldenrod')

# plot of posterior of predictor of presence plus mean line
pred.est <- laply(test, function(x) colMeans(x$pred))
pred.mean <- data.frame(x = seq(ntime), y = colMeans(pred.est))
pred.melt <- melt(pred.est)
pred.gg <- ggplot(pred.melt, aes(x = Var2, y = value, group = Var1)) 
pred.gg <- pred.gg + geom_line(alpha = 0.1, size = 0.1)
pred.gg <- pred.gg + geom_line(data = pred.mean, 
                             mapping = aes(x = x, y = y, group = NULL), 
                             colour = 'blue')

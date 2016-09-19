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
source('../R/trait_setup.r')
source('../R/sim_from_model.r')
source('../R/advi_post.r')
source('../R/make_plots.r')
source('../data/data_dump/trait_info.data.R')
sight.implied <- sight
source('../data/data_dump/trait_w_gaps.data.R')
sight.obs <- sight


#
theme_set(theme_bw())
cbp.long <- c('#000000', '#004949', '#009292', '#FF7DB6', '#FFB6DB', 
              '#490092', '#006DDB', '#B66DFF', '#6DB6FF', '#B6DBFF', 
              '#920000', '#924900', '#DBD100', '#24FF24', '#FFFF6D')
grab <- laply(seq(5), function(x) seq(from = x, to = length(cbp.long), by = 5))
cbp.long <- cbp.long[t(grab)][-1]
#
nsim <- 100
samp <- sample(1001, nsim)


###########
# advi
post <- list.files('../data/mcmc_out', pattern = 'advi',
                   full.names = TRUE)

# fit w/ implied presences and horseshoe priors
fit1 <- read_one_stan_csv(post[1])
ext1 <- post.advi(fit1)

fit2 <- read_one_stan_csv(post[2])
ext2 <- post.advi(fit2)

# posterior inference plots for advi results
make.plots(ext1 = ext1, name = 'basic', name.name = name.name, 
           group = TRUE)
make.plots(ext1 = ext2, name = 'full', name.name = name.name, 
           group = TRUE, sampling = TRUE)


## full Bayes
#post <- list.files('../data/mcmc_out', pattern = '[0-9]', full.names = TRUE)
#fit <- read_stan_csv(post)
#stan_rhat(fit)
#stan_ess(fit)
#ext <- rstan::extract(fit, permuted = TRUE)
#make.plots(ext1 = ext, name = 'basic_mcmc', name.name = name.name)


# sample, group level, individual level
# send ext$pred through the simulator
ntax <- N
ntime <- T
sim.obs <- sim.imp <- list()
for(ii in seq(nsim)) {
  sim.imp[[ii]] <- model.simulation(N, T, phi = ext1$phi[samp], 
                                    pred = ext1$pred[samp[ii], , ])
  sim.obs[[ii]] <- model.simulation(N, T, phi = ext2$phi[samp],
                                    pred = ext2$pred[samp[ii], , ], 
                                    p = ext2$p[samp[ii], ])
}

meanocc.obs <- mean(rowSums(sight.obs))
meanocc.imp <- mean(rowSums(sight.implied))
meanocc.simobs <- laply(sim.obs, function(x) mean(rowSums(x$y)))
meanocc.simimp <- laply(sim.imp, function(x) mean(rowSums(x$z)))

mos <- data.frame(x = c(meanocc.simobs, meanocc.simimp), 
                  y = c(rep('Full', nsim), rep('Basic', nsim)))
obs <- data.frame(x = c(meanocc.obs, meanocc.imp), 
                  y = c('Full', 'Basic'))
ocplot <- ggplot(mos, aes(x = x)) + geom_histogram()
ocplot <- ocplot + geom_vline(data = obs, mapping = aes(xintercept = x))
ocplot <- ocplot + facet_grid( ~ y)
ocplot <- ocplot + labs(x = 'Mean number of observations per species',
                        y = 'Posterior predictive simulations')
ggsave(filename = '../doc/figure/pred_occ.png', plot = ocplot,
       width = 4, height = 3)

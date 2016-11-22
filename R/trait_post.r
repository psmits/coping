library(plyr)
library(coda)
library(arm)
library(stringr)
library(reshape2)
library(ggplot2)
library(rstan)
library(grid)
library(scales)
library(pROC)
library(Metrics)
source('../R/multiclass_roc.r')
source('../R/trait_setup.r')
source('../R/sim_from_model.r')
source('../R/advi_post.r')
source('../R/make_plots.r')
source('../data/data_dump/trait_w_gaps_revamp.data.R')
sight.obs <- sight


#
theme_set(theme_bw())
theme_update(axis.text = element_text(size = 12),
             axis.title = element_text(size = 15),
             legend.text = element_text(size = 8),
             legend.title = element_text(size = 10),
             legend.key.size = unit(0.75, 'cm'),
             strip.text = element_text(size = 12))

cbp.long <- c('#000000', '#004949', '#009292', '#FF7DB6', '#FFB6DB', 
              '#490092', '#006DDB', '#B66DFF', '#6DB6FF', '#B6DBFF', 
              '#920000', '#924900', '#DBD100', '#24FF24', '#FFFF6D')

grab <- laply(seq(5), function(x) seq(from = x, to = length(cbp.long), by = 5))
cbp.long <- cbp.long[t(grab)][-1]
#
nsim <- 100
samp <- sample(1001, nsim)

rms <- which(inter == names(which(table(inter) == 1)))
inter <- inter[-rms]
break.inter <- str_split(inter, '\\.')
ecotype <- Reduce(rbind, break.inter)

ecotrans <- Reduce(rbind, str_split(levels(as.factor(inter)), '\\.'))

###########
# advi
post <- list.files('../data/mcmc_out', pattern = 'advi',
                   full.names = TRUE)

# fit w/ implied presences and horseshoe priors
fit1 <- read_one_stan_csv(post[1])
ext1 <- post.advi(fit1)

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
sim.obs <- list()
for(ii in seq(nsim)) {
  sim.obs[[ii]] <- model.simulation(ntax, ntime, phi = ext1$phi[samp[ii]], 
                                    pred = ext1$pred[samp[ii], , ],
                                    p = ext1$p[samp[ii], , ])
}

meanocc.obs <- mean(rowSums(sight.obs))
meanocc.simobs <- laply(sim.obs, function(x) mean(rowSums(x$y)))

mos <- data.frame(x = meanocc.simobs, 
                  y = rep('Full', nsim))
obs <- data.frame(x = meanocc.obs, 
                  y = rep('Full', nsim))
ocplot <- ggplot(mos, aes(x = x)) + geom_histogram()
ocplot <- ocplot + geom_vline(data = obs, mapping = aes(xintercept = x), 
                              colour = 'blue', size = 1.5)
ocplot <- ocplot + labs(x = 'Mean obs per species',
                        y = 'Post. pred. simulations')
ggsave(filename = '../doc/figure/pred_occ.png', plot = ocplot,
       width = 4, height = 3)

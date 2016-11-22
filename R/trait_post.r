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


# analysis of model fit
ntax <- N
ntime <- T
post.pred(ext1, ntax, ntime, sight.obs, nsim, samp)


# analysis of the posterior
vis.post(ext1, ecotype, ecotrans, mass, cbp.long)




## full Bayes
#post <- list.files('../data/mcmc_out', pattern = '[0-9]', full.names = TRUE)
#fit <- read_stan_csv(post)
#stan_rhat(fit)
#stan_ess(fit)
#ext <- rstan::extract(fit, permuted = TRUE)
#make.plots(ext1 = ext, name = 'basic_mcmc', name.name = name.name)


# sample, group level, individual level

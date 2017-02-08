library(parallel)
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
source('../R/estimate_div.r')
source('../R/advi_post.r')
source('../R/make_plots.r')
source('../data/data_dump/trait_w_gaps_revamp.data.R')
sight.obs <- sight


#
theme_set(theme_bw())
theme_update(axis.text = element_text(size = 10),
             axis.title = element_text(size = 12),
             legend.text = element_text(size = 5),
             legend.title = element_text(size = 8),
             legend.key.size = unit(0.75, 'cm'),
             strip.text = element_text(size = 9))

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
ecotype <- rbind(ecotype, t(replicate(M - N, c('augment', 'augment'))))

ecotrans <- Reduce(rbind, str_split(levels(as.factor(inter)), '\\.'))
ecotrans <- rbind(ecotrans, c('augment', 'augment'))

ntax <- N
ntime <- T

# observed time bins
time.stop <- unique(occur$true.bin)
b <- range(time.stop)
b <- seq(b[1], b[2], by = 2)
time.start.stop <- as.matrix(cbind(b - 2, b))

############
## advi
post <- list.files('../data/mcmc_out', pattern = 'advi',
                   full.names = TRUE)

# just presence
fit1 <- read_one_stan_csv(post[2])
ext1 <- post.advi(fit1)
# analysis of model fit
post.pred(ext1, ntax = M, ntime = T, sight.obs = sighta, nsim, samp)
# analysis of the posterior
vis.post(ext1, ecotype, ecotrans, mass, 
         cbp.long, time.start.stop, ecoprob = TRUE)


# full birth-death
fit2 <- read_one_stan_csv(post[1])
ext2 <- post.advi(fit2)
# posterior predictive checks
#   need to develop more
post.pred(ext2, ntax = M, ntime = T, sight.obs = sighta, nsim, samp, bd = TRUE)
# visualize posterior estimates
vis.bdpost(ext2, ecotype, ecotrans, mass, cbp.long, ecoprob = TRUE)
# estimate standing diversity given posterior
post.div <- diversity.distribution(sighta, ext2, nsim) # 

source('../R/div_plot.r')  # update this to work as functions, not just source

source('../R/prob_calc.r')  # important posterior probabilities and related



############
## full Bayes
#post <- list.files('../data/mcmc_out', pattern = '[0-9]', full.names = TRUE)
#fit <- read_stan_csv(post)
##stan_rhat(fit)
#ext <- rstan::extract(fit, permuted = TRUE)
##x <- model.simulation(ntax, ntime, 
##                      ext$phi[1], 
##                      ext$pred[1, , ], 
##                      ext$p[1, , ], 
##                      death = TRUE) 
## the issue is that everything needs to occur min 1
##   need data augmentation to make occurs "bigger"
#post.pred(ext, ntax, ntime, sight.obs, nsim, samp)  # posterior pred check
#vis.post(ext, ecotype, ecotrans, mass, cbp.long)    # make some plots

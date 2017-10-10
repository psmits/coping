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
library(ggjoy)
source('../R/borrow_plotcorr.r')
source('../R/multiclass_roc.r')
#source('../R/trait_setup.r')
source('../R/sim_from_model.r')
source('../R/estimate_div.r')
source('../R/advi_post.r')
source('../R/visual_mass.r')
source('../R/make_plots.r')
source('../data/data_dump/trait_w_gaps_NALMA.data.R')
#source('../data/data_dump/trait_w_gaps_revamp.data.R')
sight.obs <- sight
base::load('../data/trait_setup_run.rdata')

#
theme_set(theme_bw())
theme_update(axis.text = element_text(size = 10),
             axis.title = element_text(size = 12),
             legend.text = element_text(size = 5),
             legend.title = element_text(size = 8),
             legend.key.size = unit(0.75, 'cm'),
             strip.text = element_text(size = 8))

cbp.long <- c('#000000', '#004949', '#009292', '#FF7DB6', '#FFB6DB', 
              '#490092', '#006DDB', '#B66DFF', '#6DB6FF', '#B6DBFF', 
              '#920000', '#924900', '#DBD100', '#24FF24', '#FFFF6D')

grab <- laply(seq(5), function(x) seq(from = x, to = length(cbp.long), by = 5))
cbp.long <- cbp.long[t(grab)][-1]
#
nsim <- 100
samp <- sample(1001, nsim)
ecoprob <- TRUE

ecotype <- Reduce(rbind, str_split(inter, '\\.'))
ecotrans <- Reduce(rbind, str_split(levels(as.factor(inter)), '\\.'))
ntax <- N
ntime <- T

# observed time bins
time.stop <- unique(occur$true.bin)
if(bin == '2My') {
  b <- range(time.stop)
  b <- seq(b[1], b[2], by = 2)
  time.start.stop <- as.matrix(cbind(b - 2, b))
} else if(bin == 'NALMA') {
  time.start.stop <- cbind(c(1.8, time.stop[-length(time.stop)]), time.stop)
}



############
## advi
#post <- list.files('../data/mcmc_out', pattern = 'sensible_advi_NALMA',
#                   full.names = TRUE)
post <- list.files('../data/mcmc_out', pattern = 'rwprior_advi_NALMA',
                    full.names = TRUE)

# full birth-death
fit2 <- read_one_stan_csv(post)
ext2 <- post.advi(fit2)

# posterior predictive checks
#   need to develop more
post.pred(ext2, ntax = N, ntime = T, sight.obs = sight, nsim, samp, bd = TRUE)
# visualize posterior estimates
# this is a side-effects function
vis.bdpost(ext2 = ext2, ecotype = ecotype, ecotrans = ecotrans, 
           mass = mass, cbp.long = cbp.long, 
           time.start.stop = time.start.stop, ecoprob = ecoprob, 
           order.cypher = order.cypher)

## estimate standing diversity given posterior
post.div <- diversity.distribution(sight, ext2, nsim) # 

source('../R/div_plot.r')  # update this to work as functions, not just source

source('../R/prob_calc.r')  # important posterior probabilities and related

source('../R/cor_plot.r')  # plot and inspect correlation matrix from b+d



##############
### full Bayes
##post <- list.files('../data/mcmc_out', pattern = '[0-9]', full.names = TRUE)
##fit <- read_stan_csv(post)
###stan_rhat(fit)
##ext <- rstan::extract(fit, permuted = TRUE)
###x <- model.simulation(ntax, ntime, 
###                      ext$phi[1], 
###                      ext$pred[1, , ], 
###                      ext$p[1, , ], 
###                      death = TRUE) 
### the issue is that everything needs to occur min 1
###   need data augmentation to make occurs "bigger"
##post.pred(ext, ntax, ntime, sight.obs, nsim, samp)  # posterior pred check
##vis.post(ext, ecotype, ecotrans, mass, cbp.long)    # make some plots

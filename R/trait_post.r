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
ecotype <- rbind(ecotype, t(replicate(M - N, c('augment', 'augment'))))

ecotrans <- Reduce(rbind, str_split(levels(as.factor(inter)), '\\.'))
ecotrans <- rbind(ecotrans, c('augment', 'augment'))

ntax <- N
ntime <- T


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
vis.post(ext1, ecotype, ecotrans, mass, cbp.long, ecoprob = TRUE)


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

# total diversity not broken by grouping
diversity <- Map(function(x, y) cbind(div = x, time = seq(y)), 
                 llply(post.div, colSums), ncol(sight.obs))
diversity <- Map(function(x, y) cbind(x, sim = y),
                 diversity, seq(nsim))
diversity <- data.frame(Reduce(rbind, diversity))
diversity$div <- log(diversity$div + 1)
divgg <- ggplot(diversity, aes(x = time, y = div, group = sim))
divgg <- divgg + geom_line(alpha = 0.1)

# count number of gains going to t to t+1
#   ask how many 0 -> 1 in interval
gains <- list()
for(jj in seq(nsim)) {
  oo <- matrix(ncol = T - 1, nrow = M)
  for(kk in seq(M)) {
    for(ii in 2:(ncol(oo) + 1)) {
      oo[kk, ii - 1] <- post.div[[jj]][kk, ii - 1] == 0 & 
        post.div[[jj]][kk, ii] == 1
    }
  }
  gains[[jj]] <- oo
}
gains <- lapply(gains, colSums)

# count number of losses going to t to t+1
#   ask how many 1 -> 0 in interval
loss <- list()
for(jj in seq(nsim)) {
  oo <- matrix(ncol = T - 1, nrow = M)
  for(kk in seq(M)) {
    for(ii in 2:(ncol(oo) + 1)) {
      oo[kk, ii - 1] <- post.div[[jj]][kk, ii - 1] == 1 & 
        post.div[[jj]][kk, ii] == 0
    }
  }
  loss [[jj]] <- oo
}
loss <- lapply(loss, colSums)

# i have to figure this out
# i've no idea what's going on here
div <- llply(post.div, colSums)  # need this for rate calculation
growth.rate <- Map(function(x, y, a) (x/a[-(length(a))]) - (y/a[-(length(a))]), 
                   gains, loss, div)
growth.rate <- Map(function(x) data.frame(growth = x, time = seq(length(x))), 
                   growth.rate) 
growth.rate <- Map(function(x, y) cbind(x, sim = y), growth.rate, seq(nsim))
growth.rate <- Reduce(rbind, growth.rate)
growgg <- ggplot(growth.rate, aes(x = time, y = growth, group = sim))
growgg <- growgg + geom_line(alpha = 0.1)


# per capita growth rate (gains per 2 My, loss per 2 My)
grow.sums <- apply(Reduce(rbind, growth.rate), 2, summary)


# break diversity down by ecotype
# plot ideas
#   grid of diversity estimates as lines
#   relative diversity of ecotypes 
#     median estimates
#     allows there to be fractional species
eec <- factor(as.character(interaction(ecotype[, 1], ecotype[, 2])))
et <- as.character(unique(eec))

div.byeco <- list()
for(jj in seq(length(post.div))) {
  byeco <- list()
  for(ii in seq(length(et))) {
    gr <- eec == et[ii]
    byeco[[ii]] <- post.div[[jj]][gr, ]
  }
  div.byeco[[jj]] <- byeco
}

div.eco <- llply(div.byeco, function(x) Reduce(rbind, llply(x, colSums)))
div.eco <- llply(div.eco, function(x) 
                   suppressWarnings(melt(data.frame(x, et))))
div.eco <- Map(function(x, y) cbind(x, sim = y), div.eco, seq(nsim))

div.eco <- llply(div.eco, function(x) {
                   o <- str_extract_all(as.character(x[, 2]), 
                                        '[0-9]+$', simplify = TRUE)
                   x[, 2] <- o
                   x})
div.eco <- Reduce(rbind, div.eco)
div.eco <- cbind(div.eco, str_split(div.eco[, 1], '\\.', simplify = TRUE))
names(div.eco) <- c('et', 'time', 'diversity', 'sim', 'eco_1', 'eco_2')

div.eco <- div.eco[div.eco$eco_1 != 'augment', ]
div.eco$diversity <- log(div.eco$diversity + 1)
degg <- ggplot(div.eco, aes(x = time, y = diversity, group = sim))
degg <- degg + geom_line(alpha = 0.1)
degg <- degg + facet_grid(eco_1 ~ eco_2)

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

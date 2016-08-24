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
source('../data/data_dump/trait_info.data.R')
sight.implied <- sight
source('../data/data_dump/trait_w_gaps.data.R')
sight.obs <- sight


#
theme_set(theme_minimal())
cbp.long <- c('#000000', '#004949', '#009292', '#FF7DB6', '#FFB6DB', 
              '#490092', '#006DDB', '#B66DFF', '#6DB6FF', '#B6DBFF', 
              '#920000', '#924900', '#DBD100', '#24FF24', '#FFFF6D')
grab <- laply(seq(5), function(x) seq(from = x, to = length(cbp.long), by = 5))
cbp.long <- cbp.long[t(grab)][-1]
#
nsim <- 1000


###########
# advi
post <- list.files('../data/mcmc_out', pattern = 'advi',
                   full.names = TRUE)

# fit w/ implied presences and horseshoe priors
fit <- read_one_stan_csv(post[1])
ext <- post.advi(fit)
samp <- sample(1001, nsim)



############
## full Bayes
#post <- list.files('../data/mcmc_out', pattern = '[0-9]', full.names = TRUE)

## horseshoe priors
#fit <- read_stan_csv(post)
#ext <- extract(fit, permuted = TRUE)
#samp <- sample(4000, nsim)




# sample, group level, individual level
# send ext$pred through the simulator
#   ntax <- N
#   ntime <- T
#   pred <- ext$pred[1, , ] #?
#   model.simulation(ntax, ntime, pred)
#sim <- list()
#for(ii in seq(nsim)) {
#  sim[[ii]] <- model.simulation(N, T, ext$pred[samp[ii], , ])
#}
#
#meanocc.obs <- mean(rowSums(sight.obs))
#meanocc.imp <- mean(rowSums(sight.implied))
#meanocc.sim <- laply(sim, function(x) mean(rowSums(x$z)))

es <- list()
for(ii in seq(nrow(u))) {
  byd <- list()
  for(dd in seq(D)) {
    byd[[dd]] <- apply(ext$gamma[samp, , dd], 1, function(x) u[ii, ] %*% x)
  }
  es[[ii]] <- byd
}


# this whole section changes in an interesting way 
# diet and life get combined into one graph
# because i'm looking at interaction
# make it faceted diet ~ life
#   diet columns, life rows
#   drop all non-observed interactions

# plot the parameters estimated in the model
# break the binary factors up
save.cept <- diet <- move <- list()
for(ii in seq(T)) {
  cept <- ext$beta[samp, ii, c(1, 3:22)]

  tt <- es[[ii]][c(1, 3:22)]
  for(jj in ncol(cept)) {
    cept[, jj] <- cept[, jj] - tt[[jj]]
  }

  cept <- invlogit(cept)
  save.cept[[ii]] <- cept
  
  diet[[ii]] <- cept[, c(1:4)]

  move[[ii]] <- cept[, c(1, 5:9)]
}

cept.nam <- name.name[-2]
cept.nam[1] <- 'dietomni:lifescansorial'
cept.nam[2] <- 'dietinsect:lifescansorial'
cept.nam[3] <- 'dietcarni:lifescansorial'
cept.nam[4] <- 'dietherb:lifescansorial'

cept.nam[5] <- 'dietherb:lifearboreal'
cept.nam[6] <- 'dietomni:lifedigitigrade'
cept.nam[7] <- 'dietomni:lifefossorial'
cept.nam[8] <- 'dietomni:lifeunguligrade'
cept.nam[9] <- 'dietomni:lifeplantigrade'


# big graph
save.cept <- save.cept[-c(1, T)]
save.cept <- llply(save.cept, function(x) 
              apply(x, 2, function(y) 
                    quantile(y, c(0.1, 0.25, 0.5, 0.75, 0.9))))
save.cept <- llply(save.cept, function(x) {
              colnames(x) <- cept.nam
              x})
save.cept <- llply(save.cept, t)

save.cept <- Map(function(x, y) data.frame(x, 
                 type = Reduce(rbind, str_split(rownames(save.cept[[1]]), ':'))
                 , time = y), 
                 save.cept, seq(length(save.cept)))
save.cept <- llply(save.cept, function(x) {
              names(x) <- c('low', 'lowmed', 'med', 'highmed', 'high', 
                            'diet', 'move', 'time')
              x})
save.cept <- data.frame(Reduce(rbind, save.cept))

# plot of relative probability of occurrence based on diet
#   controling for expected occurrence probability due to environmental condition
ceptprob <- ggplot(save.cept, aes(x = time, y = med))
ceptprob <- ceptprob + geom_ribbon(aes(ymax = high, ymin = low), 
                                   alpha = 0.25)
ceptprob <- ceptprob + geom_ribbon(aes(ymax = highmed, ymin = lowmed), 
                                   alpha = 0.5)
ceptprob <- ceptprob + facet_grid(diet ~ move)
ceptprob <- ceptprob + scale_fill_manual(values = cbp.long)
ceptprob <- ceptprob + labs(x = 'Time', 
                            y = 'Probability of occurring relative to average')
ggsave(filename = '../doc/figure/cept_occur_prob.png', plot = ceptprob,
       width = 8, height = 6)
#dietprob <- dietprob + geom_ribbon(data = dietdupe, 
#                                 aes(ymax = highmed, ymin = lowmed,
#                                     fill = NULL, group = group), 
#                                 alpha = 0.2, colour = 'grey')



# now for gamma
byindiv <- list()
for(ii in seq(D)) {
  bygroup <- list()
  for(jj in seq(U)) {
    bygroup[[jj]] <- quantile(ext$gamma[samp, jj, ii], 
                              c(0.10, 0.25, 0.5, 0.75, 0.90))
    names(bygroup[[jj]]) <- c('low', 'lowmid', 'med', 'highmed', 'high')
  }
  bygroup <- Reduce(rbind, bygroup)
  byindiv[[ii]] <- data.frame(bygroup, group = seq(U), indiv = ii)
}
melted <- Reduce(rbind, byindiv)
melted$group <- mapvalues(melted$group, unique(melted$group), 
                          c('intercept/phase 1', 'mean temp', 'range temp', 
                            'phase 2', 'phase3'))
melted$group <- factor(melted$group,
                       levels = c('intercept/phase 1', 'mean temp', 'range temp', 
                                  'phase 2', 'phase3'))
melted$indiv <- mapvalues(melted$indiv, unique(melted$indiv), name.name)
melted$indiv <- factor(melted$indiv, levels = name.name)
# iter, indiv, group, value
gamma.plot <- ggplot(melted, aes(x = factor(indiv), y = med))
gamma.plot <- gamma.plot + geom_hline(yintercept = 0, colour = 'grey')
gamma.plot <- gamma.plot + geom_linerange(aes(ymax = highmed, ymin = lowmid), 
                                          size = 2)
gamma.plot <- gamma.plot + geom_pointrange(aes(ymax = high, ymin = low))
gamma.plot <- gamma.plot + facet_grid(group ~ .)
gamma.plot <- gamma.plot + coord_flip()
gamma.plot <- gamma.plot + labs(x = 'Coefficient estimate (log odds scale)',
                                y = 'Individual-level effect')


# shrinkage plots
# for each group level predictor
lamb.shrink <- apply(ext$lambda, 2, function(x) 
                     quantile(x, c(0.1, 0.25, 0.5, 0.75, 0.9)))
# additional shrinkage for individual level predictors
phi.shrink <- list()
for(ii in seq(U)) { # this will need to update
  phi.shrink[[ii]] <- apply(ext$phi[, , ii], 2, function(x)
                            quantile(x, c(0.1, 0.25, 0.5, 0.75, 0.9)))
}
#lamb.shrink
#phi.shrink
# these are estimates of variance
# lamb * phi is the variance of around gamma_u,d

# group level predictor shrinkage
colnames(lamb.shrink) <- name.name
lamb.shrink <- data.frame(t(lamb.shrink))
names(lamb.shrink) <- c('low', 'lowmed', 'med', 'hghmed', 'hgh')
lamb.shrink$var <- rownames(lamb.shrink)
lamb.shrink$var <- factor(lamb.shrink$var, levels = name.name)
group.shrink <- ggplot(lamb.shrink, aes(x = var, y = med))
group.shrink <- group.shrink + geom_point(size = 5)
group.shrink <- group.shrink + geom_linerange(mapping = aes(ymax = hghmed, 
                                                            ymin = lowmed),
                                               size = 1.75)
group.shrink <- group.shrink + geom_linerange(mapping = aes(ymax = hgh, 
                                                             ymin = low),
                                              size = 0.5)

# phi

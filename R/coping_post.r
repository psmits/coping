library(sn)
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

# bring in the data
source('../R/coping_setup.r')
dat.sum <- laply(split(data$trait, data$year), 
                 function(x) quantile(x, c(0.25, 0.5, 0.75)))
dat.mm <- laply(split(data$trait, data$year), mean)
dat.ss <- laply(split(data$trait, data$year), sd)

# mcmc results
post <- list.files('../data/mcmc_out', 
                   pattern = 'coping', 
                   full.names = TRUE)
fit <- read_stan_csv(post)
ext <- extract(fit, permuted = TRUE)

# posterior predictive replications
pp <- replicate(nsim, sim.series(data = data, ext), simplify = FALSE)


expect <- llply(pp, function(x) laply(x, mean))
expect <- Reduce(cbind, expect)

# basic posterior predictive checks
post.quant <- llply(pp, function(y) 
                    laply(y, function(x) quantile(x, c(0.25, 0.5, 0.75))))
post.sd <- Reduce(cbind, llply(pp, function(y) laply(y, sd)))

qex <- qsd <- q25 <- q50 <- q75 <- c()
for(ii in seq(data$T)) {
  qex[ii] <- sum(expect[ii, ] > dat.mm[ii]) / nsim
  qsd[ii] <- sum(post.sd[ii, ] > dat.ss[ii]) / nsim
  q25[ii] <- sum(laply(post.quant, 
                       function(x) x[ii, 1]) > dat.sum[ii, 1]) / nsim
  q50[ii] <- sum(laply(post.quant, 
                       function(x) x[ii, 2]) > dat.sum[ii, 2]) / nsim
  q75[ii] <- sum(laply(post.quant, 
                       function(x) x[ii, 3]) > dat.sum[ii, 3]) / nsim
}

df <- data.frame(year = data$year, trait = data$trait)

expp <- list()
for(ii in seq(nsim)) {
  expp[[ii]] <- data.frame(sim = ii, year = seq(nrow(expect)),
                           val = expect[, sample(ncol(expect), 1)])
}
expp <- Reduce(rbind, expp)

df2 <- data.frame(year = seq(nrow(expect)), 
                  expect.o = dat.mm)

# data plus observed mean and 1000 estimates of expectation
data.view <- ggplot(df, aes(x = year, y = trait))
data.view <- data.view + geom_point(alpha = (nsim / 10 ^ (nchar(nsim) + 1)) * 5)
data.view <- data.view + geom_line(data = expp, 
                                   mapping = aes(y = val, group = sim), 
                                   alpha = (nsim / 10 ^ (nchar(nsim) + 1)) * 5,
                                   colour = 'blue')
data.view <- data.view + geom_line(data = df2, mapping = aes(y = expect.o), 
                                   colour = 'black', size = 1.5)
data.view <- data.view + labs(x = 'Time', y = 'log body mass')
ggsave(filename = '../doc/figure/data_expect.png', plot = data.view,
       width = 10, height = 5)


# skewness through time
skpp <- list()
for(ii in seq(nsim)) {
  skpp[[ii]] <- data.frame(sim = ii, year = seq(ncol(ext$skew)), 
                           val = ext$skew[sample(nrow(ext$skew), 1), ])
}
skpp <- Reduce(rbind, skpp)
skmm <- data.frame(year = seq(ncol(ext$skew)), val = colMeans(ext$skew))

skew.plot <- ggplot(skpp, aes(x = year, y = val, group = sim))
skew.plot <- skew.plot + geom_line(alpha = (nsim / 10 ^ (nchar(nsim) + 1)) * 5, 
                                   colour = 'blue')
skew.plot <- skew.plot + geom_line(data = skmm, 
                                   mapping = aes(x = year, 
                                                 y = val, 
                                                 group = NULL),
                                   colour = 'black', size = 2)
skew.plot <- skew.plot + labs(x = 'Time', y = 'alpha (skewness)')
ggsave(filename = '../doc/figure/skew_series.png', plot = skew.plot,
       width = 10, height = 5)



## plot what that type of skewness looks like
#xi = mean(ext$loc)
#omega = mean(exp(ext$scale_mu))
#alpha = mean(ext$skew_mu)
#
#eep <- data.frame(sim = seq(ncol(expect)), val = colMeans(expect))
#
#x <- data.frame(x = seq(from = -2.3, to = 15, by = 0.01))
#theo.skew <- ggplot(x, aes(x = x))
#theo.skew <- theo.skew + stat_function(fun = dsn, 
#                                       args = list(xi = xi, 
#                                                   omega = omega, 
#                                                   alpha = alpha))
#theo.skew <- theo.skew + geom_vline(data = eep, 
#                                    aes(xintercept = val), 
#                                    colour = 'blue', alpha = 0.05)
#theo.skew <- theo.skew + geom_vline(xintercept = mean(data$trait), 
#                                    colour = 'black')
#theo.skew <- theo.skew + labs(y = 'Density', x = 'log body mass')
#ggsave(filename = '../doc/figure/skew_theo.png', plot = theo.skew,
#       width = 7, height = 5)

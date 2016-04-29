library(arm)
library(plyr)
library(reshape2)
library(ggplot2)
library(grid)

ntax <- 1000
ntime <- 32


y <- z <- pred <- matrix(NA, nrow = ntax, ncol = ntime)
p <- c()

p.mu <- rnorm(1, 0, 1)
p.sigma <- abs(dt(1, df = 1))
p.eff <- rnorm(ntime, 0, 1)

inter.mu <- rnorm(1, -2, 2)
inter.sigma <- abs(dt(1, df = 1))
inter.eff <- rnorm(ntime, 0, 1)


for(nn in 1:ntax) {
  p[1] <- invlogit(p.mu + p.sigma * p.eff[1])
  pred[nn, 1] <- invlogit(inter.mu + inter.sigma * inter.eff[1])
  z[nn, 1] <- rbinom(1, 1, prob = pred[nn, 1])
  y[nn, 1] <- rbinom(1, 1, prob = z[nn, 1] * p[1])
  for(tt in 2:ntime) {
    p[tt] <- invlogit(p.mu + p.sigma * p.eff[tt])
    pred[nn, tt] <- invlogit(inter.mu + inter.sigma * inter.eff[tt])

    z[nn, tt] <- rbinom(1, 1, prob = z[nn, tt - 1] * pred[nn, tt] + 
                        (prod(1 - z[nn, 1:tt - 1]) * pred[nn, tt]))
    y[nn, tt] <- rbinom(1, 1, prob = z[nn, tt] * p[tt])
  }
}

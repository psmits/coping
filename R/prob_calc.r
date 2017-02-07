# need to calculate a bunch of posterior probabilites
#   group-level effects
#     P(eff of temp < 0)
#     P(phase a > phase b)
#     P(rate_t > e[rate])
#   diversity peaks/valleys vs average diversity

# probability that group-level on origin is greater than 0
# order is phase 3, temp mean, temp range, phase 2, phase 1
od <- c(5, 4, 1, 2, 3)
out.occur <- out.origin <- out.surv <- list()
for(ii in seq(nrow(ecotrans))) {
  out.occur[[ii]] <- apply(ext1$gamma[, 2:3, ii], 2, function(x) 
                           sum(x > 0) / length(x))
  out.origin[[ii]] <- apply(ext2$o_gamma[, 2:3, ii], 2, function(x) 
                            sum(x > 0) / length(x))
  out.surv[[ii]] <- apply(ext2$s_gamma[, 2:3, ii], 2, function(x) 
                          sum(x > 0) / length(x))
}
out.occur <- Reduce(rbind, out.occur)

out.origin <- Reduce(rbind, out.origin)
out.surv <- Reduce(rbind, out.surv)


# calculate differences amoungst plant phases
plant.oc <- plant.or <- plant.su <- list()
for(ii in seq(nrow(ecotrans))) {
  tt <- ext2$gamma[, c(5, 4, 1), ii]
  plantp <- c()
  plantp[1] <- sum(tt[, 1] - tt[, 2] > 0) / nrow(tt)
  plantp[2] <- sum(tt[, 2] - tt[, 3] > 0) / nrow(tt)
  plantp[3] <- sum(tt[, 1] - tt[, 3] > 0) / nrow(tt)
  plant.oc[[ii]] <- plantp

  tt <- ext2$o_gamma[, c(5, 4, 1), ii]
  plantp <- c()
  plantp[1] <- sum(tt[, 1] - tt[, 2] > 0) / nrow(tt)
  plantp[2] <- sum(tt[, 2] - tt[, 3] > 0) / nrow(tt)
  plantp[3] <- sum(tt[, 1] - tt[, 3] > 0) / nrow(tt)
  plant.or[[ii]] <- plantp

  tt <- ext2$s_gamma[, c(5, 4, 1), ii]
  plantp <- c()
  plantp[1] <- sum(tt[, 1] - tt[, 2] > 0) / nrow(tt)
  plantp[2] <- sum(tt[, 2] - tt[, 3] > 0) / nrow(tt)
  plantp[3] <- sum(tt[, 1] - tt[, 3] > 0) / nrow(tt)
  plant.su[[ii]] <- plantp
}
plant.oc <- Reduce(rbind, plant.oc)
plant.or <- Reduce(rbind, plant.or)
plant.su <- Reduce(rbind, plant.su)


# diversity and diversification
# calculate average diversity at every time point
div.calc <- Reduce(rbind, llply(post.div, colSums))
avg.div <- rowMeans(div.calc)

pdiv <- c()
for(ii in seq(ncol(div.calc))) {
  pdiv[ii] <- sum(div.calc[, ii] > avg.div) / nrow(div.calc)
}

# diversification rate
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

div <- llply(post.div, colSums)  # need this for rate calculation
growth.rate <- Map(function(x, y, a) (x/a[-(length(a))]) - (y/a[-(length(a))]), 
                   gains, loss, div)
avg.grow <- laply(growth.rate, mean)
growth.rate <- Reduce(rbind, growth.rate)

# can then calculate p that t > average
prate <- c()
for(ii in seq(ncol(growth.rate))) {
  prate[ii] <- sum(growth.rate[, ii] > avg.grow) / nrow(growth.rate)
}

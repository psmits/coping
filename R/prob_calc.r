library(xtable)

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

# make fancy tables
ecotrans[, 1] <- mapvalues(ecotrans[, 1], 
                           unique(ecotrans[, 1]),
                           c('carnivore', 'herbivore', 
                             'insectivore', 'omnivore'))
rown <- apply(ecotrans, 1, function(x) Reduce(paste, rev(x)))

rownames(out.occur) <- rownames(out.origin) <- rownames(out.surv) <- rown
colnames(out.occur) <- colnames(out.origin) <- colnames(out.surv) <- 
  c('P(gamma_{temp mean} > 0)', 'P(gamma_{temp range} > 0)')

# make/print tables
out.occur.tab <- xtable(out.occur, label = 'tab:occur_temp', digits = 3)
print.xtable(x = out.occur.tab, file = '../doc/occur_temp_raw.tex')

out.origin.tab <- xtable(out.origin, label = 'tab:origin_temp', digits = 3)
print.xtable(x = out.origin.tab, file = '../doc/origin_temp_raw.tex')

out.surv.tab <- xtable(out.surv, label = 'tab:surv_temp', digits = 3)
print.xtable(x = out.surv.tab, file = '../doc/surv_temp_raw.tex')



# calculate differences amoungst plant phases
plant.oc <- plant.or <- plant.su <- list()
for(ii in seq(nrow(ecotrans))) {
  tt <- ext1$gamma[, c(5, 4, 1), ii]
  tt[, 2] <- tt[, 2] + tt[, 1]
  tt[, 3] <- tt[, 3] + tt[, 1]
  plantp <- c()
  plantp[1] <- sum(tt[, 1] - tt[, 2] > 0) / nrow(tt)
  plantp[2] <- sum(tt[, 2] - tt[, 3] > 0) / nrow(tt)
  plantp[3] <- sum(tt[, 1] - tt[, 3] > 0) / nrow(tt)
  plant.oc[[ii]] <- plantp

  tt <- ext2$o_gamma[, c(5, 4, 1), ii]
  tt[, 2] <- tt[, 2] + tt[, 1]
  tt[, 3] <- tt[, 3] + tt[, 1]
  plantp <- c()
  plantp[1] <- sum(tt[, 1] - tt[, 2] > 0) / nrow(tt)
  plantp[2] <- sum(tt[, 2] - tt[, 3] > 0) / nrow(tt)
  plantp[3] <- sum(tt[, 1] - tt[, 3] > 0) / nrow(tt)
  plant.or[[ii]] <- plantp

  tt <- ext2$s_gamma[, c(5, 4, 1), ii]
  tt[, 2] <- tt[, 2] + tt[, 1]
  tt[, 3] <- tt[, 3] + tt[, 1]
  plantp <- c()
  plantp[1] <- sum(tt[, 1] - tt[, 2] > 0) / nrow(tt)
  plantp[2] <- sum(tt[, 2] - tt[, 3] > 0) / nrow(tt)
  plantp[3] <- sum(tt[, 1] - tt[, 3] > 0) / nrow(tt)
  plant.su[[ii]] <- plantp
}
plant.oc <- Reduce(rbind, plant.oc)
plant.or <- Reduce(rbind, plant.or)
plant.su <- Reduce(rbind, plant.su)

rownames(plant.oc) <- rownames(plant.or) <- rownames(plant.su) <- rown
colnames(plant.oc) <- colnames(plant.or) <- colnames(plant.su) <- 
  c('P(Phase 1 > Phase 2)', 'P(Phase 2 > Phase 3)', 'P(Phase 1 > Phase 3)')

plant.oc.tab <- xtable(plant.oc, label = 'tab:occur_plant', digits = 3)
print.xtable(x = plant.oc.tab, file = '../doc/occur_plant_raw.tex')

plant.or.tab <- xtable(plant.or, label = 'tab:origin_plant', digits = 3)
print.xtable(x = plant.oc.tab, file = '../doc/origin_plant_raw.tex')

plant.su.tab <- xtable(plant.su, label = 'tab:surv_plant', digits = 3)
print.xtable(x = plant.su.tab, file = '../doc/surv_plant_raw.tex')



# diversity and diversification
# calculate average diversity at every time point
div.calc <- Reduce(rbind, llply(post.div, colSums))
avg.div <- rowMeans(div.calc)

pdiv <- c()
for(ii in seq(ncol(div.calc))) {
  pdiv[ii] <- sum(div.calc[, ii] > avg.div) / nrow(div.calc)
}

div.peaks <- data.frame(time = rev(sort(unique(diversity$time))), prob = pdiv)
names(div.peaks) <- c('Time (Mya)', 'P(N^{stand}_{t} > bar{N^{stand}})')
div.peaks.tab <- xtable(div.peaks, label = 'tab:div_peak')
print.xtable(x = div.peaks.tab, file = '../doc/div_peak_raw.tex')


# diversification rate
# count number of gains going to t to t+1
#   ask how many 0 -> 1 in interval
gains <- list()
for(jj in seq(nsim)) {
  oo <- matrix(ncol = T - 1, nrow = N)
  for(kk in seq(N)) {
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
  oo <- matrix(ncol = T - 1, nrow = N)
  for(kk in seq(N)) {
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

rat.peaks <- data.frame(time = rev(sort(unique(diversity$time))[-31]), prob = prate)
names(rat.peaks) <- c('Time (Mya)', 'P(D^{rate}_{t} > bar{D^{rate}})')
rat.peaks.tab <- xtable(rat.peaks, label = 'tab:rate_peak')
print.xtable(x = rat.peaks.tab, file = '../doc/rate_peak_raw.tex')

library(xtable)

# need to calculate a bunch of posterior probabilites
#   group-level effects
#     P(eff of temp < 0)
#     P(phase a > phase b)
#     P(rate_t > e[rate])
#   diversity peaks/valleys vs average diversity

# probability that group-level on origin is greater than 0
# order is phase (Mi.Pl), Eo.Mi, Pa.Eo, temp
out.origin <- out.surv <- list()
for(ii in seq(nrow(ecotrans))) {
  out.origin[[ii]] <- sum(ext2$o_gamma[, 3, ii] > 0) / 
    length(ext2$o_gamma[, 3, ii])
  out.surv[[ii]] <- sum(ext2$s_gamma[, 3, ii] > 0) / 
    length(ext2$s_gamma[, 3, ii])
}
out.origin <- Reduce(rbind, out.origin)
out.surv <- Reduce(rbind, out.surv)

# make fancy tables
ecotrans[, 1] <- mapvalues(ecotrans[, 1], 
                           unique(ecotrans[, 1]),
                           c('carnivore', 'herbivore', 
                             'insectivore', 'omnivore'))
rown <- apply(ecotrans, 1, function(x) Reduce(paste, rev(x)))

rownames(out.origin) <- rownames(out.surv) <- rown
colnames(out.origin) <- colnames(out.surv) <- 
  c('P(gamma_{temp mean} > 0)')

# make/print tables
out.origin.tab <- xtable(out.origin, label = 'tab:origin_temp', digits = 3)
print.xtable(x = out.origin.tab, file = '../doc/origin_temp_raw.tex')

out.surv.tab <- xtable(out.surv, label = 'tab:surv_temp', digits = 3)
print.xtable(x = out.surv.tab, file = '../doc/surv_temp_raw.tex')



# calculate differences amoungst plant phases
plant.or <- plant.su <- list()
for(ii in seq(nrow(ecotrans))) {
  tt <- ext2$o_gamma[, 1:2, ii]
  plantp <- c()
  plantp[1] <- sum(tt[, 1] > 0) / nrow(tt)
  plantp[2] <- sum(tt[, 2] > 0) / nrow(tt)
  plantp[3] <- sum(tt[, 1] - tt[, 2] > 0) / nrow(tt)
  plant.or[[ii]] <- plantp

  tt <- ext2$s_gamma[, 1:2, ii]
  plantp <- c()
  plantp[1] <- sum(tt[, 1] > 0) / nrow(tt)
  plantp[2] <- sum(tt[, 2] > 0) / nrow(tt)
  plantp[3] <- sum(tt[, 1] - tt[, 2] > 0) / nrow(tt)
  plant.su[[ii]] <- plantp
}
plant.or <- Reduce(rbind, plant.or)
plant.su <- Reduce(rbind, plant.su)

rownames(plant.or) <- rownames(plant.su) <- rown
colnames(plant.or) <- colnames(plant.su) <- 
  c('P(Eo.Mi > 0)', 'P(Pa.Eo > 0)', 'P(Eo.Mi > Pa.Eo)')

plant.or.tab <- xtable(plant.or, label = 'tab:origin_plant', digits = 3)
print.xtable(x = plant.or.tab, file = '../doc/origin_plant_raw.tex')

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

tm <- match(diversity$time, nalma$ma)
diversity$nalma <- nalma$interval[tm]
diversity$nalma <- factor(diversity$nalma, levels = nalma$interval)

div.peaks <- data.frame(time = unique(diversity$nalma), prob = pdiv)
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

dl <- length(unique(diversity$time))
rat.peaks <- data.frame(time = unique(diversity$nalma)[-dl], prob = prate)
names(rat.peaks) <- c('Time (Mya)', 'P(D^{rate}_{t} > bar{D^{rate}})')
rat.peaks.tab <- xtable(rat.peaks, label = 'tab:rate_peak')
print.xtable(x = rat.peaks.tab, file = '../doc/rate_peak_raw.tex')

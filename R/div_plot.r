# total diversity not broken by grouping
diversity <- Map(function(x, y) cbind(div = x, time = seq(y)), 
                 llply(post.div, colSums), ncol(sight.obs))
diversity <- Map(function(x, y) cbind(x, sim = y),
                 diversity, seq(nsim))
diversity <- data.frame(Reduce(rbind, diversity))
diversity$div <- log(diversity$div + 1)

diversity$time <- mapvalues(diversity$time, 
                            from = unique(diversity$time), 
                            to = rev(time.start.stop[, 2]))


div.calc <- Reduce(rbind, llply(post.div, colSums))
avg.div <- log(rowMeans(div.calc) + 1)
avg.div <- quantile(avg.div, c(0.1, 0.5, 0.9))
names(avg.div) <- c('low', 'mid', 'high')

avg.div <- data.frame(rbind(c(avg.div, time = min(diversity$time)), 
                            c(avg.div, time = max(diversity$time))))

divgg <- ggplot(diversity, aes(x = time, y = div, group = sim))
divgg <- divgg + geom_line(data = avg.div,
                           mapping = aes(x = time, y = mid, group = NULL),
                           colour = 'blue',
                           alpha = 0.5)
divgg <- divgg + geom_ribbon(data = avg.div,
                             mapping = aes(x = time, 
                                           ymin = low, 
                                           y = NULL,
                                           ymax = high, 
                                           group = NULL),
                             fill = 'blue',
                             alpha = 0.25)
divgg <- divgg + geom_line(alpha = 0.1)
divgg <- divgg + labs(x = 'Time (My)', 
                      y = expression(log~N^{textstyle('stand')}))
ggsave(filename = '../doc/figure/log_diversity.png', plot = divgg,
       width = 6, height = 4)

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
avg.grow <- laply(growth.rate, mean)
growth.rate <- Map(function(x) data.frame(growth = x, time = seq(length(x))), 
                   growth.rate) 
growth.rate <- Map(function(x, y) cbind(x, sim = y), growth.rate, seq(nsim))
growth.rate <- Reduce(rbind, growth.rate)

growth.rate$time <- mapvalues(growth.rate$time, 
                              from = unique(growth.rate$time), 
                              to = rev(time.start.stop[-T, 2]))

avg.grow <- quantile(avg.grow, c(0.1, 0.5, 0.9))
names(avg.grow) <- c('low', 'mid', 'high')

avg.grow <- data.frame(rbind(c(avg.grow, time = min(growth.rate$time)),
                             c(avg.grow, time = max(growth.rate$time))))

growgg <- ggplot(growth.rate, aes(x = time, y = growth, group = sim))
growgg <- growgg + geom_line(data = avg.grow,
                           mapping = aes(x = time, y = mid, group = NULL),
                           colour = 'blue',
                           alpha = 0.5)
growgg <- growgg + geom_ribbon(data = avg.grow,
                             mapping = aes(x = time, 
                                           ymin = low, 
                                           y = NULL,
                                           ymax = high, 
                                           group = NULL),
                             fill = 'blue',
                             alpha = 0.25)
growgg <- growgg + geom_line(alpha = 0.1)
growgg <- growgg + labs(x = 'Time (My)', 
                        y = 'diversification rate (species/species/time unit')
ggsave(filename = '../doc/figure/div_rate.png', plot = growgg,
       width = 6, height = 4)


# per capita growth rate (gains per 2 My, loss per 2 My)
#grow.sums <- apply(Reduce(rbind, growth.rate), 2, summary)
birth.rate <- Map(function(x, a) x / a[-length(a)], gains, div)
avg.birth <- laply(birth.rate, mean)
birth.rate <- lapply(birth.rate, function(x) 
                     cbind(rate = x, time = seq(length(x))))
birth.rate <- Map(function(x, y) data.frame(x, sim = y), 
                  x = birth.rate, y = seq(length(birth.rate)))
birth.rate <- Reduce(rbind, birth.rate)

birth.rate$time <- mapvalues(birth.rate$time, 
                             from = unique(birth.rate$time), 
                             to = rev(time.start.stop[-T, 2]))

avg.birth <- quantile(avg.birth, c(0.1, 0.5, 0.9))
names(avg.birth) <- c('low', 'mid', 'high')

avg.birth <- data.frame(rbind(c(avg.birth, time = min(birth.rate$time)),
                             c(avg.birth, time = max(birth.rate$time))))

birthgg <- ggplot(birth.rate, aes(x = time, y = rate, group = sim))
birthgg <- birthgg + geom_line(data = avg.birth,
                               mapping = aes(x = time, y = mid, group = NULL),
                               colour = 'blue',
                               alpha = 0.5)
birthgg <- birthgg + geom_ribbon(data = avg.birth,
                                 mapping = aes(x = time, 
                                               ymin = low, 
                                               y = NULL,
                                               ymax = high, 
                                               group = NULL),
                                 fill = 'blue',
                                 alpha = 0.25)
birthgg <- birthgg + geom_line(alpha = 0.1)
birthgg <- birthgg + labs(x = 'Time (My)',
                          y = 'origination rate (originations/species/time unit')
ggsave(filename = '../doc/figure/orig_rate.png', plot = birthgg,
       width = 6, height = 4)



death.rate <- Map(function(x, a) x / a[-length(a)], loss, div)
avg.death <- laply(death.rate, mean)
death.rate <- lapply(death.rate, function(x) 
                     cbind(rate = x, time = seq(length(x))))
death.rate <- Map(function(x, y) data.frame(x, sim = y), 
                  x = death.rate, y = seq(length(death.rate)))
death.rate <- Reduce(rbind, death.rate)

death.rate$time <- mapvalues(death.rate$time, 
                             from = unique(death.rate$time), 
                             to = rev(time.start.stop[-T, 2]))

avg.death <- quantile(avg.death, c(0.1, 0.5, 0.9))
names(avg.death) <- c('low', 'mid', 'high')

avg.death <- data.frame(rbind(c(avg.death, time = min(death.rate$time)),
                             c(avg.death, time = max(death.rate$time))))

deathgg <- ggplot(death.rate, aes(x = time, y = rate, group = sim))
deathgg <- deathgg + geom_line(data = avg.death,
                           mapping = aes(x = time, y = mid, group = NULL),
                           colour = 'blue',
                           alpha = 0.5)
deathgg <- deathgg + geom_ribbon(data = avg.death,
                             mapping = aes(x = time, 
                                           ymin = low, 
                                           y = NULL,
                                           ymax = high, 
                                           group = NULL),
                             fill = 'blue',
                             alpha = 0.25)
deathgg <- deathgg + geom_line(alpha = 0.1)
deathgg <- deathgg + labs(x = 'Time (My)',
                          y = 'extinction rate (extinctions/species/time unit')
ggsave(filename = '../doc/figure/death_rate.png', plot = deathgg,
       width = 6, height = 4)



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
div.eco$time <- as.numeric(div.eco$time)

div.eco$time <- mapvalues(div.eco$time, 
                          from = unique(div.eco$time), 
                          to = rev(time.start.stop[, 2]))

degg <- ggplot(div.eco, aes(x = time, y = diversity, group = sim))
degg <- degg + geom_line(alpha = 0.1)
degg <- degg + facet_grid(eco_1 ~ eco_2)
degg <- degg + labs(x = 'Time (My)', 
                    y = expression(log~N^{textstyle('stand')}))
ggsave(filename = '../doc/figure/ecotype_diversity.png', plot = degg,
       width = 6, height = 4)

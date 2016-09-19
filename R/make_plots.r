#' Make the plots for the coping project
#'
#' @param ext1 output from rstan::extract of model fit
#' @param name string associated with outputs
#' @param name.name vector of regression coefficient names
#' @param group logical are there any group level covariates?
#' @param sampling logical is sampling parameterized?
#' @param nsim integer how many random pulls from the posterior?
make.plots <- function(ext1, 
                       name = 'basic', 
                       name.name, 
                       group = TRUE, 
                       sampling = FALSE, 
                       nsim = 1000) {
  samp <- sample(1001, nsim)

  if(sampling) {
    ss <- t(apply(ext2$p, 2, function(x) 
                  quantile(x, c(0.1, 0.25, 0.5, 0.75, 0.9))))
    ss <- data.frame(cbind(seq(nrow(ss)), ss))
    names(ss) <- c('time', 'low', 'lowmed', 'med', 'highmed', 'high')
    fosplot <- ggplot(ss, aes(x = time, y = med))
    fosplot <- fosplot + geom_ribbon(aes(ymax = high, ymin = low), 
                                     alpha = 0.2)
    fosplot <- fosplot + geom_ribbon(aes(ymax = highmed, ymin = lowmed), 
                                     alpha = 0.4)
    fosplot <- fosplot + geom_line(size = 1.5)
    fosplot <- fosplot + labs(x = 'Time (My)', y = 'Sampling probability')
    ggsave(filename = paste0('../doc/figure/samp_prob_', name, '.png'),
           plot = fosplot, width = 4, height = 3)
  }


  # set of the predictors
  es <- list()
  for(ii in seq(nrow(u))) {
    byd <- list()
    for(dd in seq(D)) {
      if(group) {
        byd[[dd]] <- apply(ext1$gamma[samp, , dd], 1, function(x) u[ii, ] %*% x)
      } else {
        byd[[dd]] <- ext1$gamma[samp, dd]
      }
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
  for(ii in seq(T - 1)) {
    cept <- ext1$beta[samp, ii, c(1, 3:22)]

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
  save.cept$diet <- mapvalues(save.cept$diet, unique(save.cept$diet), 
                              c('omnivore', 'insectivore', 'carnivore', 'herbivore'))
  save.cept$move <- str_replace(save.cept$move, 'life', '')

  # plot of relative probability of occurrence based on diet
  #   controling for expected occurrence probability due to environmental condition
  ceptprob <- ggplot(save.cept, aes(x = time, y = med))
  ceptprob <- ceptprob + geom_ribbon(aes(ymax = high, ymin = low), 
                                     alpha = 0.2)
  ceptprob <- ceptprob + geom_ribbon(aes(ymax = highmed, ymin = lowmed), 
                                     alpha = 0.4)
  ceptprob <- ceptprob + geom_line(size = 1.5)
  ceptprob <- ceptprob + geom_hline(yintercept = 0.5)
  ceptprob <- ceptprob + facet_grid(diet ~ move)
  ceptprob <- ceptprob + scale_fill_manual(values = cbp.long)
  ceptprob <- ceptprob + labs(x = 'Time', 
                              y = 'Probability of occurring relative to average')
  ggsave(filename = paste0('../doc/figure/cept_occur_prob_', name, '.png'), 
         plot = ceptprob, width = 8, height = 6)
  #dietprob <- dietprob + geom_ribbon(data = dietdupe, 
  #                                 aes(ymax = highmed, ymin = lowmed,
  #                                     fill = NULL, group = group), 
  #                                 alpha = 0.2, colour = 'grey')

  masseff <- ext1$beta[, 2, ]
  masseff <- adply(masseff, 2, function(y) 
                   quantile(y, c(0.1, 0.25, 0.5, 0.75, 0.9)))
  names(masseff) <- c('bin', 'low', 'lowmed', 'med', 'highmed', 'high')
  masseff$bin <- as.numeric(as.character(masseff$bin))
  massplot <- ggplot(masseff, aes(x = bin, y = med))
  massplot <- massplot + geom_ribbon(aes(ymax = high, ymin = low),
                                     alpha = 0.2)
  massplot <- massplot + geom_ribbon(aes(ymax = highmed, ymin = lowmed),
                                     alpha = 0.4)
  massplot <- massplot + geom_line(size = 1.5)
  massplot <- massplot + labs(x = 'Time', 
                              y = 'log-odds of occurrence / \nstandard deviation of mass (g)')
  ggsave(filename = paste0('../doc/figure/mass_eff_', name, '.png'),
         plot = massplot, width = 4, height = 3)


  # variance of effect 
  eff.var <- list()
  for(ii in seq(22)) {
    eff.var[[ii]] <- apply(ext1$beta[, , ii], 1, sd)
  }
  name.name[-2] <- cept.nam
  names(eff.var) <- name.name
  eff.var <- melt(eff.var)
  eff.var <- eff.var[eff.var$L1 != 'mass', ]

  eff.var <- data.frame(eff.var, Reduce(rbind, str_split(eff.var$L1, ':')))
  eff.var$X1 <- str_replace(eff.var$X1, 'diet', '')
  eff.var$X1 <- mapvalues(eff.var$X1, unique(eff.var$X1), 
                          c('omnivore', 'insectivore', 'carnivore', 'herbivore'))
  eff.var$X2 <- str_replace(eff.var$X2, 'life', '')

  varplot <- ggplot(eff.var, aes(x = value, y = ..density..))
  varplot <- varplot + geom_histogram()
  varplot <- varplot + facet_grid(X1 ~ X2)
  varplot <- varplot + labs(x = 'Stdev Prob. of ecotype given occurrence', y = 'Probability density')
  ggsave(filename = paste0('../doc/figure/sd_occur_prob_', name, '.png'),
         plot = varplot, width = 8, height = 6)




  if(group) {
    # now for gamma/group level covariates
    byindiv <- list()
    for(ii in seq(D)) {
      bygroup <- list()
      for(jj in seq(U)) {
        bygroup[[jj]] <- quantile(ext1$gamma[samp, jj, ii], 
                                  c(0.10, 0.25, 0.5, 0.75, 0.90))
        names(bygroup[[jj]]) <- c('low', 'lowmid', 'med', 'highmed', 'high')
      }
      bygroup <- Reduce(rbind, bygroup)
      byindiv[[ii]] <- data.frame(bygroup, group = seq(U), indiv = ii)
    }
    melted <- Reduce(rbind, byindiv)
    melted$group <- mapvalues(melted$group, unique(melted$group), 
                              c('phase 1', 'mean temp', 'range temp', 
                                'phase 2', 'phase 3'))
    melted$group <- factor(melted$group,
                           levels = c('phase 1', 'phase 2', 'phase 3', 
                                      'mean temp', 'range temp'))
    melted$indiv <- mapvalues(melted$indiv, unique(melted$indiv), name.name)
    melted$indiv <- factor(melted$indiv, levels = sort(name.name))
    # iter, indiv, group, value
    gamma.plot <- ggplot(melted, aes(x = factor(indiv), y = med))
    gamma.plot <- gamma.plot + geom_hline(yintercept = 0, colour = 'grey')
    gamma.plot <- gamma.plot + geom_linerange(aes(ymax = highmed, ymin = lowmid), 
                                              size = 2)
    gamma.plot <- gamma.plot + geom_pointrange(aes(ymax = high, ymin = low))
    gamma.plot <- gamma.plot + facet_grid(~ group)
    gamma.plot <- gamma.plot + coord_flip()
    gamma.plot <- gamma.plot + labs(x = 'Coefficient estimate (log odds scale)',
                                    y = 'Individual-level effect')
    ggsave(filename = paste0('../doc/figure/gamma_est_', name, '.png'),
           plot = gamma.plot, width = 8, height = 6)
  }
}



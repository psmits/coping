#' Visualize key aspects of the posterior distribution for the presence model
#'
#' @param ext1 extract from stanfit object
#' @param ecotype translate ecotype coding
#' @param ecotrans translate and break ecotype code
#' @param mass species mass estimates
#' @param cbp.long color pallette information for ggplot
#' @param ecoprob present in probability form
vis.post <- function(ext1, 
                     ecotype, 
                     ecotrans, 
                     mass, 
                     cbp.long, 
                     time.start.stop,
                     ecoprob = TRUE) {

  # log-odds of occurrence associated with ecotype
  #   controlling for mass
  #   given group-level effects

  am <- melt(ext1$a)  # sim, time, state, value

  # translate the ecotype code into words
  suppressWarnings(am <- cbind(am, ecotrans[am$Var3, ]))
  names(am) <- c('sim', 'time', 'state', 'value', 'state_1', 'state_2')

  if(ecoprob) am$value <- invlogit(am$value)

  this.time <- time.start.stop[-nrow(time.start.stop), ]

  am$time <- mapvalues(x = am$time, 
                       from = unique(am$time), 
                       to = rev(this.time[, 2]))

  am$state_1 <- as.character(am$state_1)
  am$state_1 <- mapvalues(am$state_1, 
                          from = unique(am$state_1), 
                          to = c('carnivore', 'herbivore', 'insectivore', 
                                 'omnivore'))

  amplot <- ggplot(am, aes(x = time, y = value, group = sim))
  amplot <- amplot + geom_line(alpha = 0.01)
  amplot <- amplot + facet_grid(state_1 ~ state_2)
  amplot <- amplot + labs(x = 'Time (Mya)', 
                          y = 'Probability of occurrence')
  amplot <- amplot + scale_x_reverse()
  ggsave(filename = '../doc/figure/ecotype_occurrence.png', plot = amplot,
         width = 6, height = 4)


  # effect of group-level on individual-level
  tt <- ext1$gamma
  for(ii in seq(18)) {
    tt[, 4, ii] <- tt[, 1, ii] + tt[, 4, ii]
    tt[, 5, ii] <- tt[, 1, ii] + tt[, 5, ii]
  }
  #gm <- melt(ext1$gamma)  # sim, pred, ecotype, value
  gm <- melt(tt)  # sim, pred, ecotype, value

  # translate the ecotype code into words
  suppressWarnings(gm <- cbind(gm, ecotrans[gm$Var3, ]))
  names(gm) <- c('sim', 'predictor', 'ecotype', 'value', 'state_1', 'state_2')

  gm$predictor <- mapvalues(gm$predictor, 
                            from = unique(gm$predictor), 
                            to = c('phase 3', 'temp mean', 'temp range', 
                                   'phase 2', 'phase 1'))

  gm$state_1 <- as.character(gm$state_1)
  gm$state_1 <- mapvalues(gm$state_1, 
                          from = unique(gm$state_1), 
                          to = c('carnivore', 'herbivore', 'insectivore', 
                                 'omnivore'))

  gmplot <- ggplot(gm, aes(x = predictor, y = value, group = predictor))
  gmplot <- gmplot + geom_hline(yintercept = 0)
  gmplot <- gmplot + geom_violin()
  gmplot <- gmplot + facet_grid(state_1 ~ state_2)
  gmplot <- gmplot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  gmplot <- gmplot + labs(x = 'Predictor variable', 
                          y = 'Effect on log-odds of occurrence')
  ggsave(filename = '../doc/figure/group_on_ecotype.png', plot = gmplot,
         width = 6, height = 4)


  # variation in log-odds of occurrence associated with ecotype
  tm <- melt(ext1$tau)  # sample, ecotype, value

  suppressWarnings(tm <- cbind(tm, ecotrans[tm$Var2, ]))
  names(tm) <- c('sim', 'ecotype', 'value', 'state_1', 'state_2')
  tmplot <- ggplot(tm, aes(x = value, y = ..density..))
  tmplot <- tmplot + geom_histogram()
  tmplot <- tmplot + facet_grid(state_1 ~ state_2)
  tmplot <- tmplot + labs(x = 'Estimated standard deviation of distribution of ecotype log-odds of occurrence', 
                          y = 'Probability density')
  ggsave(filename = '../doc/figure/stdev_ecotype_occurrence.png', plot = tmplot,
         width = 6, height = 4)

  # difference from mean log-odds observation due to time
  pm <- melt(ext1$alpha_time)
  pm$Var2 <- mapvalues(x = pm$Var2, 
                       from = unique(pm$Var2), 
                       to = rev(time.start.stop[, 2]))

  pm$prob <- invlogit(pm$value)
  if(!is.null(pm$Var1)) {
    pmplot <- ggplot(pm, aes(x = Var2, y = prob, group = Var1))
  } else if (is.null(pm$Var1)) {
    pmplot <- ggplot(pm, aes(x = Var2, y = prob, group = iterations))
  }
  pmplot <- pmplot + geom_line(alpha = 0.01)
  pmplot <- pmplot + scale_x_reverse()
  pmplot <- pmplot + labs(x = 'Time (Mya)', 
                          y = 'Difference from mean log-odds observation')
  ggsave(filename = '../doc/figure/prob_preservation.png', plot = pmplot,
         width = 3, height = 2)


  # effect of mass on sampling
  mass.counter <- data.frame(x = seq(from = min(mass), 
                                     to = max(mass), 
                                     by = 0.1))
  mass.df <- data.frame(x = mass)
  mass.samp <- function(m) {
    i <- sample(length(ext1$alpha_0), 1)
    invlogit(ext1$alpha_1[i] * m + ext1$alpha_0[i])
  }
  mass_on_samp <- ggplot(mass.counter, aes(x = x))
  for(ii in seq(100)) {
    mass_on_samp <- mass_on_samp + stat_function(fun = mass.samp,
                                                 alpha = 0.1)
  }
  mass_on_samp <- mass_on_samp + geom_rug(data = mass.df, mapping = aes(x = x))
  mass_on_samp <- mass_on_samp + labs(x = 'Rescaled log mass (g)', 
                                      y = 'Probability of sampling')
  mass_on_samp <- mass_on_samp + theme(axis.title = element_text(size = 8),
                                       axis.text = element_text(size = 6))
  ggsave(filename = '../doc/figure/mass_on_samp.png', plot = mass_on_samp,
         width = 3, height = 2)
  # figure out on to plot on un-rescaled axis


  # effect of mass on observation
  # the idea is
  #   grid based on ecotype
  #   each panel has three counter-factuals
  #     one for each phase
  #   TODO rug for observed mass values
  mass.obs <- function(m, state, phase) {
    # intercept
    if(phase == 3) {
      cept <- median(ext1$gamma[, 1, state])
    } else if (phase == 2) {
      cept <- median(ext1$gamma[, 1, state]) + median(ext1$gamma[, 4, state])
    } else if (phase == 1) {
      cept <- median(ext1$gamma[, 1, state]) + median(ext1$gamma[, 5, state])
    }

    # bring it all together
    invlogit(cept + median(ext1$b_1) * m + median(ext1$b_2) * (m^2))
  }

  mass.counter <- seq(from = min(mass), to = max(mass), by = 0.001)
  out <- list()
  for(ii in seq(3)) {
    oo <- list()
    for(jj in seq(18)) {
      oo[[jj]] <- cbind(mass = mass.counter, 
                        est = mass.obs(mass.counter, jj, ii))
    }
    out[[ii]] <- oo
  }

  out <- llply(out, function(a) 
               Reduce(rbind, Map(function(x, y) cbind(x, y), 
                                 a, seq(length(a)))))

  out <- Reduce(rbind, Map(function(x, y) cbind(x, y), out, seq(3)))
  suppressWarnings(out <- data.frame(out, ecotrans[out[, 3], ]))
  names(out) <- c('mass', 'est', 'state', 'phase', 'state_1', 'state_2')

  out$state_1 <- as.character(out$state_1)
  out$state_1 <- mapvalues(out$state_1, 
                           from = unique(out$state_1), 
                           to = c('carnivore', 'herbivore', 'insectivore', 
                                  'omnivore'))

  mass_on_pres <- ggplot(out, aes(x = mass, y = est, colour = factor(phase)))
  mass_on_pres <- mass_on_pres + geom_line(size = 1.25)
  mass_on_pres <- mass_on_pres + facet_grid(state_1 ~ state_2)
  mass_on_pres <- mass_on_pres + labs(x = 'Rescaled log mass (g)',
                                      y = 'Probability of occurrence')
  mass_on_pres <- mass_on_pres + theme(axis.title = element_text(size = 10),
                                       axis.text = element_text(size = 8))
  mass_on_pres <- mass_on_pres + 
    scale_x_continuous(breaks = round(seq(min(out$mass), max(out$mass), 
                                          by = 1), 1))
  ggsave(filename = '../doc/figure/mass_on_pres.png', plot = mass_on_pres,
         width = 6, height = 5)

}





#' Visualize key aspects of the posterior distribution for the birth-death model
#'
#' @param ext1 extract from stanfit object
#' @param ecotype translate ecotype coding
#' @param ecotrans translate and break ecotype code
#' @param mass species mass estimates
#' @param cbp.long color pallette information for ggplot
#' @param ecoprob present in probability form
vis.bdpost <- function(ext2, 
                       ecotype, 
                       ecotrans, 
                       mass, 
                       cbp.long, 
                       time.start.stop,
                       order.cypher,
                       ecoprob = TRUE) {

  # log-odds of occurrence associated with ecotype
  #   controlling for mass
  #   given group-level effects

  ecotype.plot <- function(am, ecotrans, time.start.stop, survival = FALSE) {
    # origination
    # translate the ecotype code into words
    suppressWarnings(am <- cbind(am, ecotrans[am$Var3, ]))
    names(am) <- c('sim', 'time', 'state', 'value', 'state_1', 'state_2')
    if(ecoprob) am$value <- invlogit(am$value)

    if(survival) {
      orig.time <- time.start.stop[-nrow(time.start.stop), ]
    } else {
      orig.time <- time.start.stop
    }

    am$time <- mapvalues(x = am$time, 
                         from = unique(am$time), 
                         to = rev(orig.time[, 2]))

    am$state_1 <- as.character(am$state_1)
    am$state_1 <- mapvalues(am$state_1, 
                            from = unique(am$state_1), 
                            to = c('carnivore', 'herbivore', 'insectivore', 
                                   'omnivore'))

    amplot <- ggplot(am, aes(x = time, y = value, group = sim))
    amplot <- amplot + geom_line(alpha = 0.01)
    amplot <- amplot + facet_grid(state_1 ~ state_2)
    amplot <- amplot + scale_x_reverse()
    amplot
  }

  # origin
  am <- melt(ext2$o_a)  # sim, time, state, value
  amplot <- ecotype.plot(am, ecotrans, time.start.stop)
  amplot <- amplot + labs(x = 'Time (Mya)', 
                          y = 'Probability of originating')
  ggsave(filename = '../doc/figure/ecotype_origin_bd.png', plot = amplot,
         width = 6, height = 4)

  # survival
  am <- melt(ext2$s_a)  # sim, time, state, value
  amplot <- ecotype.plot(am, ecotrans, time.start.stop, survival = TRUE)
  amplot <- amplot + labs(x = 'Time (Mya)', 
                          y = 'Probability of survival')
  ggsave(filename = '../doc/figure/ecotype_survival_bd.png', plot = amplot,
         width = 6, height = 4)

  # preservation
  am <- melt(ext2$p_a)  # sim, time, state, value
  amplot <- ecotype.plot(am, ecotrans, time.start.stop)
  amplot <- amplot + labs(x = 'Time (Mya)', 
                          y = 'Probability of observation')
  ggsave(filename = '../doc/figure/ecotype_preserve_bd.png', plot = amplot,
         width = 6, height = 4)



  # effect of orders on origination
  ordeff.plot <- function(ord.eff, order.cypher) {
    names(ord.eff) <- c('state', 'value')
    ord.eff$state <- order.cypher[ord.eff$state]
    ordgg <- ggplot(ord.eff, aes(x = value, y = state))
    ordgg <- ordgg + geom_joy(rel_min_height = 0.01)
    ordgg <- ordgg + theme_joy()
    ordgg <- ordgg + scale_x_continuous(expand = c(0.01, 0))
    ordgg <- ordgg + scale_y_discrete(expand = c(0.01, 0))
    ordgg
  }

  ord.eff <- melt(ext2$o_ordeff)[, -1]
  ordgg <- ordeff.plot(ord.eff, order.cypher)
  ordgg <- ordgg + labs(x = 'Effect on origin',
                        y = 'Order')
  ggsave(filename = '../doc/figure/order_origin_bd.png', plot = ordgg,
         width = 6, height = 4)


  # effect of orders on survival 
  ord.eff <- melt(ext2$s_ordeff)[, -1]
  ordgg <- ordeff.plot(ord.eff, order.cypher)
  ordgg <- ordgg + labs(x = 'Effect on survival',
                        y = 'Order')
  ggsave(filename = '../doc/figure/order_survival_bd.png', plot = ordgg,
         width = 6, height = 4)

  # effect of orders on survival 
  ord.eff <- melt(ext2$p_ordeff)[, -1]
  ordgg <- ordeff.plot(ord.eff, order.cypher)
  ordgg <- ordgg + labs(x = 'Effect on observation',
                        y = 'Order')
  ggsave(filename = '../doc/figure/order_preserve_bd.png', plot = ordgg,
         width = 6, height = 4)



  # effect of group-level on individual-level
  groupeff.plot <- function(tt, ecotrans) {
    # origination
    for(ii in seq(18)) {
      tt[, 4, ii] <- tt[, 1, ii] + tt[, 4, ii]
      tt[, 5, ii] <- tt[, 1, ii] + tt[, 5, ii]
    }
    gm <- melt(tt)  # sim, pred, ecotype, value

    # translate the ecotype code into words
    suppressWarnings(gm <- cbind(gm, ecotrans[gm$Var3, ]))
    names(gm) <- c('sim', 'predictor', 'ecotype', 'value', 'state_1', 'state_2')

    gm$predictor <- mapvalues(gm$predictor, 
                              from = unique(gm$predictor), 
                              to = c('phase 3', 'temp mean', 'temp range', 
                                     'phase 2', 'phase 1'))

    gm$state_1 <- as.character(gm$state_1)
    gm$state_1 <- mapvalues(gm$state_1, 
                            from = unique(gm$state_1), 
                            to = c('carnivore', 'herbivore', 'insectivore', 
                                   'omnivore'))

    gmplot <- ggplot(gm, aes(x = predictor, y = value, group = predictor))
    gmplot <- gmplot + geom_hline(yintercept = 0)
    gmplot <- gmplot + geom_violin()
    gmplot <- gmplot + facet_grid(state_1 ~ state_2)
    gmplot <- gmplot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
    gmplot
  }

  tt <- ext2$o_gamma
  gmplot <- groupeff.plot(tt, ecotrans)
  gmplot <- gmplot + labs(x = 'Predictor variable', 
                          y = 'Change to log-odds of originating')
  ggsave(filename = '../doc/figure/group_on_origin_bd.png', plot = gmplot,
         width = 6, height = 4)

  # survival
  tt <- ext2$s_gamma
  gmplot <- groupeff.plot(tt, ecotrans)
  gmplot <- gmplot + labs(x = 'Predictor variable', 
                          y = 'Change to log-odds surviving')
  ggsave(filename = '../doc/figure/group_on_survival_bd.png', plot = gmplot,
         width = 6, height = 4)


  # observation
  tt <- ext2$p_gamma
  gmplot <- groupeff.plot(tt, ecotrans)
  gmplot <- gmplot + labs(x = 'Predictor variable', 
                          y = 'Change to log-odds of preserving')
  ggsave(filename = '../doc/figure/group_on_preserve_bd.png', plot = gmplot,
         width = 6, height = 4)



  # variation in log-odds of occurrence associated with ecotype
  # origination

  ecotypevar.plot <- function(tm, ecotrans) {
    suppressWarnings(tm <- cbind(tm, ecotrans[tm$Var2, ]))
    names(tm) <- c('sim', 'ecotype', 'value', 'state_1', 'state_2')
    tmplot <- ggplot(tm, aes(x = value, y = ..density..))
    tmplot <- tmplot + geom_histogram()
    tmplot <- tmplot + facet_grid(state_1 ~ state_2)
    tmplot
  }
  tm <- melt(ext2$o_tau)  # sample, ecotype, value
  tmplot <- ecotypevar.plot(tm, ecotrans)
  tmplot <- tmplot + labs(x = 'Estimated standard deviation of ecotype log-odds of originating', 
                          y = 'Probability density')
  ggsave(filename = '../doc/figure/stdev_ecotype_origin_bd.png', plot = tmplot,
         width = 6, height = 4)

  # survival
  tm <- melt(ext2$s_tau)  # sample, ecotype, value
  tmplot <- ecotypevar.plot(tm, ecotrans)
  tmplot <- tmplot + labs(x = 'Estimated standard deviation of ecotype log-odds of suriviving', 
                          y = 'Probability density')
  ggsave(filename = '../doc/figure/stdev_ecotype_survival_bd.png', plot = tmplot,
         width = 6, height = 4)


  #  # difference from mean log-odds observation due to time
  #  pm <- melt(ext2$alpha_time)
  #  pm$Var2 <- mapvalues(x = pm$Var2, 
  #                       from = unique(pm$Var2), 
  #                       to = rev(time.start.stop[, 2]))
  #
  #  pm$prob <- invlogit(pm$value)
  #  if(!is.null(pm$Var1)) {
  #    pmplot <- ggplot(pm, aes(x = Var2, y = prob, group = Var1))
  #  } else if (is.null(pm$Var1)) {
  #    pmplot <- ggplot(pm, aes(x = Var2, y = prob, group = iterations))
  #  }
  #  pmplot <- pmplot + geom_line(alpha = 0.01)
  #  pmplot <- pmplot + scale_x_reverse()
  #  pmplot <- pmplot + labs(x = 'Time (Mya)', 
  #                          y = 'Difference from mean log-odds observation')
  #  ggsave(filename = '../doc/figure/prob_preservation_bd.png', plot = pmplot,
  #         width = 3, height = 2)
  #
  #
  #  # effect of mass on sampling
  #  mass.counter <- data.frame(x = seq(from = min(mass), 
  #                                     to = max(mass), 
  #                                     by = 0.1))
  #  mass.df <- data.frame(x = mass)
  #  mass.samp <- function(m) {
  #    i <- sample(length(ext2$alpha_0), 1)
  #    invlogit(ext2$alpha_1[i] * m + ext2$alpha_0[i])
  #  }
  #  mass_on_samp <- ggplot(mass.counter, aes(x = x))
  #  for(ii in seq(100)) {
  #    mass_on_samp <- mass_on_samp + stat_function(fun = mass.samp,
  #                                                 alpha = 0.1)
  #  }
  #  mass_on_samp <- mass_on_samp + geom_rug(data = mass.df, mapping = aes(x = x))
  #  mass_on_samp <- mass_on_samp + labs(x = 'Rescaled log mass (g)', 
  #                                      y = 'Probability of sampling')
  #  mass_on_samp <- mass_on_samp + theme(axis.title = element_text(size = 8),
  #                                       axis.text = element_text(size = 6))
  #  ggsave(filename = '../doc/figure/mass_on_samp_bd.png', plot = mass_on_samp,
  #         width = 3, height = 2)
  #  # figure out on to plot on un-rescaled axis


  # effect of mass on observation
  # the idea is
  #   grid based on ecotype
  #   each panel has three counter-factuals
  #     one for each phase
  #   TODO rug for observed mass value
  mass.obs <- function(m, state, phase, type, ext2) {
    # intercept
    if(phase == 3) {
      if(type == 'origin') {
        cept <- median(ext2$o_gamma[, 1, state])
      } else if(type == 'survival') {
        cept <- median(ext2$s_gamma[, 1, state])
      } else if(type == 'preserve') {
        cept <- median(ext2$p_gamma[, 1, state])
      }
    } else if (phase == 2) {
      if(type == 'origin') {
        cept <- median(ext2$o_gamma[, 1, state]) + 
          median(ext2$o_gamma[, 4, state])
      } else if(type == 'survival') {
        cept <- median(ext2$s_gamma[, 1, state]) + 
          median(ext2$s_gamma[, 4, state])
      } else if(type == 'preserve') {
        cept <- median(ext2$p_gamma[, 1, state]) + 
          median(ext2$p_gamma[, 4, state])
      }
    } else if (phase == 1) {
      if(type == 'origin') {
        cept <- median(ext2$o_gamma[, 1, state]) + 
          median(ext2$o_gamma[, 5, state])
      } else if(type == 'survival') {
        cept <- median(ext2$s_gamma[, 1, state]) + 
          median(ext2$s_gamma[, 5, state])
      } else if(type == 'preserve') {
        cept <- median(ext2$p_gamma[, 1, state]) + 
          median(ext2$p_gamma[, 5, state])
      }
    }

    # bring it all together
    if(type == 'origin') {
      hold <- invlogit(cept + median(ext2$o_b_1) * m + 
                       median(ext2$o_b_2) * (m^2))
    } else if(type == 'survival') {
      hold <- invlogit(cept + median(ext2$s_b_1) * m + 
                       median(ext2$s_b_2) * (m^2))
    } else if(type == 'preserve') {
      hold <- invlogit(cept + median(ext2$p_b_1) * m + 
                       median(ext2$p_b_2) * (m^2))
    }
  }

  # origin
  massorigin.plot <- function(mass, ecotrans, type = 'origin', ext2) {
    mass.counter <- seq(from = min(mass), to = max(mass), by = 0.001)
    out <- list()
    for(ii in seq(3)) {
      oo <- list()
      for(jj in seq(18)) {
        oo[[jj]] <- cbind(mass = mass.counter, 
                          est = mass.obs(mass.counter, state = jj, phase = ii, 
                                         type = type, ext2 = ext2))
      }
      out[[ii]] <- oo
    }

    out <- llply(out, function(a) 
                 Reduce(rbind, Map(function(x, y) cbind(x, y), 
                                   a, seq(length(a)))))

    out <- Reduce(rbind, Map(function(x, y) cbind(x, y), out, seq(3)))
    suppressWarnings(out <- data.frame(out, ecotrans[out[, 3], ]))
    names(out) <- c('mass', 'est', 'state', 'phase', 'state_1', 'state_2')

    out$state_1 <- as.character(out$state_1)
    out$state_1 <- mapvalues(out$state_1, 
                             from = unique(out$state_1), 
                             to = c('carnivore', 'herbivore', 'insectivore', 
                                    'omnivore'))

    mass_on_pres <- ggplot(out, aes(x = mass, y = est, colour = factor(phase)))
    mass_on_pres <- mass_on_pres + geom_line(size = 1.25)
    mass_on_pres <- mass_on_pres + facet_grid(state_1 ~ state_2)
    mass_on_pres <- mass_on_pres + theme(axis.title = element_text(size = 10),
                                         axis.text = element_text(size = 8))
    mass_on_pres <- mass_on_pres + 
      scale_x_continuous(breaks = round(seq(min(out$mass), max(out$mass), 
                                            by = 1), 1))
    mass_on_pres
  }

  mass_on_pres <- massorigin.plot(mass, ecotrans, type = 'origin', ext2 = ext2)
  mass_on_pres <- mass_on_pres + labs(x = 'Rescaled log mass (g)',
                                      y = 'Probability of originating')
  ggsave(filename = '../doc/figure/mass_on_origin_bd.png', plot = mass_on_pres,
         width = 6, height = 4)

  # survival
  mass_on_pres <- massorigin.plot(mass, ecotrans, type = 'survival', ext2)
  mass_on_pres <- mass_on_pres + labs(x = 'Rescaled log mass (g)',
                                      y = 'Probability of survival')
  ggsave(filename = '../doc/figure/mass_on_surv_bd.png', plot = mass_on_pres,
         width = 6, height = 4)
  
  mass_on_pres <- massorigin.plot(mass, ecotrans, type = 'preserve', ext2)
  mass_on_pres <- mass_on_pres + labs(x = 'Rescaled log mass (g)',
                                      y = 'Probability of observation')
  ggsave(filename = '../doc/figure/mass_on_pres_bd.png', plot = mass_on_pres,
         width = 6, height = 4)

}

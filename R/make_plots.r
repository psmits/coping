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
    # logodds group with line for lod_odds based only climate 
    
    # calculate intercept
    # 18 possible ecotypes
    # 18 time points
    # time are rows, ecotypes are columns
    #ecotrans1 <- apply(ecotrans, 2, rev)
    comboname <- apply(ecotrans, 1, function(x) paste(x, collapse = '_'))
    if(!survival) {
      med <- ufull %*% apply(ext2$o_gamma, 2:3, median)
    } else if(survival) {
      med <- u %*% apply(ext2$s_gamma, 2:3, median)
    }
    colnames(med) <- comboname
    med <- melt(med)
    names(med) <- c('time', 'state', 'value')

    suppressWarnings(am <- cbind(am, ecotrans[am$Var3, ]))
    names(am) <- c('sim', 'time', 'state', 'value', 'state_1', 'state_2')

    # median line based on phase and temp
    if(ecoprob) am$value <- invlogit(am$value)

    if(survival) {
      orig.time <- time.start.stop[-nrow(time.start.stop), ]
    } else {
      orig.time <- time.start.stop
    }

    am$time <- mapvalues(x = am$time, 
                         from = unique(am$time), 
                         to = rev(orig.time[, 2]))
    med$time <- mapvalues(x = med$time,
                          from = unique(med$time),
                          to = rev(orig.time[, 2]))

    am$state_1 <- as.character(am$state_1)
    am$state_1 <- mapvalues(am$state_1, 
                            from = unique(am$state_1), 
                            to = c('carnivore', 'herbivore', 'insectivore', 
                                   'omnivore'))

    med$state <- as.character(med$state)
    med <- cbind(med, Reduce(rbind, str_split(med$state, '\\_')))
    names(med) <- c('time', 'state', 'value', 'state_1', 'state_2')
    med$state_1 <- as.character(med$state_1)
    med$state_1 <- mapvalues(med$state_1, 
                             from = unique(med$state_1), 
                             to = c('carnivore', 'herbivore', 'insectivore', 
                                    'omnivore'))

    if(ecoprob) med$value <- invlogit(med$value)

    amplot <- ggplot(am, aes(x = time, y = value, group = sim))
    amplot <- amplot + geom_line(data = med, 
                                 mapping = aes(x = time, y = value, group = NULL), 
                                 colour = 'blue', size = 1)
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
                        y = 'Taxonomic Order')
  ggsave(filename = '../doc/figure/order_origin_bd.png', plot = ordgg,
         width = 6, height = 4)


  # effect of orders on survival 
  ord.eff <- melt(ext2$s_ordeff)[, -1]
  ordgg <- ordeff.plot(ord.eff, order.cypher)
  ordgg <- ordgg + labs(x = 'Effect on survival',
                        y = 'Taxonomic Order')
  ggsave(filename = '../doc/figure/order_survival_bd.png', plot = ordgg,
         width = 6, height = 4)



  # effect of group-level on individual-level
  groupeff.plot <- function(tt, ecotrans) {
    # origination
    for(ii in seq(18)) {
      tt[, 2, ii] <- tt[, 1, ii] + tt[, 2, ii]
      tt[, 3, ii] <- tt[, 1, ii] + tt[, 3, ii]
    }
    gm <- melt(tt)  # sim, pred, ecotype, value

    # translate the ecotype code into words
    suppressWarnings(gm <- cbind(gm, ecotrans[gm$Var3, ]))
    names(gm) <- c('sim', 'predictor', 'ecotype', 'value', 'state_1', 'state_2')

    gm$predictor <- mapvalues(gm$predictor, 
                              from = unique(gm$predictor), 
                              to = c('phase 1', 'phase 2', 
                                     'phase 3', 'mean temp.'))

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

  # effect of mass on observation
  # the idea is
  #   grid based on ecotype
  #   each panel has three counter-factuals
  #     one for each phase
  #   TODO rug for observed mass value
  mass.obs <- function(m, state, phase, type, ext2) {
    # intercept
    if(phase == 1) {
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
          median(ext2$o_gamma[, 2, state])
      } else if(type == 'survival') {
        cept <- median(ext2$s_gamma[, 1, state]) + 
          median(ext2$s_gamma[, 2, state])
      } else if(type == 'preserve') {
        cept <- median(ext2$p_gamma[, 1, state]) + 
          median(ext2$p_gamma[, 2, state])
      }
    } else if (phase == 3) {
      if(type == 'origin') {
        cept <- median(ext2$o_gamma[, 1, state]) + 
          median(ext2$o_gamma[, 3, state])
      } else if(type == 'survival') {
        cept <- median(ext2$s_gamma[, 1, state]) + 
          median(ext2$s_gamma[, 3, state])
      } else if(type == 'preserve') {
        cept <- median(ext2$p_gamma[, 1, state]) + 
          median(ext2$p_gamma[, 3, state])
      }
    }

    # bring it all together
    if(type == 'origin') {
      hold <- invlogit(cept + median(ext2$o_b_1) * m)
    } else if(type == 'survival') {
      hold <- invlogit(cept + median(ext2$s_b_1) * m)
    } else if(type == 'preserve') {
      hold <- invlogit(cept + median(ext2$p_b_1) * m)
    }
    hold
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
  
  ## observation
  # this has a slightly different graph now
  #mass_on_pres <- massorigin.plot(mass, ecotrans, type = 'preserve', ext2)
  #mass_on_pres <- mass_on_pres + labs(x = 'Rescaled log mass (g)',
  #                                    y = 'Probability of survival')
  #ggsave(filename = '../doc/figure/mass_on_pres_bd.png', plot = mass_on_pres,
  #       width = 6, height = 4)

  mass.counter <- data.frame(x = seq(from = min(mass), to = max(mass), by = 0.001))
  mass.obs.foo <- function(x) 
    invlogit(median(ext2$p_inter) + median(ext2$p_b_1) * x)
  massobs <- ggplot(mass.counter, aes(x = x))
  massobs <- massobs + stat_function(fun = mass.obs.foo, size = 2)
  massobs <- massobs + labs(x = 'Std. from mean mass (g)',
                            y = 'Probability of observting is present')
  ggsave(filename = '../doc/figure/mass_on_pres_bd.png', plot = massobs,
         width = 6, height = 4)



  # effects of various single factors on preservation
  # and their relative variances (possible size of effects)
  ext2$p_0  # constant intercept
  ext2$p_b_1  # mass slope

  # functional group effect
  ecotrans1 <- ecotrans
  ecotrans1[, 1] <- mapvalues(ecotrans[, 1], 
                              from = unique(ecotrans[, 1]), 
                              to = c('carnivore', 'herbivore', 
                                     'insectivore', 'omnivore'))
  ecotrans1 <- apply(ecotrans1, 2, rev)
  comboname <- apply(ecotrans1, 1, function(x) paste(x, collapse = ' '))
  fe <- melt(ext2$p_funceff)[, -1]
  names(fe) <- c('state', 'value')
  fe$state <- comboname[fe$state]
  fegg <- ggplot(fe, aes(x = value, y = state))
  fegg <- fegg + geom_joy(rel_min_height = 0.01)
  fegg <- fegg + theme_joy()
  fegg <- fegg + scale_x_continuous(expand = c(0.01, 0))
  fegg <- fegg + scale_y_discrete(expand = c(0.01, 0))
  fegg <- fegg + labs(x = 'Effect on log-odds observation',
                      y = 'Functional group')
  ggsave(filename = '../doc/figure/ecotype_observation.png', plot = fegg,
         width = 6, height = 4)
 
  ## order effect
  #oe <- melt(ext2$p_ordeff)[, -1]
  #names(oe) <- c('state', 'value')
  #oe$state <- order.cypher[oe$state]
  #oegg <- ggplot(oe, aes(x = value, y = state))
  #oegg <- oegg + geom_joy(rel_min_height = 0.01)
  #oegg <- oegg + theme_joy()
  #oegg <- oegg + scale_x_continuous(expand = c(0.01, 0))
  #oegg <- oegg + scale_y_discrete(expand = c(0.01, 0))
  #oegg <- oegg + labs(x = 'Effect on log-odds of observation',
  #                    y = 'Taxonomic Order')
  #ggsave(filename = '../doc/figure/order_observation.png', plot = oegg,
  #       width = 6, height = 4)

  # time effect
  te <- melt(ext2$p_timeeff)[, -1]
  names(te) <- c('time', 'value')
  te$time <- time.start.stop[te$time, 2]
  te$time <- mapvalues(te$time, from = nalma$ma, to = nalma$interval)
  te$time <- factor(te$time, levels = unique(te$time))
  tegg <- ggplot(te, aes(x = value, y = time))
  tegg <- tegg + geom_joy(rel_min_height = 0.01)
  tegg <- tegg + theme_joy()
  tegg <- tegg + scale_x_continuous(expand = c(0.01, 0))
  tegg <- tegg + scale_y_discrete(expand = c(0.01, 0))
  tegg <- tegg + labs(x = 'Effect on log-odds of observation',
                      y = 'NALMA')
  ggsave(filename = '../doc/figure/time_observation.png', plot = tegg,
         width = 6, height = 4)

  # relative variance
  effscale <- melt(data.frame(functional = ext2$p_funcscale, 
                              #order = ext2$p_ordscale, 
                              time = ext2$p_timescale))
  names(effscale) <- c('source', 'value')
  scgg <- ggplot(effscale, aes(x = value, y = source))
  scgg <- scgg + geom_joy(rel_min_height = 0.01)
  scgg <- scgg + theme_joy()
  scgg <- scgg + scale_x_continuous(expand = c(0.01, 0))
  scgg <- scgg + scale_y_discrete(expand = c(0.01, 0))
  scgg <- scgg + labs(x = 'Standard deviation of effect',
                      y = 'Effect source')
  ggsave(filename = '../doc/figure/scales_observation.png', plot = scgg,
         width = 6, height = 4)

}

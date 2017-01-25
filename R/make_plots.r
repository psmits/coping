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
                     ecoprob = TRUE) {

  # log-odds of occurrence associated with ecotype
  #   controlling for mass
  #   given group-level effects

  am <- melt(ext1$a)  # sim, time, state, value

  # translate the ecotype code into words
  suppressWarnings(am <- cbind(am, ecotrans[am$Var3, ]))
  names(am) <- c('sim', 'time', 'state', 'value', 'state_1', 'state_2')

  if(ecoprob) am$value <- invlogit(am$value)

  amplot <- ggplot(am, aes(x = time, y = value, group = sim))
  amplot <- amplot + geom_line(alpha = 0.01)
  amplot <- amplot + facet_grid(state_1 ~ state_2)
  amplot <- amplot + labs(x = 'Time (Mya)', 
                          y = 'Probability of occurrence')
  ggsave(filename = '../doc/figure/ecotype_occurrence_aug.png', plot = amplot,
         width = 6, height = 4)
  
  sam <- am[am$state_1 != 'augment', ]
  samplot <- amplot %+% sam
  ggsave(filename = '../doc/figure/ecotype_occurrence.png', plot = samplot,
         width = 6, height = 4)



  # effect of group-level on individual-level
  gm <- melt(ext1$gamma)  # sim, pred, ecotype, value

  # translate the ecotype code into words
  suppressWarnings(gm <- cbind(gm, ecotrans[gm$Var3, ]))
  names(gm) <- c('sim', 'predictor', 'ecotype', 'value', 'state_1', 'state_2')

  gm$predictor <- mapvalues(gm$predictor, 
                            from = unique(gm$predictor), 
                            to = c('phase 3', 'temp mean', 'temp range', 
                                   'phase 2', 'phase 1'))
  gmplot <- ggplot(gm, aes(x = predictor, y = value, group = predictor))
  gmplot <- gmplot + geom_hline(yintercept = 0)
  gmplot <- gmplot + geom_violin()
  gmplot <- gmplot + facet_grid(state_1 ~ state_2)
  gmplot <- gmplot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  gmplot <- gmplot + labs(x = 'Predictor variable', 
                          y = 'Effect on log-odds of occurrence')
  ggsave(filename = '../doc/figure/group_on_ecotype_aug.png', plot = gmplot,
         width = 6, height = 4)
  
  gam <- gm[gm$state_1 != 'augment', ]
  samplot <- gmplot %+% gam
  ggsave(filename = '../doc/figure/group_on_ecotype.png', plot = samplot,
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
  ggsave(filename = '../doc/figure/stdev_ecotype_occurrence_aug.png', plot = tmplot,
         width = 6, height = 4)

  gam <- tm[tm$state_1 != 'augment', ]
  samplot <- tmplot %+% gam
  ggsave(filename = '../doc/figure/stdev_ecotype_occurrence.png', plot = samplot,
         width = 6, height = 4)

  # difference from mean log-odds observation due to time
  pm <- melt(ext1$alpha_time)

  pm$prob <- invlogit(pm$value)
  if(!is.null(pm$Var1)) {
    pmplot <- ggplot(pm, aes(x = Var2, y = prob, group = Var1))
  } else if (is.null(pm$Var1)) {
    pmplot <- ggplot(pm, aes(x = Var2, y = prob, group = iterations))
  }
  pmplot <- pmplot + geom_line(alpha = 0.01)
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

  mass_on_pres <- ggplot(out, aes(x = mass, y = est, colour = factor(phase)))
  mass_on_pres <- mass_on_pres + geom_line(size = 1.25)
  mass_on_pres <- mass_on_pres + facet_grid(state_1 ~ state_2)
  mass_on_pres <- mass_on_pres + labs(x = 'Rescaled log mass (g)',
                                      y = 'Probability of occurrence')
  ggsave(filename = '../doc/figure/mass_on_pres_aug.png', plot = mass_on_pres,
         width = 6, height = 4)

  sam <- out[out$state_1 != 'augment', ]
  samplot <- mass_on_pres %+% sam
  ggsave(filename = '../doc/figure/mass_on_pres.png', plot = samplot,
         width = 6, height = 4)


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
                     ecoprob = TRUE) {

  # log-odds of occurrence associated with ecotype
  #   controlling for mass
  #   given group-level effects

  # origination
  am <- melt(ext2$o_a)  # sim, time, state, value
  # translate the ecotype code into words
  suppressWarnings(am <- cbind(am, ecotrans[am$Var3, ]))
  names(am) <- c('sim', 'time', 'state', 'value', 'state_1', 'state_2')
  if(ecoprob) am$value <- invlogit(am$value)

  amplot <- ggplot(am, aes(x = time, y = value, group = sim))
  amplot <- amplot + geom_line(alpha = 0.01)
  amplot <- amplot + facet_grid(state_1 ~ state_2)
  amplot <- amplot + labs(x = 'Time (Mya)', 
                          y = 'Probability of originating')
  ggsave(filename = '../doc/figure/ecotype_origin_bd_aug.png', plot = amplot,
         width = 6, height = 4)
 
  sam <- am[am$state_1 != 'augment', ]
  samplot <- amplot %+% sam
  ggsave(filename = '../doc/figure/ecotype_origin_bd.png', plot = samplot,
         width = 6, height = 4)

  # survival
  am <- melt(ext2$s_a)  # sim, time, state, value
  # translate the ecotype code into words
  suppressWarnings(am <- cbind(am, ecotrans[am$Var3, ]))
  names(am) <- c('sim', 'time', 'state', 'value', 'state_1', 'state_2')
  if(ecoprob) am$value <- invlogit(am$value)

  amplot <- ggplot(am, aes(x = time, y = value, group = sim))
  amplot <- amplot + geom_line(alpha = 0.01)
  amplot <- amplot + facet_grid(state_1 ~ state_2)
  amplot <- amplot + labs(x = 'Time (Mya)', 
                          y = 'Probability of surviving')
  ggsave(filename = '../doc/figure/ecotype_survival_bd_aug.png', plot = amplot,
         width = 6, height = 4)

  sam <- am[am$state_1 != 'augment', ]
  samplot <- amplot %+% sam
  ggsave(filename = '../doc/figure/ecotype_survival_bd.png', plot = samplot,
         width = 6, height = 4)



  # effect of group-level on individual-level
  # origination
  gm <- melt(ext2$o_gamma)  # sim, pred, ecotype, value
  # translate the ecotype code into words
  suppressWarnings(gm <- cbind(gm, ecotrans[gm$Var3, ]))
  names(gm) <- c('sim', 'predictor', 'ecotype', 'value', 'state_1', 'state_2')

  gm$predictor <- mapvalues(gm$predictor, 
                            from = unique(gm$predictor), 
                            to = c('phase 3', 'temp mean', 'temp range', 
                                   'phase 2', 'phase 1'))
  gmplot <- ggplot(gm, aes(x = predictor, y = value, group = predictor))
  gmplot <- gmplot + geom_hline(yintercept = 0)
  gmplot <- gmplot + geom_violin()
  gmplot <- gmplot + facet_grid(state_1 ~ state_2)
  gmplot <- gmplot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  gmplot <- gmplot + labs(x = 'Predictor variable', 
                          y = 'Change to log-odds of originating')
  ggsave(filename = '../doc/figure/group_on_origin_bd_aug.png', plot = gmplot,
         width = 6, height = 4)
  
  gam <- gm[gm$state_1 != 'augment', ]
  samplot <- gmplot %+% gam
  ggsave(filename = '../doc/figure/group_on_origin_bd.png', plot = samplot,
         width = 6, height = 4)


  # survival
  gm <- melt(ext2$s_gamma)  # sim, pred, ecotype, value
  # translate the ecotype code into words
  suppressWarnings(gm <- cbind(gm, ecotrans[gm$Var3, ]))
  names(gm) <- c('sim', 'predictor', 'ecotype', 'value', 'state_1', 'state_2')

  gm$predictor <- mapvalues(gm$predictor, 
                            from = unique(gm$predictor), 
                            to = c('phase 3', 'temp mean', 'temp range', 
                                   'phase 2', 'phase 1'))
  gmplot <- ggplot(gm, aes(x = predictor, y = value, group = predictor))
  gmplot <- gmplot + geom_hline(yintercept = 0)
  gmplot <- gmplot + geom_violin()
  gmplot <- gmplot + facet_grid(state_1 ~ state_2)
  gmplot <- gmplot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  gmplot <- gmplot + labs(x = 'Predictor variable', 
                          y = 'Change to log-odds surviving')
  ggsave(filename = '../doc/figure/group_on_survival_bd_aug.png', plot = gmplot,
         width = 6, height = 4)
  
  gam <- gm[gm$state_1 != 'augment', ]
  samplot <- gmplot %+% gam
  ggsave(filename = '../doc/figure/group_on_survival_bd.png', plot = samplot,
         width = 6, height = 4)


  # variation in log-odds of occurrence associated with ecotype
  # origination
  tm <- melt(ext2$o_tau)  # sample, ecotype, value
  suppressWarnings(tm <- cbind(tm, ecotrans[tm$Var2, ]))
  names(tm) <- c('sim', 'ecotype', 'value', 'state_1', 'state_2')
  tmplot <- ggplot(tm, aes(x = value, y = ..density..))
  tmplot <- tmplot + geom_histogram()
  tmplot <- tmplot + facet_grid(state_1 ~ state_2)
  tmplot <- tmplot + labs(x = 'Estimated standard deviation of ecotype log-odds of originating', 
                          y = 'Probability density')
  ggsave(filename = '../doc/figure/stdev_ecotype_origin_bd_aug.png', plot = tmplot,
         width = 6, height = 4)
  
  gam <- tm[tm$state_1 != 'augment', ]
  samplot <- tmplot %+% gam
  ggsave(filename = '../doc/figure/stdev_ecotype_origin_bd.png', plot = samplot,
         width = 6, height = 4)
  
  # survival
  tm <- melt(ext2$s_tau)  # sample, ecotype, value
  suppressWarnings(tm <- cbind(tm, ecotrans[tm$Var2, ]))
  names(tm) <- c('sim', 'ecotype', 'value', 'state_1', 'state_2')
  tmplot <- ggplot(tm, aes(x = value, y = ..density..))
  tmplot <- tmplot + geom_histogram()
  tmplot <- tmplot + facet_grid(state_1 ~ state_2)
  tmplot <- tmplot + labs(x = 'Estimated standard deviation of ecotype log-odds of suriviving', 
                          y = 'Probability density')
  ggsave(filename = '../doc/figure/stdev_ecotype_survival_bd_aug.png', plot = tmplot,
         width = 6, height = 4)
  
  gam <- tm[tm$state_1 != 'augment', ]
  samplot <- tmplot %+% gam
  ggsave(filename = '../doc/figure/stdev_ecotype_survival_bd.png', plot = samplot,
         width = 6, height = 4)
  


  # difference from mean log-odds observation due to time
  pm <- melt(ext2$alpha_time)

  pm$prob <- invlogit(pm$value)
  if(!is.null(pm$Var1)) {
    pmplot <- ggplot(pm, aes(x = Var2, y = prob, group = Var1))
  } else if (is.null(pm$Var1)) {
    pmplot <- ggplot(pm, aes(x = Var2, y = prob, group = iterations))
  }
  pmplot <- pmplot + geom_line(alpha = 0.01)
  pmplot <- pmplot + labs(x = 'Time (Mya)', 
                          y = 'Difference from mean log-odds observation')
  ggsave(filename = '../doc/figure/prob_preservation_bd.png', plot = pmplot,
         width = 3, height = 2)


  # effect of mass on sampling
  mass.counter <- data.frame(x = seq(from = min(mass), 
                                     to = max(mass), 
                                     by = 0.1))
  mass.df <- data.frame(x = mass)
  mass.samp <- function(m) {
    i <- sample(length(ext2$alpha_0), 1)
    invlogit(ext2$alpha_1[i] * m + ext2$alpha_0[i])
  }
  mass_on_samp <- ggplot(mass.counter, aes(x = x))
  for(ii in seq(100)) {
    mass_on_samp <- mass_on_samp + stat_function(fun = mass.samp,
                                                 alpha = 0.1)
  }
  mass_on_samp <- mass_on_samp + geom_rug(data = mass.df, mapping = aes(x = x))
  mass_on_samp <- mass_on_samp + labs(x = 'Rescaled log mass (g)', 
                                      y = 'Probability of sampling')
  ggsave(filename = '../doc/figure/mass_on_samp_bd.png', plot = mass_on_samp,
         width = 3, height = 2)
  # figure out on to plot on un-rescaled axis


  # effect of mass on observation
  # the idea is
  #   grid based on ecotype
  #   each panel has three counter-factuals
  #     one for each phase
  #   TODO rug for observed mass value
  mass.obs <- function(m, state, phase, origin = TRUE) {
    # intercept
    if(phase == 3) {
      if(origin) {
        cept <- median(ext2$o_gamma[, 1, state])
      } else {
        cept <- median(ext2$s_gamma[, 1, state])
      }
    } else if (phase == 2) {
      if(origin) {
        cept <- median(ext2$o_gamma[, 1, state]) + 
          median(ext2$o_gamma[, 4, state])
      } else {
        cept <- median(ext2$s_gamma[, 1, state]) + 
          median(ext2$s_gamma[, 4, state])
      }
    } else if (phase == 1) {
      if(origin) {
        cept <- median(ext2$o_gamma[, 1, state]) + 
          median(ext2$o_gamma[, 5, state])
      } else {
        cept <- median(ext2$s_gamma[, 1, state]) + 
          median(ext2$s_gamma[, 5, state])
      }
    }

    # bring it all together
    if(origin) {
      hold <- invlogit(cept + median(ext2$o_b_1) * m + 
                       median(ext2$o_b_2) * (m^2))
    } else {
      hold <- invlogit(cept + median(ext2$s_b_1) * m + 
                       median(ext2$s_b_2) * (m^2))
    }
    hold
  }
  
  # origin
  mass.counter <- seq(from = min(mass), to = max(mass), by = 0.001)
  out <- list()
  for(ii in seq(3)) {
    oo <- list()
    for(jj in seq(19)) {
      oo[[jj]] <- cbind(mass = mass.counter, 
                        est = mass.obs(mass.counter, jj, ii, origin = TRUE))
    }
    out[[ii]] <- oo
  }

  out <- llply(out, function(a) 
               Reduce(rbind, Map(function(x, y) cbind(x, y), 
                                 a, seq(length(a)))))

  out <- Reduce(rbind, Map(function(x, y) cbind(x, y), out, seq(3)))
  suppressWarnings(out <- data.frame(out, ecotrans[out[, 3], ]))
  names(out) <- c('mass', 'est', 'state', 'phase', 'state_1', 'state_2')

  mass_on_pres <- ggplot(out, aes(x = mass, y = est, colour = factor(phase)))
  mass_on_pres <- mass_on_pres + geom_line(size = 1.25)
  mass_on_pres <- mass_on_pres + facet_grid(state_1 ~ state_2)
  mass_on_pres <- mass_on_pres + labs(x = 'Rescaled log mass (g)',
                                      y = 'Probability of occurrence')
  ggsave(filename = '../doc/figure/mass_on_origin_bd_aug.png', plot = mass_on_pres,
         width = 6, height = 4)

  sam <- out[out$state_1 != 'augment', ]
  samplot <- mass_on_pres %+% sam
  ggsave(filename = '../doc/figure/mass_on_origin_bd.png', plot = samplot,
         width = 6, height = 4)

  
  # survival
  mass.counter <- seq(from = min(mass), to = max(mass), by = 0.001)
  out <- list()
  for(ii in seq(3)) {
    oo <- list()
    for(jj in seq(19)) {
      oo[[jj]] <- cbind(mass = mass.counter, 
                        est = mass.obs(mass.counter, jj, ii, origin = FALSE))
    }
    out[[ii]] <- oo
  }

  out <- llply(out, function(a) 
               Reduce(rbind, Map(function(x, y) cbind(x, y), 
                                 a, seq(length(a)))))

  out <- Reduce(rbind, Map(function(x, y) cbind(x, y), out, seq(3)))
  suppressWarnings(out <- data.frame(out, ecotrans[out[, 3], ]))
  names(out) <- c('mass', 'est', 'state', 'phase', 'state_1', 'state_2')

  mass_on_pres <- ggplot(out, aes(x = mass, y = est, colour = factor(phase)))
  mass_on_pres <- mass_on_pres + geom_line(size = 1.25)
  mass_on_pres <- mass_on_pres + facet_grid(state_1 ~ state_2)
  mass_on_pres <- mass_on_pres + labs(x = 'Rescaled log mass (g)',
                                      y = 'Probability of survival')
  ggsave(filename = '../doc/figure/mass_on_surv_bd_aug.png', plot = mass_on_pres,
         width = 6, height = 4)

  sam <- out[out$state_1 != 'augment', ]
  samplot <- mass_on_pres %+% sam
  ggsave(filename = '../doc/figure/mass_on_surv_bd.png', plot = samplot,
         width = 6, height = 4)


}

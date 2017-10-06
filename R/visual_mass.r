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


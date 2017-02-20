#' Simulate from the turnover_revamp model
#'
#' Works either for prior or posterior predictive stuff.
#'
#' @param ntax int number of taxa
#' @param ntime int number of time slices
#' @param phi starting condition probability
#' @param pred matrix of predictors for N x T
#' @param p matrix N x T of prob preservation
#' @param death logical for if death is absorbing state
#' @return list of z (true), y (observed)
#' @export
model.simulation <- function(ntax, ntime, phi, pred, p = NULL, death = TRUE) {
  y <- z <- matrix(NA, nrow = ntax, ncol = ntime)

  for(nn in 1:ntax) {
    z[nn, 1] <- rbinom(1, 1, prob = phi)

    if(!is.null(p)) 
      y[nn, 1] <- rbinom(1, 1, prob = z[nn, 1] * p[nn, 1])

    prod.term <- 1 - z[nn, 1]
    for(tt in 2:ntime) {
      prod.term <- prod.term * (1 - z[nn, tt - 1])
      if(death == TRUE) {
        z[nn, tt] <- rbinom(1, 1, prob = z[nn, tt - 1] * pred[nn, tt - 1] + 
                            prod.term * pred[nn, tt - 1])
      } else {
        z[nn, tt] <- rbinom(1, 1, prob = z[nn, tt - 1] * pred[nn, tt - 1] + 
                            (1 - z[nn, tt - 1]) * pred[nn, tt - 1])
      }
      
      # preservation
      if(!is.null(p)) 
        y[nn, tt] <- rbinom(1, 1, prob = z[nn, tt] * p[nn, tt])
    }
  }
  if(is.null(p)) y <- NULL

  out <- list(z = z, y = y)
  out
}




#' Simulate from the turnover_revamp_2 model
#'
#' Works either for prior or posterior predictive stuff.
#'
#' @param ntax int number of taxa
#' @param ntime int number of time slices
#' @param phi starting condition probability
#' @param origin matrix of predictors for N x T
#' @param stay matrix of predictors for N x T
#' @param p matrix N x T of prob preservation
#' @param death logical for if death is absorbing state
#' @return list of z (true), y (observed)
#' @export
bdmodel.simulation <- function(ntax, ntime, origin, stay, p = NULL, death = TRUE) {
  y <- z <- matrix(NA, nrow = ntax, ncol = ntime)

  for(nn in 1:ntax) {
    z[nn, 1] <- rbinom(1, 1, prob = origin[nn, 1])

    if(!is.null(p)) 
      y[nn, 1] <- rbinom(1, 1, prob = z[nn, 1] * p[nn, 1])

    prod.term <- 1 - z[nn, 1]
    for(tt in 2:ntime) {
      prod.term <- prod.term * (1 - z[nn, tt - 1])
      if(death == TRUE) {
        z[nn, tt] <- rbinom(1, 1, prob = z[nn, tt - 1] * stay[nn, tt - 1] + 
                            prod.term * origin[nn, tt])
      } else {
        z[nn, tt] <- rbinom(1, 1, prob = z[nn, tt - 1] * stay[nn, tt - 1] + 
                            (1 - z[nn, tt - 1]) * origin[nn, tt])
      }
      
      # preservation
      if(!is.null(p)) 
        y[nn, tt] <- rbinom(1, 1, prob = z[nn, tt] * p[nn, tt])
    }
  }
  if(is.null(p)) y <- NULL

  out <- list(z = z, y = y)
  out
}




#' Posterior predictive checks (suite)
#'
#' @param ext1 output from stan
#' @param ntax number of species
#' @param ntime number of time points
#' @param sight.obs empirical sighting matrix
#' @param nsim number of simulations
#' @param samp
post.pred <- function(ext1, ntax, ntime, sight.obs, nsim, samp, bd = FALSE) {
  sim.obs <- list()
  for(ii in seq(nsim)) {
    if(bd) {
      sim.obs[[ii]] <- bdmodel.simulation(ntax, ntime, 
                                          origin = ext1$origin[samp[ii], , ],
                                          stay = ext1$stay[samp[ii], , ],
                                          p = ext1$p[samp[ii], , ])
    } else if(!bd) {
      sim.obs[[ii]] <- model.simulation(ntax, ntime, phi = ext1$phi[samp[ii]], 
                                        pred = ext1$pred[samp[ii], , ],
                                        p = ext1$p[samp[ii], , ])
    }
  }

  meanocc.obs <- mean(rowSums(sight.obs))
  meanocc.simobs <- laply(sim.obs, function(x) mean(rowSums(x$y)))

  mos <- data.frame(x = meanocc.simobs, 
                    y = rep('Full', nsim))
  obs <- data.frame(x = meanocc.obs, 
                    y = rep('Full', nsim))
  ocplot <- ggplot(mos, aes(x = x)) + geom_histogram()
  ocplot <- ocplot + geom_vline(data = obs, mapping = aes(xintercept = x), 
                                colour = 'blue', size = 1.5)
  ocplot <- ocplot + labs(x = 'Mean obs. per species',
                          y = 'Post. pred. simulations')

  if(bd) {
    plotname <- paste0('../doc/figure/pred_occ_bd.png') 
  } else {
    plotname <- paste0('../doc/figure/pred_occ.png') 
  }

  ggsave(filename = plotname, plot = ocplot,
         width = 3, height = 2.5)

}

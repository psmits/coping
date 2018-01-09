#' Estimate standing diversity
#'
#' This function only makes one diversity simulation. 
#' Wrap it to create a distribution of diveristy simulations.
#' Currently only works for the birth-death model.
#'
#' @param data mtraix of observed occurrences
#' @param posterior stanfit object, extracted
estimate.diversity <- function(data, posterior) {
  #data <- sighta
  #posterior <- ext2
  #nsim <- 1
  grab <- sample(length(posterior$lp__), 1) 

  # range through matrix
  z <- matrix(0, ncol = ncol(data), nrow = nrow(data))
  for(ii in seq(nrow(data))) {
    if(all(data[ii, ] != 0)) {
      rant <- range(which(data[ii, ] == 1))
      z[ii, rant[1]:rant[2]] <- 1
    }
  }

  for(ii in seq(nrow(z))) {
    if(all(data[ii, ] != 0)) {
      fad <- min(which(z[ii, ] == 1))

      if(fad > 1) {
        for(jj in seq(fad - 1)) {
          if(jj == 1) {
            z[ii, jj] <- rbinom(1, 1, prob = posterior$origin[grab, ii, jj])
          } else if(jj > 1) {
            if(z[ii, jj - 1] == 0) {
              z[ii, jj] <- rbinom(1, 1, prob = posterior$origin[grab, ii, jj])
            } else if(z[ii, jj - 1] == 1) {
              z[ii, jj] <- 1
            }
          }
        }
      }

      lad <- max(which(z[ii, ] == 1))
      if(lad < ncol(z)) {
        for(jj in (lad + 1):length(z[ii, ])) {
          if(z[ii, jj - 1] == 1) {
            z[ii, jj] <- rbinom(1, 1, prob = posterior$stay[grab, ii, jj - 1])
          } else if(z[ii, jj - 1] == 0) {
            z[ii, jj] <- 0
          }
        }
      }
    } else { 
      z[ii, 1] <- rbinom(1, 1, prob = posterior$origin[grab, ii, 1])
      prod.term <- 1 - z[ii, 1]
      for(tt in 2:(ncol(z))) {
        prod.term <- prod.term * (1 - z[ii, tt - 1])
        # never observed
        pp <- z[ii, tt - 1] * posterior$stay[grab, ii, tt - 1] + 
          prod.term * posterior$origin[grab, ii, tt]
        z[ii, tt] <- rbinom(1, 1, prob = pp)
      }
    }
  }

  return(z)
}

#' Wrapper for estimate.diveristy to give posterior of diversity from birth-death model
#'
#' @param data mtraix of observed occurrences
#' @param posterior stanfit object, extracted
#' @param nsim integer, number of posterior simulations to do
diversity.distribution <- function(data, posterior, nsim) {
  out <- vector(mode = 'list', length = nsim)
  if(is.na(detectCores())) {
    out <- Map(function(x) estimate.diversity(data, posterior), out)
  } else {
    out <- mcMap(function(x) estimate.diversity(data, posterior), 
                 out, mc.cores = detectCores())
  }
  out
}

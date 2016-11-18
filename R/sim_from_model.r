#' Simulate from the turnover_simple model
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

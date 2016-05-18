#' Simulate from the turnover_simple model
#'
#' Works either for prior or posterior predictive stuff.
#'
#' @param ntax int number of taxa
#' @param ntime int number of time slices
#' @param pred matrix of predictors for N x T
#' @param p vector length T of preservation probabilities
#' @return list of z (true), y (observed)
#' @export
model.simulation <- function(ntax, ntime, pred, p = NULL) {
  if(!is.null(p)) p <- as.numeric(p)
  
  y <- z <- matrix(NA, nrow = ntax, ncol = ntime)

  for(nn in 1:ntax) {
    z[nn, 1] <- rbinom(1, 1, prob = pred[nn, 1])
    if(!is.null(p)) y[nn, 1] <- rbinom(1, 1, prob = z[nn, 1] * p[1])
    prod.term <- 1 - z[nn, 1]
    for(tt in 2:ntime) {
      prod.term <- prod.term * (1 - z[nn, tt - 1])
      z[nn, tt] <- rbinom(1, 1, prob = z[nn, tt - 1] * pred[nn, tt] + 
                          prod.term * pred[nn, tt])
      if(!is.null(p)) y[nn, tt] <- rbinom(1, 1, prob = z[nn, tt] * p[tt])
    }
  }
  if(is.null(p)) y <- NULL

  out <- list(z = z, y = y)
  out
}

#' Simulate from the turnover_simple model
#'
#' Works either for prior or posterior predictive stuff.
#'
#' @param ntax int number of taxa
#' @param ntime int number of time slices
#' @param p.mu scalar mean of sampling
#' @param p.sigma scalar st dev of sampling
#' @param p.eff vector time effect of sampling 
#' @param inter.mu scalar mean of presence
#' @param inter.sigma scalar st dev of presence
#' @param inter.eff vector time effect of presence
#' @return list of z, y, pred, and p
#' @export
model.simulation <- function(ntax, ntime, 
                             p.mu, p.sigma, p.eff, 
                             inter.mu, inter.sigma, inter.eff) {

  y <- z <- pred <- matrix(NA, nrow = ntax, ncol = ntime)
  p <- c()

  for(nn in 1:ntax) {
    p[1] <- invlogit(p.mu + p.sigma * p.eff[1])
    pred[nn, 1] <- invlogit(inter.mu + inter.sigma * inter.eff[1])
    z[nn, 1] <- rbinom(1, 1, prob = pred[nn, 1])
    y[nn, 1] <- rbinom(1, 1, prob = z[nn, 1] * p[1])
    for(tt in 2:ntime) {
      p[tt] <- invlogit(p.mu + p.sigma * p.eff[tt])
      pred[nn, tt] <- invlogit(inter.mu + inter.sigma * inter.eff[tt])

      z[nn, tt] <- rbinom(1, 1, prob = z[nn, tt - 1] * pred[nn, tt] + 
                          (prod(1 - z[nn, 1:tt - 1]) * pred[nn, tt]))
      y[nn, tt] <- rbinom(1, 1, prob = z[nn, tt] * p[tt])
    }
  }

  out <- list(z = z, y = y, pred = pred, p = p)
  out
}

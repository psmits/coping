# simulate from the posterior predictive distribution
sim.series <- function(data, ext) {
  out <- list()
  ny <- table(data$year)
  for(ii in seq(data$T)) {
    ll <- sample(ext$loc[, ii], 1) 
    # ext$phy
    #   grab those that originated at that time
    ss <- sample(ext$scale[, ii], 1)
    aa <- sample(ext$skew[, ii], 1)

    bt <- ext$beta_cat[sample(nrow(ext$beta_cat), 1), ]

    hold <- c()
    subdata <- data$cat[data$year == ii, ]
    for(jj in seq(nrow(subdata))) {
      # what taxa are from this year?
      # this is all fucked up
      sim <- rsn(1, xi = ll + subdata[jj, ] %*% bt,
                 omega = ss, 
                 alpha = aa)
      hold[jj] <- sim
    }
    out[[ii]] <- hold
  }
  out
}


# expectation
skew.exp <- function(mu, sd, alpha) {  # expectation
  delt <- alpha / sqrt(1 + (alpha ^ 2))
  mu + (sd * delt * sqrt(2 / pi))
}
skew.wrap <- function(attr) {
  skew.exp(mu = attr[1], sd = attr[2], alpha = attr[3])
}


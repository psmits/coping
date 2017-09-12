sort.climate <- function(x, val = 'mean', bin = '2My', nalma = NULL) {

  temp.est <- x$Temperature
  temp.range <- x$Temperature.max - x$Temperature.min

  temp.time.mean <- list()
  temp.time.range <- list()
  if(bin == '2My') {
    b <- range(occur$bins)
    b <- seq(b[1], b[2], by = 2)
    b <- as.matrix(cbind(b - 2, b))
    for(ii in seq(nrow(b))) {
      temp.time.mean[[ii]] <- temp.est[x$Age >= b[ii, 1] & x$Age < b[ii, 2]]
      temp.time.range[[ii]] <- temp.range[x$Age >= b[ii, 1] & x$Age < b[ii, 2]]
    }
  } else if (bin == 'NALMA') {
		b <- cbind(c(0, nalma$ma[-nrow(nalma)]), nalma$ma)
    for(ii in seq(nrow(b))) {
      temp.time.mean[[ii]] <- temp.est[x$Age >= b[ii, 1] & x$Age < b[ii, 2]]
      temp.time.range[[ii]] <- temp.range[x$Age >= b[ii, 1] & x$Age < b[ii, 2]]
    }
  }
  if(val == 'mean') {
    out <- arm::rescale(laply(temp.time.mean, function(x) 
                              mean(x, na.rm = TRUE)))
  } else if (val == 'range') {
    out <- arm::rescale(laply(temp.time.range, function(x) 
                              mean(x, na.rm = TRUE)))
  }
  out
}

library(plyr)
library(coda)
library(arm)
library(stringr)
library(reshape2)
library(ggplot2)
library(rstan)
library(grid)
library(pROC)
library(Metrics)
source('../R/multiclass_roc.r')
source('../R/read_one_stan.r')
source('../R/sim_from_model.r')
source('../data/data_dump/trait_info.data.R')
sight.implied <- sight
source('../data/data_dump/trait_w_gaps.data.R')
sight.obs <- sight



post.advi <- function(fit) {
  name.situation <- str_split(colnames(fit), '\\.')
  which.param <- laply(name.situation, function(x) x[1])
  uni.param <- unique(which.param)

  out <- list()
  for(ii in seq(length(uni.param))) {
    matching <- which.param == uni.param[ii]
    name.m <- name.situation[matching]

    te <- length(name.m[[1]])
    if(te == 1) {
      out[[ii]] <- fit[, matching]
    } else if (te == 2) {
      tt <- fit[, matching]
      numbers <- laply(name.m, function(x) x[-1])
      dims <- max(as.numeric(numbers))
      hh <- array(NA, dim = c(1001, dims)) 
      for(jj in seq(length(numbers))) {
        posit <- as.numeric(numbers[jj])
        hh[, posit] <- tt[, jj]
      }
      out[[ii]] <- hh
    } else if (te == 3) {
      tt <- fit[, matching]
      numbers <- laply(name.m, function(x) x[-1])
      dims <- apply(numbers, 2, function(x) max(as.numeric(x)))
      hh <- array(NA, dim = c(1001, dims))
      for(jj in seq(nrow(numbers))) {
        posit <- as.numeric(numbers[jj, ])
        hh[, posit[1], posit[2]] <- tt[, jj]
      }
      out[[ii]] <- hh
    }
  }
  names(out) <- uni.param
  out
}

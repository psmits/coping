# multiclass roc
library(pROC)
library(Metrics)

allvone <- function(pred.res, obs) {
  #pred.res <- roc.setup
  #obs <- paste0('X', y[cohort == 27])
  pred <- factor(pred.res[, 1], levels = colnames(pred.res[, -1]))
  prob <- lapply(levels(pred), function(class) {
                 pp <- ifelse(pred == class, 1, 0)
                 oo <- ifelse(obs == class, 1, 0)
                 prob <- pred.res[, class]

                 ps <- Metrics::auc(oo, prob)
                 ps})
  prob <- prob[!is.nan(unlist(prob))]
  mauc <- base::colMeans(do.call(rbind, prob), na.rm = TRUE)
  mauc
}



roc.dist <- function(ext, y, cohort) {
  n <- nrow(ext$hold)
  oo <- list()
  for(ii in seq(C)) {
    oo[[ii]] <- ext$hold[sample(n, 1), cohort == ii, ]
  }
  bad <- laply(llply(oo, dim), is.null)
  oo <- oo[!bad]
  hold.prob <- llply(oo, function(x) t(apply(x, 1, softmax))) # 

  roc.hold <- c()
  for(ii in seq(length(hold.prob))) {
    a <- ii
    if(ii >= which(bad)) a <- a + 1

    roc.setup <- data.frame(pred = paste0('X', apply(hold.prob[[ii]], 
                                                     1, which.max)), 
                            hold.prob[[ii]])
    roc.hold[ii] <- allvone(roc.setup, obs = paste0('X', y[cohort == a]))
  }
  roc.hold
}

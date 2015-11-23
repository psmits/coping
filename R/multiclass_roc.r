# this code is based on other peoples work.
#
# namely modern tool making blogspot
#
# also based on caret::twoClassSummary
multiClassSummary <- function(data, lev = NULL, model = NULL) {
  require(pROC)
  require(Metrics)
  require(caret)

  #check data
  if(!all(levels(data[, 'pred']) == levels(data[, 'obs']))) {
    stop('levels of observed and predicted data do not match')
  }

  # one-vs-all stats for each class
  prob.stats <- lapply(levels(data[, 'pred']), function(class) {

                       # one-vs-all set up
                       pred <- ifelse(data[, 'pred'] == class, 1, 0)
                       obs <- ifelse(data[, 'obs'] == class, 1, 0)
                       prob <- data[, class]

                       # calc one-vs-all AUC, logLoss
                       cap.prob <- pmin(pmax(prob, 0.000001), 0.999999)
                       prob.stats <- c(auc(obs, prob), logLoss(obs, cap.prob))
                       names(prob.stats) <- c('ROC', 'logLoss')
                       return(prob.stats)})
  
  prob.stats <- do.call(rbind, prob.stats)
  rownames(prob.stats) <- paste('Class:', levels(data[, 'pred']))

  # confusion stats
  CM <- confusionMatrix(data[, 'pred'], data[, 'obs'])

  # aggregate
  class.stats <- cbind(CM$byClass, prob.stats)
  class.stats <- colMeans(class.stats)

  overall.stats <- c(CM$overall)

  # combine
  stats <- c(overall.stats, class.stats)
  stats <- stats[! names(stats) %in% c('AccuracyNull', 
                                       'Prevelance', 
                                       'Detection Prevelance')]
  names(stats) <- gsub('[[:blank:]]+', '_', names(stats))
  return(stats)

}


# multiclass roc
allvone <- function(pred.res, obs) {
  pred <- factor(pred.res[, 1])
  prob <- lapply(levels(pred), function(class) {
                 pp <- ifelse(pred == class, 1, 0)
                 oo <- ifelse(obs == class, 1, 0)
                 prob <- pred.res[, class]

                 ps <- auc(oo, prob)
                 ps})
  mauc <- colMeans(do.call(rbind, prob))
  mauc
}

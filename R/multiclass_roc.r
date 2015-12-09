# multiclass roc
library(pROC)
library(Metrics)

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

#' Make the plots for the coping project
#'
#' @param ext1 output from rstan::extract of model fit
#' @param name string associated with outputs
#' @param name.name vector of regression coefficient names
#' @param group logical are there any group level covariates?
#' @param sampling logical is sampling parameterized?
#' @param nsim integer how many random pulls from the posterior?
make.plots <- function(ext1) {

  # log-odds of occurrence associated with ecotype
  #   controlling for mass
  #   given group-level effects
  am <- melt(ext1$a)  # sim, time, state, value
  
  amplot <- ggplot(am, aes(x = Var2, y = value, group = Var1))
  amplot <- amplot + geom_line(alpha = 0.1)
  amplot <- amplot + facet_wrap(~ Var3)
  amplot <- amplot + labs(x = 'Time (Mya)', 
                          y = 'log-odds occurrence')


  # effect of group-level on individual-level
  gm <- melt(ext1$gamma)  # sim, pred, ecotype, value

  gmplot <- ggplot(gm, aes(x = Var2, y = value, group = Var2))
  gmplot <- gmplot + geom_boxplot()
  gmplot <- gmplot + facet_wrap(~ Var3)
  gmplot <- gmplot + labs(x = 'Predictor variable', 
                          y = 'Effect on log-odds of occurrence')


  # variation in log-odds of occurrence associated with ecotype
  tm <- melt(ext1$tau)  # sample, ecotype, value

  tmplot <- ggplot(tm, aes(x = value, y = ..density..))
  tmplot <- tmplot + geom_histogram()
  tmplot <- tmplot + facet_wrap(~ Var2)
  tmplot <- tmplot + labs(x = 'Estimated standard deviation of distribution of ecotype log-odds of occurrence', 
                          y = 'Probability density')

  ext1$alpha_0
  apply(ext1$alpha_time, 2, function(x) x + ext1$alpha_0)

  # difference from mean log-odds observation due to time
  pm <- melt(ext1$alpha_time)
  
  pmplot <- ggplot(pm, aes(x = Var2, y = value, group = Var1))
  pmplot <- pmplot + geom_line(alpha = 0.1)
  pmplot <- pmplot + labs(x = 'Time (Mya)', 
                          y = 'Difference from mean log-odds observation')

}

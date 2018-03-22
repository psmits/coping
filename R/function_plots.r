# transform occurrence matrix into list of functional distance matrices
#   functional distance is based on two categorical variables
#     diet
#     locomotion
#   gower dissimilarity

# columns are time
# rows are species

# what i do is aggregate by functional group
# each functional group has a count
# then as average functional diversity
# functional diversity of each 
# does this even make sense?
qs <- seq(0, 4, by = 0.1)
qq_post <- list()
for(jj in seq(nsim)) {
  qq <- list()
 
  pp <- data.frame(post.div[[jj]], ecotype)
  pp <- pp %>% 
    group_by(X1.1, X2.1) %>% 
    dplyr::summarise_if(is.numeric, sum) %>% 
    ungroup()
  dd <- daisy(pp[, 1:2], metric = 'gower')
  
  # calculate effective diversity for each of the units
  oo <- purrr::map(qs, function(x) Func2014(dd, pp[, -(1:2)], q = x))
  fd <- purrr::map(oo, ~ .x$FuncD)
  mfd <- purrr::map(oo, ~ .x$Q * .x$FuncD)
  tfd <- purrr::map2(oo, mfd, ~ .x$FuncD * .y)

  # calculate effective similarity between units
  lst_cqn <- lst_uqn <- list()
  for(nn in seq(length(qs))) {
    mat_cqn <- mat_uqn <- matrix(NA, ncol = T, nrow = T)
    for(ii in seq(T)) {
      for(kk in seq(T)) {
        hol <- Func2014(dd, pp[, c(ii + 2, kk + 2)], q = qs[nn])
        mat_cqn[ii, kk] <- hol[[3]]['FunCqN']
        mat_uqn[ii, kk] <- hol[[3]]['FunUqN']
      }
    }
    lst_cqn[[nn]] <- mat_cqn
    lst_uqn[[nn]] <- mat_uqn
  }
  names(fd) <- names(mfd) <- names(tfd) <- names(lst_cqn) <- names(lst_uqn) <- 
    paste0('q', qs)
  
  qq_post[[jj]] <- list(fd = fd,
                        mfd = mfd,
                        tfd = tfd,
                        fcqn = lst_cqn,
                        fuqn = lst_uqn)
}

get_diversity <- function(qq_post, div = 'fd') {
  funcd <- qq_post %>% purrr::map(~ purrr::transpose(.x[[div]])) %>%
    purrr::map(., ~ purrr::map(.x, flatten_dbl)) %>%
    purrr::map(., function(x) {
                 x <- purrr::map(x, ~ cbind(fd = .x, q = qs))
                 x}) %>%
    purrr::map(., function(x) {
                 x <- purrr::map2(x, seq(T), ~ cbind(.x, t = .y))
                 x}) %>%
    purrr::map(., ~ reduce(.x, rbind)) %>%
    purrr::map2(., seq(length(.)), function(x, y) {
                  x <- cbind(x, sim = y)
                  x}) %>%
    purrr::reduce(., rbind) %>%
    data.frame(.)
  funcd
}

fd <- get_diversity(qq_post, div = 'fd')
mfd <- get_diversity(qq_post, div = 'mfd')
tfd <- get_diversity(qq_post, div = 'tfd')
# comparison of hill numbers for effective diversity of functional group
# partial function evaluation
pm <- partial(merge, by = c('q', 't', 'sim'))
funcd <- reduce(list(fd, mfd, tfd), pm)
names(funcd) <- c('q', 't', 'sim', 'fd', 'mfd', 'tfd')

# plot effective functional diversity
funcd <- funcd %>%
  group_by(q, t) %>%
  dplyr::summarize(fd_med = median(fd),
                   fd_low = quantile(fd, 0.25),
                   fd_hgh = quantile(fd, 0.75),
                   mfd_med = median(mfd),
                   mfd_low = quantile(mfd, 0.25),
                   mfd_hgh = quantile(mfd, 0.75),
                   tfd_med = median(tfd),
                   tfd_low = quantile(tfd, 0.25),
                   tfd_hgh = quantile(tfd, 0.75)) %>%
  ungroup()

funcd$t <- factor(funcd$t, levels = seq(T))
fdg <- ggplot(funcd, aes(x = q, y = fd_med, group = t))
fdg <- fdg + geom_ribbon(mapping = aes(ymax = fd_hgh, ymin = fd_low, fill = t),
                         alpha = 0.2)
fdg <- fdg + geom_line(alpha = 1, mapping = aes(colour = t))





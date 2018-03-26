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
shannon <- which(qs == 1)
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
  lst_beta <- lst_cqn <- lst_uqn <- list()
  for(nn in seq(length(qs))) {
    mat_beta <- mat_cqn <- mat_uqn <- matrix(NA, ncol = T, nrow = T)
    for(ii in seq(T)) {
      for(kk in seq(T)) {
        hol <- Func2014(dd, pp[, c(ii + 2, kk + 2)], q = qs[nn])
        mat_cqn[ii, kk] <- hol[[3]]['FunCqN']
        mat_uqn[ii, kk] <- hol[[3]]['FunUqN']
        mat_beta[ii, kk] <- hol[[3]]['Beta']
      }
    }
    lst_cqn[[nn]] <- mat_cqn
    lst_uqn[[nn]] <- mat_uqn
    lst_beta[[nn]] <- mat_beta
  }
  names(fd) <- names(mfd) <- names(tfd) <- names(lst_cqn) <- names(lst_uqn) <- 
    paste0('q', qs)
  
  qq_post[[jj]] <- list(fd = fd,
                        mfd = mfd,
                        tfd = tfd,
                        betadiv = lst_beta,
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
funcg <- funcd %>%
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

funcg$t <- factor(funcg$t, levels = seq(T))
fdg <- ggplot(funcg, aes(x = q, y = fd_med, group = t))
fdg <- fdg + geom_ribbon(mapping = aes(ymax = fd_hgh, ymin = fd_low, fill = t),
                         alpha = 0.2)
fdg <- fdg + geom_line(alpha = 1, mapping = aes(colour = t))

mfdg <- ggplot(funcg, aes(x = q, y = mfd_med, group = t))
mfdg <- mfdg + geom_ribbon(mapping = aes(ymax = mfd_hgh, ymin = mfd_low, fill = t),
                           alpha = 0.2)
mfdg <- mfdg + geom_line(alpha = 1, mapping = aes(colour = t))

tfdg <- ggplot(funcg, aes(x = q, y = tfd_med, group = t))
tfdg <- tfdg + geom_ribbon(mapping = aes(ymax = tfd_hgh, ymin = tfd_low, fill = t),
                           alpha = 0.2)
tfdg <- tfdg + geom_line(alpha = 1, mapping = aes(colour = t))


# beta diverstiy
bd <- cmdscale(qq_post[[1]]$betadiv[[shannon]], k = D - 1, eig = TRUE)

# diff from t to t+1
shannon_beta <- purrr::map(qq_post, function(x) x$betadiv[[shannon]])

tt <- as.list(1:17)
tm <- purrr::map(tt, ~ .x + 1)
tm <- purrr::map2(tt, tm, ~ c(.x, .y))

beta_step <- purrr::map(shannon_beta, function(x) { 
                          purrr::map_dbl(tm, function(y) {
                                           x[y[1], y[2]]})})
beta_step <- reduce(beta_step, rbind)
colnames(beta_step) <- fct_drop(interaction(1:17, 2:18))

bs <- melt(beta_step)
bs$Var2 <- as.factor(bs$Var2)
ggplot(bs, aes(x = Var2, y = value, group = Var2)) + 
  geom_violin() +
  geom_jitter(width = .25, height = 0)

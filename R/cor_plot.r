library(corrplot)
library(igraph)

oo <- ss <- list()
for(ii in seq(dim(ext2$o_Omega)[1])) {
  oo[[ii]] <- ext2$o_Omega[ii, , ]
  ss[[ii]] <- ext2$s_Omega[ii, , ]
}
ocorr.mean <- Reduce('+', oo) / length(oo)
scorr.mean <- Reduce('+', ss) / length(ss)

ecotran <- ecotrans
ecotran[, 1] <- mapvalues(ecotran[, 1], 
                          from = unique(ecotran[, 1]), 
                          to = c('carnivore', 'herbivore', 
                                 'insectivore', 'omnivore'))
ecotran <- ecotran[, 2:1]
rownames(ocorr.mean) <- colnames(ocorr.mean) <- 
  apply(ecotran, 1, paste0, collapse = ' ')
col1 <- colorRampPalette(c('red', 'white', 'blue'))

png(file = '../doc/figure/origination_correlation.png',
    width = 1500, height = 1500)
corrplot::corrplot.mixed(ocorr.mean, 
                         lower = 'ellipse', upper = 'number', 
                         diag = 'n',
                         #col = col1(200),
                         tl.pos = 'lt',
                         tl.col = 'black',
                         tl.cex = 2,
                         cl.cex = 2,
                         number.cex = 2,
                         number.digits = 2,
                         mar = rep(0, 4))
dev.off()

rownames(scorr.mean) <- colnames(scorr.mean) <- 
  apply(ecotran, 1, paste0, collapse = ' ')

png(file = '../doc/figure/survival_correlation.png',
    width = 1500, height = 1500)
corrplot::corrplot.mixed(scorr.mean, 
                         lower = 'ellipse', upper = 'number', 
                         diag = 'n',
                         #col = col1(200),
                         tl.pos = 'lt',
                         tl.col = 'black',
                         tl.cex = 2,
                         cl.cex = 2,
                         number.cex = 2,
                         number.digits = 2)
dev.off()


# names and numbers of correlation in ecotype origination
coro.test <- Reduce('+', llply(oo, function(x) x > 0)) / length(oo)
coro.number <- coro.test
coro.test <- melt(coro.test > 0.85)
coro.test <- coro.test[coro.test[, 3] == TRUE, 1:2]
coro.test <- coro.test[coro.test[, 1] > coro.test[, 2], ]
if(nrow(coro.test) != 0) {
  coro.row <- apply(ecotran[coro.test[, 1], ], 1, paste0, collapse = ' ')
  coro.col <- apply(ecotran[coro.test[, 2], ], 1, paste0, collapse = ' ')
  coro.frame <- data.frame(coro.row, coro.col)

  cn <- c()
  for(ii in seq(nrow(coro.test))) {
    cn[ii] <- ocorr.mean[coro.test[ii, 1], coro.test[ii, 2]]
  }
  cn <- round(cn, 2)
  # make a graph of what's correlated with what
  coro.graph <- graph_from_edgelist(as.matrix(coro.frame), directed = FALSE)
  coord <- layout_in_circle(coro.graph)
  png(file = '../doc/figure/origin_sig_corr.png',
      width = 1500, height = 1500)
  plot(coro.graph, layout = coord, 
       vertex.shape = 'none', 
       vertex.label.cex = 2.5, 
       vertex.label.color = 'black',
       edge.width = 1.5,
       edge.color = 'darkgrey',
       edge.label = cn,
       edge.label.color = 'blue',
       edge.label.cex = 2,
       margin = rep(0, 4))
  dev.off()
}


# names and numbers of correlation in ecotype survival
cors.test <- Reduce('+', llply(ss, function(x) x > 0)) / length(oo)
cors.number <- cors.test
cors.test <- melt(cors.test > 0.85)
cors.test <- cors.test[cors.test[, 3] == TRUE, 1:2]
cors.test <- cors.test[cors.test[, 1] > cors.test[, 2], ]
if(nrow(cors.test) != 0) {
  #cors.row <- apply(ecotran[cors.test[, 1], ], 1, paste0, collapse = ' ')
  #cors.col <- apply(ecotran[cors.test[, 2], ], 1, paste0, collapse = ' ')
  cors.row <- paste0(ecotran[cors.test[, 1], ], collapse = ' ')
  cors.col <- paste0(ecotran[cors.test[, 2], ], collapse = ' ')
  cors.frame <- data.frame(cors.row, cors.col)

  cn <- c()
  for(ii in seq(nrow(cors.test))) {
    cn[ii] <- scorr.mean[cors.test[ii, 1], cors.test[ii, 2]]
  }
  cn <- round(cn, 2)
  # make a graph of what's correlated with what
  cors.graph <- graph_from_edgelist(as.matrix(cors.frame), directed = FALSE)
  coord <- layout_in_circle(cors.graph)
  png(file = '../doc/figure/surv_sig_corr.png',
      width = 500, height = 500)
  plot(cors.graph, layout = coord, 
       vertex.shape = 'none', 
       vertex.label.cex = 1.5, 
       vertex.label.color = 'black',
       edge.width = 1.5,
       edge.color = 'darkgrey',
       edge.label = cn,
       edge.label.color = 'blue',
       edge.label.cex = 1.5,
       margin = rep(0, 4))
  dev.off()
}

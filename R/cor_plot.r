library(ellipse)
library(corrplot)
corrplot(ext2$o_Omega[1, , ], method = 'number')
corrplot(ext2$s_Omega[1, , ], method = 'number')

oo <- ss <- list()
for(ii in seq(dim(ext2$o_Omega)[1])) {
  oo[[ii]] <- ext2$o_Omega[ii, , ]
  ss[[ii]] <- ext2$s_Omega[ii, , ]
}
ocorr.mean <- Reduce('+', oo) / length(oo)
scorr.mean <- Reduce('+', ss) / length(ss)

rownames(ocorr.mean) <- apply(ecotrans, 1, paste0, collapse = ' ')
colnames(ocorr.mean) <- apply(ecotrans, 1, paste0, collapse = ' ')
corrplot(ocorr.mean, method = 'number')

rownames(scorr.mean) <- apply(ecotrans, 1, paste0, collapse = ' ')
colnames(scorr.mean) <- apply(ecotrans, 1, paste0, collapse = ' ')
corrplot(scorr.mean, method = 'number')


# names and numbers of correlation in ecotype origination
coro.test <- Reduce('+', llply(oo, function(x) x > 0)) / length(oo)
coro.number <- coro.test
coro.test <- melt(coro.test > 0.95)
coro.test <- coro.test[coro.test[, 3] == TRUE, 1:2]
coro.test <- coro.test[coro.test[, 1] > coro.test[, 2], ]
coro.row <- apply(ecotrans[coro.test[, 1], ], 1, paste0, collapse = ' ')
coro.col <- apply(ecotrans[coro.test[, 2], ], 1, paste0, collapse = ' ')
coro.frame <- data.frame(coro.row, coro.col)

# names and numbers of correlation in ecotype survival
cors.test <- Reduce('+', llply(ss, function(x) x > 0)) / length(oo)
cors.number <- cors.test
cors.test <- melt(cors.test > 0.95)
cors.test <- cors.test[cors.test[, 3] == TRUE, 1:2]
cors.test <- cors.test[cors.test[, 1] > cors.test[, 2], ]
#cors.row <- apply(ecotrans[cors.test[, 1], ], 1, paste0, collapse = ' ')
#cors.col <- apply(ecotrans[cors.test[, 2], ], 1, paste0, collapse = ' ')
cors.row <- paste0(ecotrans[cors.test[, 1], ], collapse = ' ')
cors.col <- paste0(ecotrans[cors.test[, 2], ], collapse = ' ')
cors.frame <- data.frame(cors.row, cors.col)

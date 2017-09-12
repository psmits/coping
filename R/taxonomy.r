update.taxonomy.eol <- function(occur, key) {
  # get away from the input asap
  tt <- occur
  #tt <- occur[, c('order', 'family', 'genus')]
  
  # these parts are slow, use uniques and then put back in
  gen <- unique(tt[, 3])
  #eol.id <- llply(gen, function(x) try(get_eolid(x, key = key, rows = 1)))
  eol.id <- mclapply(gen, function(x) try(get_eolid(x, key = key, rows = 1)), 
                     mc.cores = detectCores())
    # there can be strange errors surrounding some taxa
    # by using try, i at least always get output
  #eol.class <- llply(eol.id, function(x) try(classification(x, key = key)))
  eol.class <- mclapply(eol.id, function(x) try(classification(x, key = key)), 
                        mc.cores = detectCores())

  # get the order, family information from those taxa found
  out <- llply(eol.class, function(x) {
                 if(class(x[[1]]) == 'data.frame') {
                   # make sure it is a mammal
                   if(any(x[[1]][, 1] == 'Mammalia')) {
                     a <- x[[1]]$rank %in% c('order', 'family')
                     b <- x[[1]][a, ]
                     if(nrow(b) == 1) {
                       c(NA, b)
                     } else {
                       b[, 1]
                     }
                   } else {
                     c(NA, NA)
                   }
                 } else {  # bad responses incl. errors
                   c(NA, NA)
                 }})
  out <- Reduce(rbind, out)
  out <- cbind(out, gen)
 
  # everything is good till about here
  # sometimes i haven't learned more
  no.new <- apply(out[, 1:2], 2, is.na)
  for(ii in seq(nrow(out))) {
    if(!any(is.na(out[ii, ]))) {
      upd <- tt[, 3] %in% out[ii, 3]
      tt[upd, 1:2] <- out[ii, 1:2]
    }
  }

  tt
}

# not everything is on EOL
# especially extinct things
# two stage process
#   put genera in families
#   put families in orders
update.taxonomy.extinct <- function(occur, extinct) {
  # get away from the input asap
  tt <- occur
  #tt <- occur[, c('order', 'family', 'genus')]
 
  # genera in families
  ext.match <- llply(extinct[[1]], function(x) tt[, 3] %in% x)
  
  for(ii in seq(length(ext.match))) {
    if(any(ext.match[[ii]])) {
      tt[ext.match[[ii]], 2] <- names(extinct[[1]][ii])
    }
  }


  # families in orders
  # if the family is extinct and i know it, give it an order
  ext.match <- llply(extinct[[2]], function(x) tt[, 2] %in% x)
  
  for(ii in seq(length(ext.match))) {
    if(any(ext.match[[ii]])) {
      tt[ext.match[[ii]], 1] <- names(extinct[[2]][ii])
    }
  }

  tt
}



# i need to consolidate orders
#   eg Lypotiphila == Eulypotiphila
#

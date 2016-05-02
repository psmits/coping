# process the "climate" information
zac <- read.csv('../data/2008_zachos_data.csv', stringsAsFactors = FALSE)

names(zac) <- c('loc', 'age', 'genus', 'o18', 'c13', 
                'o18.5pt', 'c13.5pt', 'comment')
zac <- zac[zac$age <= (max(occur$bins)), ]
zac <- zac[!is.na(zac$o18), ]
b <- unique(occur$bins)
b <- as.matrix(cbind(b - 2, b))
isotope <- list()
zac.cohort <- array(NA, dim = nrow(zac))
for(ii in seq(nrow(b))) {
  isotope[[ii]] <- zac$o18[zac$age > b[ii, 1] & zac$age <= b[ii, 2]]
  zac.cohort[(zac$age > b[ii, 1] & zac$age <= b[ii, 2])] <- ii
}
zac <- zac[!is.na(zac.cohort), ]
zac.cohort <- zac.cohort[!is.na(zac.cohort)]
mean.o18 <- laply(split(zac$o18, zac.cohort), mean)
mean.o18 <- arm::rescale(mean.o18)
range.o18 <- laply(split(zac$o18, zac.cohort), function(x) 
                   abs(quantile(x, 0.25) - quantile(x, 0.75)))
range.o18 <- arm::rescale(range.o18)


# how about the mg/ca data...
cram <- read.delim('../data/cramer/cramer_mgca_pref.txt', sep = '\t')
cram.temp <- read.delim('../data/cramer/cramer_temp.txt', sep = '\t')
temp.est <- cram.temp$Temperature
temp.range <- cram.temp$Temperature.max - cram.temp$Temperature.min

b <- unique(occur$bins)
b <- as.matrix(cbind(b - 2, b))
temp.time.mean <- list()
temp.time.range <- list()
for(ii in seq(nrow(b))) {
  temp.time.mean[[ii]] <- temp.est[cram.temp$Age > b[ii, 1] & 
                                   cram.temp$Age <= b[ii, 2]]
  temp.time.range[[ii]] <- temp.range[cram.temp$Age > b[ii, 1] & 
                                      cram.temp$Age <= b[ii, 2]]
}
#temp.time.mean <- laply(temp.time.mean, function(x) mean(x, na.rm = TRUE))
temp.time.mean <- arm::rescale(laply(temp.time.mean, function(x) 
                                     mean(x, na.rm = TRUE)))
temp.time.range <- arm::rescale(laply(temp.time.range, function(x) 
                                      mean(x, na.rm = TRUE)))
#plot(cram.temp$Age, cram.temp$Temperature)
#points(seq(length(temp.time.mean)) * 2, temp.time.mean, col = 'blue')

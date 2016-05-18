library(plyr)
library(reshape2)
library(stringr)

binom.make <- function(genus, species, sep = ' ') {
  paste(genus, species, sep = sep)
}

clean.occurrence <- function(dat) {
  # remove specimens that don't have time assigned
  dat <- dat[!is.na(dat$max_ma), ]
  dat <- dat[!is.na(dat$min_ma), ]
  dat$ma_mid <- apply(cbind(dat$max_ma, dat$min_ma), 1, mean)

  # make sure i'm entirely in the Cenozoic via max
  kpg <- 65.6
  dat <- dat[!dat$ma_mid > kpg, ]

  dat$name.bi <- str_replace(dat$accepted_name, 
                             pattern = ' ', replacement = '_')

  # exclude aquatic and volant nonsense
  aq <- c('Cetacea', 'Desmostylia', 'Sirenia', 'Chiroptera')
  seal <- c('Odobenidae', 'Otariidae', 'Phocidae', 'Desmatophocidae')
  badgen <- c('Enaliarctos', 'Pacificotaria', 'Pinnarctidion', 'Pteronarctos', 'Wallia')
  dat <- dat[!(dat$order %in% aq), ]
  dat <- dat[!(dat$family %in% seal), ]
  dat <- dat[!(dat$genus %in% badgen), ]
  lf <- c('amphibious', 'volant', 'aquatic', 'gliding')
  dat <- dat[!(dat$life_habit %in% lf), ]

  # diet assignments
  herb <- c('herbivore', 'grazer', 'browser', 'folivore', 'granivore')
  omm <- c('frugivore', 'omnivore')
  car <- c('carnivore')#, 'insectivore')
  insect <- c('insectivore')
  dat$diet1 <- str_extract(dat$diet, pattern = '[^,]*')
  dat$comdiet <- dat$diet1
  dat$comdiet[dat$diet1 %in% herb] <- 'herb'
  dat$comdiet[dat$diet1 %in% omm] <- 'omni'
  dat$comdiet[dat$diet1 %in% car] <- 'carni'
  dat$comdiet[dat$diet1 %in% insect] <- 'insect'

  # locomotor assignments
  tree <- c('arboreal')
  ground <- c('ground dwelling', 'saltatorial')
  dig <- c('semifossorial', 'fossorial')
  dat$comlife <- dat$life_habit
  dat$comlife[dat$life_habit %in% tree] <- 'arboreal'
  dat$comlife[dat$life_habit %in% ground] <- 'ground dwelling'
  dat$comlife[dat$life_habit %in% dig] <- 'fossorial'

  # assign every occurence to a 2 My bin
  bins <- seq(from = 0, to = 66, by = 2)
  bins <- cbind(top = bins[-1], bot = bins[-length(bins)])
  dat$bins <- rep(NA, nrow(dat))
  for (ii in seq(nrow(bins))) {
    out <- which(dat$ma_mid < bins[ii, 1] & dat$ma_mid >= bins[ii, 2])
    dat$bins[out] <- bins[ii, 1]
  }

  dat
}

library(plyr)
library(reshape2)
library(stringr)
library(rstan)
library(arm)
library(ape)
library(geiger)
library(caret)

source('../R/mung.r')

load('../data/scaled_super.rdata')
load('../data/body_mass.rdata')

dat <- read.csv('../data/mam-occs.csv', stringsAsFactors = FALSE)
posture <- read.csv('../data/posture.csv', stringsAsFactors = FALSE) 
# these specific assignments are based on Carano's papers on posture

occur <- clean.occurrence(dat)
ss <- split(occur, occur$bins)
occur <- Reduce(rbind, llply(ss, function(x) x[duplicated(x$name.bi), ]))

load('../data/update_taxonomy.rdata')
na.tax <- na.tax[!(duplicated(na.tax$name.bi)), ]
fx <- na.tax$name.bi %in% occur$name.bi[occur$comlife == 'ground dwelling']
na.fx <- na.tax[fx, ]
occur[match(na.fx$name.bi, occur$name.bi), 
      c('order_name', 'family_name')] <- na.fx[, 1:2]

stance.group <- split(posture$taxon, posture$stance)
for(ii in seq(length(stance.group))) {
  mm <- occur$family_name %in% stance.group[[ii]]
  gd <- occur$comlife == 'ground dwelling'
  occur$comlife[mm & gd] <- names(stance.group)[ii]

  mm <- occur$order_name %in% stance.group[[ii]]
  gd <- occur$comlife == 'ground dwelling'
  occur$comlife[mm & gd] <- names(stance.group)[ii] 
}

occur <- occur[occur$comlife != 'ground dwelling', ]



# shrink data to match all inputs
#   this would be the opportunity for setting up an imputation step
occur <- occur[occur$name.bi %in% na.mass$name, ]
occur$mass <- na.mass[match(occur$name.bi, na.mass$name), 2]

# make data and tree match
name <- matrix(str_replace(occur$name.bi, ' ', '_'), ncol = 1)
hot.fix <- name.check(spt, data.names = name)
occur <- occur[!(occur$name.bi %in% 
                 str_replace(hot.fix$data_not_tree, '_', ' ')), ]
spt <- ape::drop.tip(spt, hot.fix$tree_not_data)

## for testing purposes
#keep <- createDataPartition(occur$bins, p = 0.2)
#keepname <- str_replace(unique(occur$name.bi[keep[[1]]]), ' ', '_')
#hot.fix <- name.check(spt, data.names = keepname)
#spt <- ape::drop.tip(spt, hot.fix$tree_not_data)
#occur <- occur[keep[[1]], ]


# process climate information
source('../R/mung_clim.r')


# need to make the things work
by.tax <- split(occur, occur$occurrence.genus_name)

sight <- matrix(0, nrow = length(by.tax), ncol = length(unique(occur$bins)))
for(ii in seq(length(by.tax))) {
  sight[ii, (by.tax[[ii]]$bins / 2)] <- 1
}

N <- length(by.tax)
diet <- laply(by.tax, function(x) names(which.max(table(x$comdiet))))
life <- laply(by.tax, function(x) names(which.max(table(x$comlife))))
mass <- arm::rescale(log(laply(by.tax, function(x) mean(x$mass))))

# make covariates the right shape
diet <- model.matrix( ~ diet - 1)[, -1]
life <- model.matrix( ~ life - 1)[, -1]

# covariates
x <- cbind(mass, diet, life)  # for sep intercept set up
D <- ncol(x)


# temporal data
cohort <- occur$bins / 2
T <- length(unique(cohort))
cohort <- mapvalues(cohort, from = unique(cohort), to = seq(T))

# phylogenetic data
#vcv <- vcv(spt)
#vcv <- vcv / max(diag(vcv))
#U <- nrow(vcv)
#id <- as.numeric(as.factor(occur$name.bi))

# climate data
u <- cbind(mean.o18, range.o18)
U <- ncol(u)
### WARNING to include need to remove oldest bin!!!!
# u <- cbind(temp.time.mean, temp.time.range)
# U <- ncol(u)

# make the plant phase indicator
pa.eo <- c(66, 50)
eo.mi <- c(50, 16)
mi.pl <- c(16, 2)
mod <- c(2, 0)
plants <- rbind(mod, mi.pl, eo.mi, pa.eo)

phase <- array(dim = nrow(occur))
for(ii in seq(length(phase))) {
  mm <- c()
  for(jj in seq(nrow(plants))) {
    mm[jj] <- plants[jj, 1] >= occur$bins[ii] & plants[jj, 2] < occur$bins[ii]
  }
  phase[ii] <- which(mm)
}
P <- 4


# dump it out
stan_rdump(list = c('N', 'T', 'D', 'U', 'P',
                    'sight', 'x', 'u', 'phase'),
           file = '../data/data_dump/trait_info.data.R')

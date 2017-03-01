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

TESTING.x <- FALSE
TESTING.u <- FALSE



posture <- read.csv('../data/posture.csv', stringsAsFactors = FALSE) 
# these specific assignments are based on Carano's papers on posture

#dat <- read.csv('https://paleobiodb.org/data1.2/occs/list.csv?datainfo&rowcount&base_name=Mammalia&taxon_reso=species&interval=Maastrichtian,Gelasian&cc=NOA&show=class,genus,ecospace,loc,strat,stratext,lith,acconly', stringsAsFactors = FALSE, skip = 20)
dat <- read.csv('../data/pbdb_data.csv', stringsAsFactors = FALSE, skip = 21)

occur <- clean.occurrence(dat)
ss <- split(occur, occur$bins)
occur <- Reduce(rbind, llply(ss, function(x) x[duplicated(x$name.bi), ]))

load('../data/update_taxonomy.rdata')
na.tax <- na.tax[!(duplicated(na.tax$name.bi)), ]
fx <- na.tax$name.bi %in% occur$name.bi[occur$comlife == 'ground dwelling']
na.fx <- na.tax[fx, ]
occur[match(na.fx$name.bi, occur$name.bi), 
      c('order', 'family')] <- na.fx[, 1:2]

stance.group <- split(posture$taxon, posture$stance)
for(ii in seq(length(stance.group))) {
  mm <- occur$family %in% stance.group[[ii]]
  gd <- occur$comlife == 'ground dwelling'
  occur$comlife[mm & gd] <- names(stance.group)[ii]

  mm <- occur$order %in% stance.group[[ii]]
  gd <- occur$comlife == 'ground dwelling'
  occur$comlife[mm & gd] <- names(stance.group)[ii] 
}

occur <- occur[occur$comlife != 'ground dwelling', ]


# shrink data to match all inputs
#   this would be the opportunity for setting up an imputation step
na.mass$name <- str_replace(na.mass[, 1], ' ', '_')
occur <- occur[occur$name.bi %in% na.mass$name, ]
occur$mass <- na.mass[match(occur$name.bi, na.mass$name), 2]

## make data and tree match
#name <- matrix(str_replace(occur$name.bi, ' ', '_'), ncol = 1)
#hot.fix <- name.check(spt, data.names = name)
#occur <- occur[!(occur$name.bi %in% 
#                 str_replace(hot.fix$data_not_tree, '_', ' ')), ]
#spt <- ape::drop.tip(spt, hot.fix$tree_not_data)

## for testing purposes
#keep <- createDataPartition(occur$bins, p = 0.2)
#keepname <- str_replace(unique(occur$name.bi[keep[[1]]]), ' ', '_')
#hot.fix <- name.check(spt, data.names = keepname)
#spt <- ape::drop.tip(spt, hot.fix$tree_not_data)
#occur <- occur[keep[[1]], ]

# process climate information
source('../R/mung_clim.r')




occur <- occur[occur$bins != 66, ]

# save true bins
occur$true.bin <- occur$bins
# make easy bins
occur$bins <- occur$bins / 2


# !!! makes everything genus level !!!
# need to make the things work
by.tax <- split(occur, occur$name.bi)
#by.tax <- split(occur, occur$genus)

# cols go from younger to older
sight <- matrix(0, nrow = length(by.tax), ncol = max(occur$bins))
for(ii in seq(length(by.tax))) {
  sight[ii, (by.tax[[ii]]$bins)] <- 1
}
sight <- sight[, -1]


N <- length(by.tax)
diet <- laply(by.tax, function(x) names(which.max(table(x$comdiet))))
life <- laply(by.tax, function(x) names(which.max(table(x$comlife))))
mass <- arm::rescale(log(laply(by.tax, function(x) mean(x$mass))))

# covariates
diet <- factor(diet, levels = unique(diet)[c(4, 1:3)])
life <- factor(life, levels = unique(life))
inter <- interaction(diet, life, drop = TRUE)
inter <- as.character(inter)

rms <- which(inter == names(which(table(inter) == 1)))
sight <- sight[-rms, ]

x <- stats::model.matrix( ~ mass + inter)
x <- x[-rms, ]
mass <- mass[-rms]
#x <- x[, colSums(x) != 0]  
# i think this is the correct thing to do because certain inter not obs
D <- ncol(x)
name.name <- colnames(x)

# phylogenetic data
#vcv <- vcv(spt)
#vcv <- vcv / max(diag(vcv))
#U <- nrow(vcv)
#id <- as.numeric(as.factor(occur$name.bi))

# climate data
#u <- cbind(mean.o18, range.o18)
#U <- ncol(u)
# WARNING to include need to remove oldest bin!!!!
u <- cbind(temp.time.mean, temp.time.range)  # young to old
u <- u[-(nrow(u)), ]
u <- apply(u, 2, rev)  # old to young


# make the plant phase indicator
mi.pl <- c(16, 2)
eo.mi <- c(50, 16)
pa.eo <- c(66, 50)
plants <- rbind(mi.pl, eo.mi, pa.eo)

co.h <- seq(4, (ncol(sight) + 1) * 2, by = 2)
phase <- array(dim = length(co.h))
for(ii in seq(length(co.h))){
  phase[ii] <- which(plants[, 1] >= co.h[ii] & plants[, 2] < co.h[ii])
}
P <- 3

# fixed the reversed order
sight <- sight[, rev(seq(ncol(sight)))]  # from older to younger
phase <- factor(rev(phase)) # from older to younger
u <- model.matrix( ~ u + phase)
U <- ncol(u)

N <- nrow(sight)
T <- ncol(sight)

state <- as.numeric(factor(inter))
state <- state[-rms]
ufull <- u
u <- u[-1, ]
state <- mapvalues(state, 
                   from = sort(unique(state)), 
                   to = seq(length(unique(state))))

# dump out the stan data
stan_rdump(list = c('N', 'T', 'D', 'U', 
                    'sight', 'state', 'u', 'ufull', 'mass'),
           file = '../data/data_dump/trait_w_gaps.data.R')

library(plyr)
library(reshape2)
library(stringr)
library(rstan)
library(arm)
library(ape)
library(geiger)

source('../R/mung.r')

load('../data/body_mass.rdata')
load('../data/scaled_super.rdata')

dat <- read.csv('../data/mam-occs.csv', stringsAsFactors = FALSE)

occur <- clean.occurrence(dat)
occur <- occur[occur$name.bi %in% na.mass[, 1], ]

sp.oc <- split(occur, occur$name.bi)

first.occ <- laply(sp.oc, function(x) max(x$bins))
first.occ <- mapvalues(first.occ, unique(first.occ), rank(unique(first.occ)))
first.occ <- mapvalues(first.occ, sort(unique(first.occ)), 
                       rev(sort(unique(first.occ))))

sp.mass <- log(na.mass[order(na.mass[, 1]), 2])

info <- data.frame(sp = names(sp.oc), fad = first.occ, mass = sp.mass)
hot.fix <- str_replace(info$sp, pattern = '\\s', replacement = '_')
matchy <- name.check(spt, data.names = hot.fix)
spt.hot <- ape::drop.tip(spt, matchy$tree_not_data)
spt.hot <- vcv(spt.hot) / max(diag(vcv(spt.hot)))
ff <- match(hot.fix, rownames(spt.hot))
spt.hot <- spt.hot[ff, ff]

data <- list(N = nrow(info), T = length(unique(info$fad)), 
             trait = info$mass, year = info$fad, vcv = spt.hot)
with(data, {stan_rdump(list = c('N', 'T', 'trait', 'year', 'vcv'),
                       file = '../data/data_dump/coping_info.data.R')})

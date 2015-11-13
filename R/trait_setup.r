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

dat <- read.csv('../data/mam-occs.csv', stringsAsFactors = FALSE)

occur <- clean.occurrence(dat)
ss <- split(occur, occur$bins)
occur <- Reduce(rbind, llply(ss, function(x) x[duplicated(x$name.bi), ]))

# make data and tree match
name <- matrix(str_replace(occur$name.bi, ' ', '_'), ncol = 1)
hot.fix <- name.check(spt, data.names = name)
occur <- occur[!(occur$name.bi %in% 
                 str_replace(hot.fix$data_not_tree, '_', ' ')), ]
spt <- ape::drop.tip(spt, hot.fix$tree_not_data)

# for testing purposes
keep <- createDataPartition(occur$bins, p = 0.1)
keepname <- str_replace(unique(occur$name.bi[keep[[1]]]), ' ', '_')
hot.fix <- name.check(spt, data.names = keepname)
spt <- ape::drop.tip(spt, hot.fix$tree_not_data)
occur <- occur[keep[[1]], ]

# final step is name the variables
y <- as.numeric(as.factor(occur$comdiet))
K <- length(unique(y))
N <- length(y)
D <- 1
x <- matrix(1, ncol = 1, nrow = N)
cohort <- occur$bins / 2
C <- length(unique(cohort))
cohort <- mapvalues(cohort, from = unique(cohort), to = seq(C))
vcv <- vcv(spt)
vcv <- vcv / max(diag(vcv))
U <- nrow(vcv)
id <- as.numeric(as.factor(occur$name.bi))

stan_rdump(list = c('K', 'N', 'D', 'C', 'U', 'y', 'x', 'cohort', 'id', 'vcv'),
           file = '../data/data_dump/trait_info.data.R')

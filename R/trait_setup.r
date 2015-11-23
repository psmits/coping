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

occur <- clean.occurrence(dat)
ss <- split(occur, occur$bins)
occur <- Reduce(rbind, llply(ss, function(x) x[duplicated(x$name.bi), ]))

occur <- occur[occur$name.bi %in% na.mass$name, ]
occur$mass <- na.mass[match(occur$name.bi, na.mass$name), 2]

# make data and tree match
name <- matrix(str_replace(occur$name.bi, ' ', '_'), ncol = 1)
hot.fix <- name.check(spt, data.names = name)
occur <- occur[!(occur$name.bi %in% 
                 str_replace(hot.fix$data_not_tree, '_', ' ')), ]
spt <- ape::drop.tip(spt, hot.fix$tree_not_data)

# for testing purposes
keep <- createDataPartition(occur$bins, p = 0.5)
keepname <- str_replace(unique(occur$name.bi[keep[[1]]]), ' ', '_')
hot.fix <- name.check(spt, data.names = keepname)
spt <- ape::drop.tip(spt, hot.fix$tree_not_data)
occur <- occur[keep[[1]], ]

occur$mass <- scale(log(occur$mass))

# final step is name the variables
y <- as.numeric(interaction(occur$comdiet, occur$comlife))
K <- length(unique(y))
N <- length(y)
D <- 2
x <- matrix(1, ncol = 1, nrow = N)
x <- cbind(x, occur$mass)
cohort <- occur$bins / 2
C <- length(unique(cohort))
cohort <- mapvalues(cohort, from = unique(cohort), to = seq(C))
#vcv <- vcv(spt)
#vcv <- vcv / max(diag(vcv))
#U <- nrow(vcv)
#id <- as.numeric(as.factor(occur$name.bi))

stan_rdump(list = c('K', 'N', 'D', 'C', 'y', 'x', 'cohort'),
           file = '../data/data_dump/trait_info.data.R')

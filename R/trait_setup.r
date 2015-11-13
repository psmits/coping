library(plyr)
library(reshape2)
library(stringr)
library(rstan)
library(arm)
library(ape)
library(geiger)
library(caret)

source('../R/mung.r')

#load('../data/scaled_super.rdata')

dat <- read.csv('../data/mam-occs.csv', stringsAsFactors = FALSE)

occur <- clean.occurrence(dat)
ss <- split(occur, occur$bins)
occur <- Reduce(rbind, llply(ss, function(x) x[duplicated(x$name.bi), ]))
keep <- createDataPartition(occur$bins, p = 0.1)
occur <- occur[keep[[1]], ]

y <- as.numeric(as.factor(occur$comdiet))
K <- length(unique(y))
N <- length(y)
D <- 1
x <- matrix(1, ncol = 1, nrow = N)
cohort <- occur$bins / 2
C <- length(unique(cohort))

stan_rdump(list = c('K', 'N', 'D', 'C', 'y', 'x', 'cohort'),
           file = '../data/data_dump/trait_info.data.R')

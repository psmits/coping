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
keep <- createDataPartition(occur$bins, p = 0.4)
keepname <- str_replace(unique(occur$name.bi[keep[[1]]]), ' ', '_')
hot.fix <- name.check(spt, data.names = keepname)
spt <- ape::drop.tip(spt, hot.fix$tree_not_data)
occur <- occur[keep[[1]], ]


# process climate information
source('../R/mung_clim.r')


# split by genus to make summaries just in case
bygen <- split(occur, occur$occurrence.genus_name)
gen.trait <- data.frame(life = laply(bygen, function(x) 
                                     names(which.max(table(x$comlife)))),
                        diet = laply(bygen, function(x) 
                                     names(which.max(table(x$comdiet)))),
                        mass = laply(bygen, function(x) 
                                     mean(x$mass)))
for(ii in seq(length(bygen))) {
  bygen[[ii]][, c('comlife', 'comdiet', 'mass')] <- gen.trait[ii, ]
}
bygen <- llply(bygen, function(x) x[!(duplicated(x$bins)), ])
#occur <- Reduce(rbind, bygen)

# make covariates the right shape
occur$mass <- arm::rescale(log(occur$mass))
diet <- model.matrix( ~ occur$comdiet - 1)[, -1]


# final step is name the variables
y <- as.numeric(as.factor(occur$comlife))
K <- length(unique(y))
N <- length(y)
x <- matrix(1, ncol = 1, nrow = N)
x <- cbind(x, occur$mass, diet)
D <- ncol(x)
cohort <- occur$bins / 2
C <- length(unique(cohort))
cohort <- mapvalues(cohort, from = unique(cohort), to = seq(C))
#vcv <- vcv(spt)
#vcv <- vcv / max(diag(vcv))
#U <- nrow(vcv)
#id <- as.numeric(as.factor(occur$name.bi))
isoval <- mean.o18
isorang <- range.o18
### WARNING to include need to remove oldest bin!!!!
# tempval <- temp.time.mean
# isorang <- temp.time.range

stan_rdump(list = c('K', 'N', 'D', 'C', 'y', 'x', 
                    'cohort', 'isoval', 'isorang'),
           file = '../data/data_dump/trait_info.data.R')

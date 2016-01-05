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
occur <- occur[occur$name.bi %in% na.mass$name, ]
occur$mass <- na.mass[match(occur$name.bi, na.mass$name), 2]

# make data and tree match
name <- matrix(str_replace(occur$name.bi, ' ', '_'), ncol = 1)
hot.fix <- name.check(spt, data.names = name)
occur <- occur[!(occur$name.bi %in% 
                 str_replace(hot.fix$data_not_tree, '_', ' ')), ]
spt <- ape::drop.tip(spt, hot.fix$tree_not_data)

# for testing purposes
keep <- createDataPartition(occur$bins, p = 0.2)
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

occur$comlife <- as.factor(occur$comlife)
occur$comlife <- factor(occur$comlife, 
                        levels = c(levels(occur$comlife)[-4], 
                                   levels(occur$comlife)[4]))
# final step is name the variables
y <- as.numeric(occur$comlife)
K <- length(unique(y))
N <- length(y)

# covariates
x <- cbind(occur$mass, diet)  # for sep intercept set up
D <- ncol(x)

# temporal data
cohort <- occur$bins / 2
C <- length(unique(cohort))
cohort <- mapvalues(cohort, from = unique(cohort), to = seq(C))

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

# dump it out
stan_rdump(list = c('K', 'N', 'D', 'C', 'U', 'y', 'x', 'u', 'cohort'),
           file = '../data/data_dump/trait_info.data.R')

library(ape)
library(plyr)
library(phytools)
library(stringr)
library(paleotree)
library(phangorn)

load('../data/update_taxonomy.rdata')
source('../R/phylo_gen.r')
source('../R/na_mung.r')

raia.tree <- read.tree('../data/raia_tree.txt')
tom.tree <- read.nexus('../data/tomiya_tree.nex')

hal.list <- list.files('../data/halliday', full.names = TRUE)
# which of these is recommended by halliday?
halliday.tree <- read.nexus(hal.list[length(hal.list)])
# this is a genus tree...fuck
#   assign each genus a random species

fam.tree <- read.nexus('../data/meredith.nex')
super.tree <- read.nexus('../data/bininda_emonds.nex')


# make taxonomy tree of what i have
# north america
no.fam <- which(na.tax$family_name == '')
no.ord <- which(na.tax$order_name == '')
rms <- unique(c(no.fam, no.ord))
clean.na <- na.tax[-rms, ]
uni.tax <- unique(clean.na[, 1:4])

new.tax <- replace.taxonomy(dat, uni.tax[, 1:3])
# some things don't have known orders and families. 
# so i want them to be polytomy at the "root"
new.tax$order_name[new.tax$order_name == 'Artiodactyla'] <- 'Cetartiodactyla'

# split by orders
by.order <- split(new.tax, new.tax$order_name)
order.tree <- list()
for(ii in seq(length(by.order))) {
  check <- by.order[[ii]]
  check <- check[, c('family_name', 'occurrence.genus_name', 'name.bi')]
  check$family_name <- as.factor(check$family_name)
  check$occurrence.genus_name <- as.factor(check$occurrence.genus_name)
  check$name.bi <- str_replace(check$name.bi, ' ', '_')
  check$name.bi <- as.factor(check$name.bi)
  check <- unique(check)

  # screen situations with only 1 family, causes problems
  if(length(levels(check[, 1])) == 1) {
    temp <- as.phylo.formula(~ occurrence.genus_name/
                             name.bi, 
                             data = check)
  } else {
    temp <- as.phylo.formula(~ family_name/
                             occurrence.genus_name/
                             name.bi, 
                             data = check)
  }
  temp <- collapse.singles(temp)

  # make the unknown family node "disappear"
  temp <- unitLengthTree(temp)
  roo <- getMRCA(temp, check$name.bi[check$family_name == ''])
  temp$edge.length[apply(temp$edge, 1, function(x) any(x == roo))] <- 0
  temp <- di2multi(temp)

  order.tree[[ii]] <- temp
}
names(order.tree) <- names(by.order)
ntip <- laply(order.tree, function(x) length(x$tip.label))
order.tree <- order.tree[ntip != 1]
# combine the orders into a huge tree
taxon.tree <- Reduce(bind.tree, order.tree)

species.trees <- list(raia.tree, super.tree[[1]], taxon.tree)
class(species.trees) <- 'multiPhylo'

big.tree <- mrp.supertree(species.trees)


# get rid of the stupid tips
if(class(big.tree) == 'multiPhylo') {
  spt <- big.tree[[1]] 
} else {
  spt <- big.tree
}

fad <- ddply(dat, .(name.bi), summarize, max(bins))
lad <- ddply(dat, .(name.bi), summarize, min(bins))
fad[, 1] <- str_replace(fad[, 1], ' ', '_')
lad[, 1] <- str_replace(lad[, 1], ' ', '_')
datmat <- cbind(fad[, 2], lad[, 2])
rownames(datmat) <- fad[, 1]

check <- name.check(spt, datmat)
spt <- drop.tip(spt, check$tree_not_data)
datmat <- datmat[!(rownames(datmat) %in% check$data_not_tree), ]
datmat <- datmat[match(spt$tip.label, rownames(datmat)), ]

spt <- timeLadderTree(spt, timeData = datmat)
spt <- timePaleoPhy(spt, timeData = datmat, 
                         type = 'mbl', vartime = 0.1)
# how do i really want to scale this?
save(spt, file = '../data/scaled_super.rdata')

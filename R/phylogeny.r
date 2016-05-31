library(ape)
library(plyr)
library(phytools)
library(stringr)
library(paleotree)
library(phangorn)
library(geiger)

load('../data/update_taxonomy.rdata')
source('../R/phylo_gen.r')
source('../R/taxon_names.r')
source('../R/mung.r')

dat <- read.csv('https://paleobiodb.org/data1.2/occs/list.csv?datainfo&rowcount&base_name=Mammalia&taxon_reso=species&interval=Maastrichtian,Gelasian&cc=NOA&show=class,genus,ecospace,strat,stratext,lith,acconly', stringsAsFactors = FALSE, skip = 20)
occur <- clean.occurrence(dat)

raia.tree <- read.tree('../data/raia_tree.txt')
tom.tree <- read.nexus('../data/tomiya_tree.nex')

hal.list <- list.files('../data/halliday', full.names = TRUE)
# which of these is recommended by halliday?
halliday.tree <- read.nexus(hal.list[length(hal.list)])
dat$
# this is a genus tree...fuck
#   assign each genus a random species

fam.tree <- read.nexus('../data/meredith.nex')
super.tree <- read.nexus('../data/bininda_emonds.nex')


make.genus <- function(phy) {
  new.tip <- str_extract(phy$tip.label, pattern = '[^\\_]*')
  dups <- duplicated(new.tip)
  phy <- drop.tip(phy, phy$tip.label[dups])
  phy$tip.label <- str_extract(phy$tip.label, pattern = '[^\\_]*')
  phy
}

raia.tree <- make.genus(raia.tree)
tom.tree <- make.genus(tom.tree)
super.tree <- llply(super.tree, make.genus)


# make taxonomy tree of what i have
# north america
no.fam <- which(na.tax$family_name == '')
no.ord <- which(na.tax$order_name == '')
rms <- unique(c(no.fam, no.ord))
clean.na <- na.tax[-rms, ]
uni.tax <- unique(clean.na[, 1:4])

new.tax <- replace.taxonomy(occur, uni.tax[, 1:3])

# some things don't have known orders and families. 
# so i want them to be polytomy at the "root"
new.tax$order_name[new.tax$order_name == 'Artiodactyla'] <- 'Cetartiodactyla'
new.tax <- new.tax[!(duplicated(new.tax$name.bi)), ]

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
  if(length(temp$tip.label) > 1) {
    temp <- collapse.singles(temp)
  }

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
taxon.tree <- make.genus(taxon.tree)

tom.tree$tip.label


genus.trees <- list(raia.tree,
                    tom.tree,
                    halliday.tree[[1]],
                    super.tree[[1]], 
                    taxon.tree)
class(genus.trees) <- 'multiPhylo'

big.tree <- mrp.supertree(genus.trees)

# get rid of the stupid tips
if(class(big.tree) == 'multiPhylo') {
  spt <- big.tree[[1]] 
} else {
  spt <- big.tree
}

fad <- ddply(occur, .(occurrence.genus_name), summarize, max(bins))
lad <- ddply(occur, .(occurrence.genus_name), summarize, min(bins))
fad[, 1] <- str_replace(fad[, 1], ' ', '_')
lad[, 1] <- str_replace(lad[, 1], ' ', '_')
datmat <- cbind(fad[, 2], lad[, 2])
rownames(datmat) <- fad[, 1]

check <- name.check(spt, datmat)

scale.tree <- function(phy, dat) {
  check <- name.check(phy, dat)
  spt <- drop.tip(phy, check$tree_not_data)
  datmat <- dat[!(rownames(datmat) %in% check$data_not_tree), ]
  datmat <- datmat[match(spt$tip.label, rownames(datmat)), ]
  
  #spt <- timeLadderTree(spt, timeData = datmat)
  spt <- timePaleoPhy(spt, timeData = datmat, 
                           type = 'mbl', vartime = 0.1, 
                           randres = TRUE,
                           ntrees = 100)
}
spt <- scale.tree(spt, datmat)

# how do i really want to scale this?
save(spt, file = '/home/psmits/coping/data/scaled_super.rdata')

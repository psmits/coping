library(xtable)

source('../R/mung.r')

dat <- read.csv('../data/pbdb_data.csv', stringsAsFactors = FALSE, skip = 21)
occur <- clean.occurrence(dat)
posture <- read.csv('../data/posture.csv', stringsAsFactors = FALSE) 

load('../data/update_taxonomy.rdata')
na.tax <- na.tax[!(duplicated(na.tax$name.bi)), ]
fx <- na.tax$name.bi %in% occur$name.bi[occur$comlife == 'ground dwelling']
na.fx <- na.tax[fx, ]
occur[match(na.fx$name.bi, occur$name.bi), 
      c('order', 'family')] <- na.fx[, 1:2]

mats <- posture[, 1] %in% unique(na.tax$order_name) | 
        posture[, 1] %in% unique(na.tax$family_name)

posture <- posture[mats, ]

posture$family <- posture$taxon
posture$order <- posture$taxon

posture$family[!(posture[, 1] %in% unique(na.tax$family_name))] <- NA
posture$order[!(posture[, 1] %in% unique(na.tax$order_name))] <- NA

posture <- posture[order(posture$taxon), ]
posture <- posture[, c(4, 3, 2)]

post.tab <- xtable(posture)
print.xtable(x = post.tab, 
             file = '../doc/posture_raw.tex', 
             include.rownames = FALSE)

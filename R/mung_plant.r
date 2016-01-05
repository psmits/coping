library(plyr)
library(neotoma)

gpids <- get_table(table.name='GeoPoliticalUnits')
usa <- gpids[gpids$GeoPoliticalName == 'United States', ]
can <- gpids[gpids$GeoPoliticalName == 'Canada', ]
mex <- gpids[gpids$GeoPoliticalName == 'Mexico', ]
pol.dst.usa <- get_dataset(datasettype = 'pollen', 
                           gpid = usa[, 1], ageold = 66 * 10^6)
pol.dst.can <- get_dataset(datasettype = 'pollen', 
                           gpid = can[, 1], ageold = 66 * 10^6)
pol.dst.mex <- get_dataset(datasettype = 'pollen', 
                           gpid = mex[, 1], ageold = 66 * 10^6)
pol.usa <- get_download(pol.dst.usa)

#plant <- read.csv('../data/plant_occ.csv', stringsAsFactors = FALSE)
#b <- cbind(seq(from = 0, to = 64, by = 2), seq(from = 2, to = 66, by = 2))
#mem <- array(dim = nrow(plant))
#
#for(ii in seq(nrow(b))) {
#  mem[plant$ma_mid >= b[ii, 1] & plant$ma_mid < b[ii, 2]] <- ii
#}
#plant <- plant[!(is.na(mem)), ]
#mem <- mem[!(is.na(mem))]

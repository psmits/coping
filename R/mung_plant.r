library(plyr)

plant <- read.csv('../data/plant_occ.csv', stringsAsFactors = FALSE)
b <- cbind(seq(from = 0, to = 64, by = 2), seq(from = 2, to = 66, by = 2))
mem <- array(dim = nrow(plant))

for(ii in seq(nrow(b))) {
  mem[plant$ma_mid >= b[ii, 1] & plant$ma_mid < b[ii, 2]] <- ii
}
plant <- plant[!(is.na(mem)), ]
mem <- mem[!(is.na(mem))]

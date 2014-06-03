traits <- read.csv("Traits2006J.csv", row.names = 1)
community <- read.csv("Comm2006J.csv", row.names = 1)
require(betapart)
beta.results2006J <- functional.beta.pair(x = community, traits = as.matrix(traits))
save(beta.results2006J, file = "betaresults2006J.RData")

?read.csv

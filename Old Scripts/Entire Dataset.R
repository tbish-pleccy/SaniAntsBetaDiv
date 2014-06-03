head(presabs)

library(betapart)
library(sp)

head(site.info)
site.info$s.id <- as.numeric(site.info$s.id)
rownames(site.info) <- site.info$s.id
site.info <- site.info[order(site.info$s.id), ]

presabs <- comm
presabs <- subset(presabs, select = - c(s.id, Site, rep, year, seasn))
rownames(presabs) <- comm$s.id
presabs[presabs > 0] <- 1
presabs <- presabs[, order(colnames(presabs))] # alphabetically order!
# Only use communities with 5 or more species
#presabs <- presabs[as.vector(which(rowSums(presabs) > 4)), ]
site.info$ric <- rowSums(presabs)

comms2use <- site.info[which(site.info$ric > 4), "s.id"]

## ele dist

elev <- site.info[site.info$s.id %in% comms2use, c("altitude", "s.id")]
elev.dist <- dist(elev$altitude)

## spa dist

spat <- site.info[site.info$s.id %in% comms2use, c("x", "y")]
spat.dist <- as.dist(spDists(as.matrix(spat), longlat = TRUE))

## yearly dist

years <- site.info[site.info$s.id %in% comms2use, c("year", "s.id")]
year.dist <- dist(years$year)

## seasonal dist

seas <- site.info[site.info$s.id %in% comms2use, c("seasn", "s.id")]
seas$seasn <- ifelse(seas$seasn == "J", 0, 1)
seas.dist <- dist(seas$seasn)

## comm dist
comm.restricted <- presabs[rownames(presabs) %in% comms2use, ]

comm.all <- beta.pair(comm.restricted)
str(comm.all)
comm.sor <- comm.all$beta.sor
comm.sim <- comm.all$beta.sim
comm.nes <- comm.all$beta.nes

## func dist
func.all <- functional.beta.pair(x = comm.restricted, traits = pcoa.axes)
func.sor <- func.all$beta.sor
func.sim <- func.all$beta.sim
func.nes <- func.all$beta.sne



traits.temp <- pcoa.axes[rownames(pcoa.axes) %in% present.species, ]
fun.dist <- functional.beta.pair(x = c.pa, traits = traits.temp)






test1 <- MRM(comm.sim ~ elev.dist + spat.dist + year.dist + seas.dist, nperm = 1000)
test1

library(hier.part)
hier.part(y = as.vector(comm.sim), xcan = data.frame(as.vector(elev.dist), 
                                              as.vector(spat.dist), 
                                              as.vector(year.dist), 
                                              as.vector(seas.dist))) 

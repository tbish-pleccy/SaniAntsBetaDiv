# Species only present in either season

head(comm)

spec.frame <- NULL

for (i in unique(comm$Site)){
  for (j in unique(comm$year)){
    communityJ <- comm[comm$seasn == "J" & comm$Site == i & comm$year == j, ]
    communityJ <- subset(communityJ, select = - c(s.id, Site, rep, year, seasn))
    rownames(communityJ) <- comm[comm$seasn == "J" & comm$Site == i & comm$year == j, "s.id"]
    communityJ <- communityJ[, order(colnames(communityJ))] # alphabetically order!
    communityJ <- communityJ[as.vector(which(rowSums(communityJ) > 0)), ] # remove empty communities
    absent.speciesJ <- which(colSums(communityJ) == 0)
    communityJ <- communityJ[, - absent.speciesJ]  # remove empty species
    #head(communityJ)
    #length(communityJ)
    
    length(communityJ[, which(colSums(communityJ) > 10)])
    
    communityS <- comm[comm$seasn == "S" & comm$Site == i & comm$year == j, ]
    communityS <- subset(communityS, select = - c(s.id, Site, rep, year, seasn))
    rownames(communityS) <- comm[comm$seasn == "S" & comm$Site == i & comm$year == j, "s.id"]
    communityS <- communityS[, order(colnames(communityS))] # alphabetically order!
    communityS <- communityS[as.vector(which(rowSums(communityS) > 0)), ] # remove empty communities
    absent.speciesS <- which(colSums(communityS) == 0)
    communityS <- communityS[, - absent.speciesS]  # remove empty species
    #head(communityS)
    #length(communityS)
    
    
    janpres <- (data.frame(cbind(names(communityJ), 1)))
    names(janpres) <- c("species", "J")
    janpres$J <- as.numeric(janpres$J)
    seppres <- (data.frame(cbind(names(communityS), 1)))
    names(seppres) <- c("species", "S")
    seppres$S <- as.numeric(seppres$S)
    
    seapres <- merge(janpres, seppres, by = "species", all.x = TRUE, all.y = TRUE)
    seapres[is.na(seapres)] <- 0
    
    S.specialists <- nrow(seapres[seapres$J == 0, ])
    J.specialists <- nrow(seapres[seapres$S == 0, ])
    comm.tot <- comm[comm$Site == i & comm$year == j, ]
    comm.tot <- subset(comm.tot, select = - c(s.id, Site, rep, year, seasn))
    absent.species <- which(colSums(comm.tot) == 0)
    comm.tot <- comm.tot[, - absent.species]  # remove empty species
    
    pool.size <- length(comm.tot)
    Site <- i
    year <- j
    temp <- data.frame(Site, year, pool.size, J.specialists, S.specialists)
    spec.frame <- rbind(spec.frame, temp)
    
  }
  
  
}

spec.frame
spec.frame$prop.J.spec <- (spec.frame$J.specialists / spec.frame$pool.size) * 100
tapply(spec.frame$prop.J.spec, list(spec.frame$Site), mean)
par(mfrow = c(1, 1))
plot(tapply(spec.frame$prop.J.spec, list(spec.frame$Site), mean))
head(site.infoS)
spec.frame <- merge(spec.frame, site.infoS)
plot(tapply(spec.frame$prop.J.spec, list(spec.frame$altitude), mean))

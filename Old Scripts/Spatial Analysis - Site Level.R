# Spatial Analysis with Partial Mantels
head(comm)
commS <- subset(comm, select = - c(s.id, rep))
names(commS)
#aggregate(commS, by = list(commS$Site), FUN = sum)
library(reshape)
#commS <- melt(commS, id = c("Site", "year", "seasn"))
#commS <- cast(commS, Site + year + seasn ~ variable, value = "value", FUN = sum)
commS <- aggregate(. ~ Site + year + seasn, commS, sum)
head(commS)
nrow(commS)

head(site.infoS)
site.infoS <- subset(site.info, select = - c(s.id, rep))
site.infoS$altitude <- as.numeric(site.infoS$altitude)
site.infoS$pixels <- as.numeric(site.infoS$pixels)

site.infoS <- aggregate(. ~ Site + year + seasn, site.infoS, mean)
str(site.infoS)

commS$s.id <- paste(commS$year, commS$seasn, commS$Site, sep = "")
rownames(commS) <- commS$s.id
site.infoS$s.id <- paste(site.infoS$year, site.infoS$seasn, site.infoS$Site, sep = "")
communityS<- subset(commS, select = - c(s.id, Site, year, seasn))

# Get packages
library(FD)
library(betapart)
library(vegan)
library(boot)
library(ecodist)
library(sp)

# Create empty objects

taxonomic.partial.mantels.compS <- NULL
functional.partial.mantels.compS <- NULL
comm.dissimS <- list()
fun.dissimS <- list()
ele.dissimS <- list()
spa.dissimS <- list()

i <- 2009
j <- "J"
# Start loop
starttime <- Sys.time()

for (i in unique(commS$year)){
  for (j in unique(commS$seasn)){
    c.ab <- commS[commS$year == i & commS$seasn == j, ]  # c.ab is abundance data
    
    c.pa <- subset(c.ab, select = - c(s.id, Site, year, seasn))
    rownames(c.pa) <- c.ab$s.id
    c.pa[c.pa > 0] <- 1  # from here, c.pa is presence/absence data
    c.pa <- c.pa[ which(rowSums(c.pa) > 5), ]  # use only comms with >5 species
    absent.species <- which(colSums(c.pa) == 0)
    c.pa <- c.pa[, - absent.species]  # remove empty species
    c.pa <- c.pa[as.vector(which(rowSums(c.pa) > 0)), ]  # remove empty communities
    c.pa <- c.pa[, order(colnames(c.pa))]  # alphabetically order!
    present.species <- as.vector(names(c.pa))  # remeber present species
    
    
    # Create appropriate site info frame.
    site.info.temp <- site.infoS[site.infoS$year == i & site.infoS$seasn == j, ]
    rownames(site.info.temp) <- site.info.temp$s.id
    #site.info.temp$s.id <- as.numeric(site.info.temp$s.id)
    site.info.temp <- site.info.temp[order(site.info.temp$s.id), ]
    
    #rownames(c.pa) <- as.numeric(rownames(c.pa))
    site.info.temp <- site.info.temp[which(site.info.temp$s.id %in% 
                                             (rownames(c.pa))), ]
    
    # Create altitudinal distance matrix
    altitudes <- data.frame(site.info.temp[, "altitude"])
    rownames(altitudes) <- rownames(site.info.temp)
    ele.dist <- dist(altitudes)
    
    ele.dissimS[[paste(i, j, sep = "")]] <- ele.dist
    
    # create spatial distance matrix
    coord <- as.matrix(subset(site.info.temp, select = c(x, y)))
    spa.dist <- as.dist(spDists(coord, longlat = TRUE))
    
    spa.dissimS[[paste(i, j, sep = "")]] <- spa.dist
    
    # Create community distance matrices
    comm.dist <- beta.pair(c.pa)
    comm.dist.bsor <- comm.dist$beta.sor
    comm.dist.bsim <- comm.dist$beta.sim
    comm.dist.bnes <- comm.dist$beta.sne
    
    comm.dissimS[[paste(i, j, "bsor", sep = "")]] <- comm.dist.bsor
    comm.dissimS[[paste(i, j, "bsim", sep = "")]] <- comm.dist.bsim
    comm.dissimS[[paste(i, j, "bnes", sep = "")]] <- comm.dist.bnes    
    
    # Functional beta diversity
    traits.temp <- pcoa.axes[rownames(pcoa.axes) %in% present.species, ]
    fun.dist <- functional.beta.pair(x = c.pa, traits = traits.temp)
    fun.dist.bsor <- fun.dist$funct.beta.sor
    fun.dist.bsim <- fun.dist$funct.beta.sim
    fun.dist.bsne <- fun.dist$funct.beta.sne
    
    fun.dissimS[[paste(i, j, "bsor", sep = "")]] <- fun.dist.bsor
    fun.dissimS[[paste(i, j, "bsim", sep = "")]] <- fun.dist.bsim
    fun.dissimS[[paste(i, j, "bsne", sep = "")]] <- fun.dist.bsne
    
    #taxonomic partial mantels
    tax.sor.ele <- mantel(comm.dist.bsor ~ ele.dist + spa.dist, 
                          nperm = 1000, nboot = 1000)
    tax.sor.spa <- mantel(comm.dist.bsor ~ spa.dist + ele.dist, 
                          nperm = 1000, nboot = 1000)
    
    tax.sim.ele <- mantel(comm.dist.bsim ~ ele.dist + spa.dist, 
                          nperm = 1000, nboot = 1000)
    tax.sim.spa <- mantel(comm.dist.bsim ~ spa.dist + ele.dist, 
                          nperm = 1000, nboot = 1000)
    
    tax.nes.ele <- mantel(comm.dist.bnes ~ ele.dist + spa.dist, 
                          nperm = 1000, nboot = 1000)
    tax.nes.spa <- mantel(comm.dist.bnes ~ spa.dist + ele.dist, 
                          nperm = 1000, nboot = 1000)
    
    taxonomic.partial.mantels <- data.frame(rbind(tax.sor.ele, tax.sor.spa, 
                                                  tax.sim.ele, tax.sim.spa, 
                                                  tax.nes.ele, tax.nes.spa))
    
    taxonomic.partial.mantels$year <- i
    taxonomic.partial.mantels$seasn <- j
    taxonomic.partial.mantels$test <- rownames(taxonomic.partial.mantels)
    taxonomic.partial.mantels$div <- 
      substring(taxonomic.partial.mantels$test, 1, 3)
    taxonomic.partial.mantels$index <- 
      substring(taxonomic.partial.mantels$test, 5, 7)
    taxonomic.partial.mantels$gradient <- 
      substring(taxonomic.partial.mantels$test, 9, 11)
    
    taxonomic.partial.mantels.compS <- rbind(taxonomic.partial.mantels.compS, 
                                            taxonomic.partial.mantels)
    rownames(taxonomic.partial.mantels.compS) <- 
      seq(1, nrow(taxonomic.partial.mantels.compS), 1)
    
    #functional partial mantels
    fun.sor.ele <- mantel(fun.dist.bsor ~ ele.dist + spa.dist, 
                          nperm = 1000, nboot = 1000)
    fun.sor.spa <- mantel(fun.dist.bsor ~ spa.dist + ele.dist, 
                          nperm = 1000, nboot = 1000)
    
    fun.sim.ele <- mantel(fun.dist.bsim ~ ele.dist + spa.dist, 
                          nperm = 1000, nboot = 1000)
    fun.sim.spa <- mantel(fun.dist.bsim ~ spa.dist + ele.dist, 
                          nperm = 1000, nboot = 1000)
    
    fun.sne.ele <- mantel(fun.dist.bsne ~ ele.dist + spa.dist, 
                          nperm = 1000, nboot = 1000)
    fun.sne.spa <- mantel(fun.dist.bsne ~ spa.dist + ele.dist, 
                          nperm = 1000, nboot = 1000)
    
    functional.partial.mantels <- data.frame(rbind(fun.sor.ele, fun.sor.spa, 
                                                   fun.sim.ele, fun.sim.spa, 
                                                   fun.sne.ele, fun.sne.spa))
    functional.partial.mantels$year <- i
    functional.partial.mantels$seasn <- j
    functional.partial.mantels$test <- rownames(functional.partial.mantels)
    functional.partial.mantels$div <- 
      substring(functional.partial.mantels$test, 1, 3)
    functional.partial.mantels$index <- 
      substring(functional.partial.mantels$test, 5, 7)
    functional.partial.mantels$gradient <- 
      substring(functional.partial.mantels$test, 9, 11)
    
    functional.partial.mantels.compS <- rbind(functional.partial.mantels.compS, 
                                             functional.partial.mantels)
    rownames(functional.partial.mantels.compS) <- 
      seq(1, nrow(functional.partial.mantels.compS), 1)
  }
  
  
}


taxonomic.partial.mantels.compS
functional.partial.mantels.compS

endtime <- Sys.time()
endtime - starttime  # 45 mins


taxonomic.partial.mantels.compS[taxonomic.partial.mantels.compS$seasn == "J" & 
                                 taxonomic.partial.mantels.compS$gradient == "ele" &
                                 taxonomic.partial.mantels.compS$index == "sor", ]
functional.partial.mantels.compS[functional.partial.mantels.compS$seasn == "J" & 
                                  functional.partial.mantels.compS$gradient == "ele" & 
                                  functional.partial.mantels.compS$index == "sor", ]


taxonomic.partial.mantels.compS[taxonomic.partial.mantels.compS$seasn == "J" & 
                                 taxonomic.partial.mantels.compS$gradient == "ele" &
                                 taxonomic.partial.mantels.compS$index == "sim", ]
functional.partial.mantels.compS[functional.partial.mantels.compS$seasn == "J" & 
                                  functional.partial.mantels.compS$gradient == "ele" & 
                                  functional.partial.mantels.compS$index == "sim", ]


taxonomic.partial.mantels.compS[taxonomic.partial.mantels.compS$seasn == "J" & 
                                 taxonomic.partial.mantels.compS$gradient == "ele" &
                                 taxonomic.partial.mantels.compS$index == "nes", ]
functional.partial.mantels.compS[functional.partial.mantels.compS$seasn == "J" & 
                                  functional.partial.mantels.compS$gradient == "ele" & 
                                  functional.partial.mantels.compS$index == "sne", ]



taxonomic.partial.mantels.compS[taxonomic.partial.mantels.compS$seasn == "S" & 
                                 taxonomic.partial.mantels.compS$gradient == "ele" &
                                 taxonomic.partial.mantels.compS$index == "sor", ]
functional.partial.mantels.compS[functional.partial.mantels.compS$seasn == "S" & 
                                  functional.partial.mantels.compS$gradient == "ele" & 
                                  functional.partial.mantels.compS$index == "sor", ]

taxonomic.partial.mantels.compS[taxonomic.partial.mantels.compS$seasn == "S" & 
                                 taxonomic.partial.mantels.compS$gradient == "ele" &
                                 taxonomic.partial.mantels.compS$index == "sim", ]
functional.partial.mantels.compS[functional.partial.mantels.compS$seasn == "S" & 
                                  functional.partial.mantels.compS$gradient == "ele" & 
                                  functional.partial.mantels.compS$index == "sim", ]

taxonomic.partial.mantels.compS[taxonomic.partial.mantels.compS$seasn == "S" & 
                                 taxonomic.partial.mantels.compS$gradient == "ele" &
                                 taxonomic.partial.mantels.compS$index == "nes", ]
functional.partial.mantels.compS[functional.partial.mantels.compS$seasn == "S" & 
                                  functional.partial.mantels.compS$gradient == "ele" & 
                                  functional.partial.mantels.compS$index == "sne", ]




par(mfrow = c(3, 4))
metrics <- c("bsor", "bsim", "bnes")
for (k in unique(metrics)){
  for(m in unique(c("tax", "fun"))){
    if (m == "tax"){
      dissim <- comm.dissimS
      partial.mantels.comp <- taxonomic.partial.mantels.compS
    }
    else {
      dissim <- fun.dissimS
      partial.mantels.comp <- functional.partial.mantels.compS
    }
    
    
    if(k == "bnes" && m == "fun"){
      k <- "bsne"
    }
    else {
      k <- k
    }
    
    
    for (j in unique(commS$seasn)){
      ylabel <- paste("Community Dissimilarity Residuals:", k, sep = " ")
      plot(x = c(-500, 500), y = c(-0.25, 0.25), ylab = ylabel, xlab = "Elevational Distance Residuals", type = "n")
      
      for (i in unique(commS$year)){
        c.dist <- dissim[[paste(i, j, k, sep = "")]]
        e.dist <- ele.dissimS[[paste(i, j, sep = "")]]
        s.dist <- spa.dissimS[[paste(i, j, sep = "")]]
        
        residual_vectors <- cbind(c.dist, e.dist, s.dist)
        resid_vecs <- as.data.frame(residual_vectors)
        
        comm_residualsG <- resid(lm(c.dist ~ e.dist, data = resid_vecs))
        eco_residuals <- resid(lm(s.dist ~ e.dist, data = resid_vecs))
        comm_residualsE <- resid(lm(c.dist ~ s.dist, data = resid_vecs))
        geo_residuals <- resid(lm(e.dist ~ s.dist, data = resid_vecs))
        
        ind <- substring(k, 2, 4)
        
        
        
        
        
        elep <- partial.mantels.comp[partial.mantels.comp$index == ind & 
                                       partial.mantels.comp$gradient == "ele" & 
                                       partial.mantels.comp$year == i & 
                                       partial.mantels.comp$seasn == j, "pval3"]
        
        
        
        
        linecol <- ifelse(elep < 0.05, "red", "blue")
        eline <- ifelse(elep < 0.05, 1, 3)
        abline(lm(comm_residualsE~-1+geo_residuals), col = linecol, lty = eline)
        
        
      }
    }
  }
}


comm.dissimS
communityS[c("2012J7", "2012J1"), - (which(colSums(communityS[c("2012J7", "2012J1"), ]) == 0))]
communityS[c("2012S2", "2012S1"), - (which(colSums(communityS[c("2012S2", "2012S1"), ]) == 0))]
communityS[c("2010S2", "2010S1"), - (which(colSums(communityS[c("2010S2", "2010S1"), ]) == 0))]




par(mfrow = c(2, 7))
for (i in unique(commS$year)){
  for (j in unique(commS$seasn)){
    dat <- commS[commS$year == i & commS$seasn == j, ]
    dat <- subset(dat, select = - c(s.id, Site, year, seasn))
    present <- as.vector(names(dat[, - (which(colSums(dat) == 0))]))
    y2 <- y[which(rownames(y) %in% present), ]
    ran.temp <- rangesize(y2)
    hist(ran.temp, main = paste(i, j), xlim = c(0, 2200), breaks = 8)
  }

}

ran <- rangesize(y)

mids
hist(ran)


site.info.temp <- site.info.temp[which(site.info.temp$s.id %in% 
                                         (rownames(c.pa))), ]
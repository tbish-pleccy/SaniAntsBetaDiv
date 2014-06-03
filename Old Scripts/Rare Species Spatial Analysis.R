head(commS)

hist(colSums(communityS), breaks = seq(0, 7500, 100))
# So lets call the rare species those with an abundance less than 100 over
# the entire dataset?

rare.species <- as.vector(names(which(colSums(communityS) < 1000)))
str(rare.species)

hist(colSums(communityS[, - which(colnames(communityS)  %in% rare.species)]), breaks = seq(0, 7500, 100))


#########################################################

# Get packages
library(FD)
library(betapart)
library(vegan)
library(boot)
library(ecodist)
library(sp)

# Create empty objects

taxonomic.partial.mantels.compR <- NULL
#functional.partial.mantels.compS <- NULL
rare.partial.mantels.compR <- NULL
comm.dissimR <- list()
rare.dissimR <- list()
#fun.dissimR <- list()
ele.dissimR <- list()
spa.dissimR <- list()
ele.dissim.rR <- list()
spa.dissim.rR <- list()

starttime <- Sys.time()

for (i in unique(comm$year)){
  for (j in unique(comm$seasn)){
    c.ab <- comm[comm$year == i & comm$seasn == j, ]  # c.ab is abundance data
    
    c.pa <- subset(c.ab, select = - c(s.id, Site, year, seasn, rep))
    rownames(c.pa) <- c.ab$s.id
    c.pa[c.pa > 0] <- 1  # from here, c.pa is presence/absence data
    c.pa <- c.pa[ which(rowSums(c.pa) > 5), ]  # use only comms with >5 species
    absent.species <- which(colSums(c.pa) == 0)
    c.pa <- c.pa[, - absent.species]  # remove empty species
    c.pa <- c.pa[as.vector(which(rowSums(c.pa) > 0)), ]  # remove empty communities
    c.pa <- c.pa[, order(colnames(c.pa))]  # alphabetically order!
    present.species <- as.vector(names(c.pa))  # remeber present species
    
    
    # Create appropriate site info frame.
    site.info.temp <- site.info[site.info$year == i & site.info$seasn == j, ]
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
    
    ele.dissimR[[paste(i, j, sep = "")]] <- ele.dist
    
    # create spatial distance matrix
    coord <- as.matrix(subset(site.info.temp, select = c(x, y)))
    spa.dist <- as.dist(spDists(coord, longlat = TRUE))
    
    spa.dissimR[[paste(i, j, sep = "")]] <- spa.dist
    
    # Create community distance matrices
    comm.dist <- beta.pair(c.pa)
    comm.dist.bsor <- comm.dist$beta.sor
    comm.dist.bsim <- comm.dist$beta.sim
    comm.dist.bnes <- comm.dist$beta.sne
    
    comm.dissimR[[paste(i, j, "bsor", sep = "")]] <- comm.dist.bsor
    comm.dissimR[[paste(i, j, "bsim", sep = "")]] <- comm.dist.bsim
    comm.dissimR[[paste(i, j, "bnes", sep = "")]] <- comm.dist.bnes    
    
    #     # Functional beta diversity
    #     traits.temp <- pcoa.axes[rownames(pcoa.axes) %in% present.species, ]
    #     fun.dist <- functional.beta.pair(x = c.pa, traits = traits.temp)
    #     fun.dist.bsor <- fun.dist$funct.beta.sor
    #     fun.dist.bsim <- fun.dist$funct.beta.sim
    #     fun.dist.bsne <- fun.dist$funct.beta.sne
    #     
    #     fun.dissimR[[paste(i, j, "bsor", sep = "")]] <- fun.dist.bsor
    #     fun.dissimR[[paste(i, j, "bsim", sep = "")]] <- fun.dist.bsim
    #     fun.dissimR[[paste(i, j, "bsne", sep = "")]] <- fun.dist.bsne
    
    
    # create common species only distance matrix
    c.ab2 <- subset(c.ab, select = - c(s.id, Site, year, seasn, rep))
    c.ab2 <- c.ab2[, - which(colSums(c.ab2) == 0)]  # remove empty species
    
    rare.species2 <- names(which(colSums(c.ab2)/sum(colSums(c.ab2)) < 0.04))
    
    c.pa.r <- c.pa[, which(colnames(c.pa) %in% rare.species2)]
    c.pa.r <- c.pa.r[as.vector(which(rowSums(c.pa.r) > 0)), ]  # remove empty communities
    # Create altitudinal distance matrix for rare
    altitudes.r <- data.frame(site.info.temp[which(rownames(site.info.temp)
                                                   %in% rownames(c.pa.r)), "altitude"])
    ele.dist.r <- dist(altitudes.r)
    ele.dissim.rR[[paste(i, j, sep = "")]] <- ele.dist.r
    # create spatial distance matrix for rare
    coord <- as.matrix(subset(site.info.temp, select = c(x, y)))
    coord.r <- coord[which(rownames(coord) %in% rownames(c.pa.r)), ]
    spa.dist.r <- as.dist(spDists(coord.r, longlat = TRUE))
    
    spa.dissim.rR[[paste(i, j, sep = "")]] <- spa.dist.r
    
    
    rare.dist <- beta.pair(c.pa.r)
    rare.dist.bsor <- rare.dist$beta.sor
    rare.dist.bsim <- rare.dist$beta.sim
    rare.dist.bnes <- rare.dist$beta.sne
    
    rare.dissimR[[paste(i, j, "bsor", sep = "")]] <- rare.dist.bsor
    rare.dissimR[[paste(i, j, "bsim", sep = "")]] <- rare.dist.bsim
    rare.dissimR[[paste(i, j, "bnes", sep = "")]] <- rare.dist.bnes
    
    # rare only partial mantels
    rar.sor.ele <- mantel(comm.dist.bsor ~ ele.dist.r + rare.dist.bsor, 
                          nperm = 1000, nboot = 1000)
    rar.sor.spa <- mantel(comm.dist.bsor ~ spa.dist.r + rare.dist.bsor, 
                          nperm = 1000, nboot = 1000)
    
    rar.sim.ele <- mantel(comm.dist.bsim ~ ele.dist.r + rare.dist.bsim, 
                          nperm = 1000, nboot = 1000)
    rar.sim.spa <- mantel(comm.dist.bsim ~ spa.dist.r + rare.dist.bsim, 
                          nperm = 1000, nboot = 1000)
    
    rar.nes.ele <- mantel(comm.dist.bnes ~ ele.dist.r + rare.dist.bnes, 
                          nperm = 1000, nboot = 1000)
    rar.nes.spa <- mantel(comm.dist.bnes ~ spa.dist.r + rare.dist.bnes, 
                          nperm = 1000, nboot = 1000)
    
    rare.partial.mantels <- data.frame(rbind(rar.sor.ele, rar.sor.spa, 
                                             rar.sim.ele, rar.sim.spa, 
                                             rar.nes.ele, rar.nes.spa))
    
    rare.partial.mantels$year <- i
    rare.partial.mantels$seasn <- j
    rare.partial.mantels$test <- rownames(rare.partial.mantels)
    rare.partial.mantels$div <- 
      substring(rare.partial.mantels$test, 1, 3)
    rare.partial.mantels$index <- 
      substring(rare.partial.mantels$test, 5, 7)
    rare.partial.mantels$gradient <- 
      substring(rare.partial.mantels$test, 9, 11)
    
    rare.partial.mantels.compR <- rbind(rare.partial.mantels.compR, 
                                        rare.partial.mantels)
    rownames(rare.partial.mantels.compR) <- 
      seq(1, nrow(rare.partial.mantels.compR), 1)
    
    
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
    
    taxonomic.partial.mantels.compR <- rbind(taxonomic.partial.mantels.compR, 
                                             taxonomic.partial.mantels)
    rownames(taxonomic.partial.mantels.compR) <- 
      seq(1, nrow(taxonomic.partial.mantels.compR), 1)
    
    #     #functional partial mantels
    #     fun.sor.ele <- mantel(fun.dist.bsor ~ ele.dist + spa.dist, 
    #                           nperm = 1000, nboot = 1000)
    #     fun.sor.spa <- mantel(fun.dist.bsor ~ spa.dist + ele.dist, 
    #                           nperm = 1000, nboot = 1000)
    #     
    #     fun.sim.ele <- mantel(fun.dist.bsim ~ ele.dist + spa.dist, 
    #                           nperm = 1000, nboot = 1000)
    #     fun.sim.spa <- mantel(fun.dist.bsim ~ spa.dist + ele.dist, 
    #                           nperm = 1000, nboot = 1000)
    #     
    #     fun.sne.ele <- mantel(fun.dist.bsne ~ ele.dist + spa.dist, 
    #                           nperm = 1000, nboot = 1000)
    #     fun.sne.spa <- mantel(fun.dist.bsne ~ spa.dist + ele.dist, 
    #                           nperm = 1000, nboot = 1000)
    #     
    #     functional.partial.mantels <- data.frame(rbind(fun.sor.ele, fun.sor.spa, 
    #                                                    fun.sim.ele, fun.sim.spa, 
    #                                                    fun.sne.ele, fun.sne.spa))
    #     functional.partial.mantels$year <- i
    #     functional.partial.mantels$seasn <- j
    #     functional.partial.mantels$test <- rownames(functional.partial.mantels)
    #     functional.partial.mantels$div <- 
    #       substring(functional.partial.mantels$test, 1, 3)
    #     functional.partial.mantels$index <- 
    #       substring(functional.partial.mantels$test, 5, 7)
    #     functional.partial.mantels$gradient <- 
    #       substring(functional.partial.mantels$test, 9, 11)
    #     
    #     functional.partial.mantels.compS <- rbind(functional.partial.mantels.compS, 
    #                                               functional.partial.mantels)
    #     rownames(functional.partial.mantels.compS) <- 
    #       seq(1, nrow(functional.partial.mantels.compS), 1)
  }
  
  
}

#comm.dissimR
#rare.dissimR

taxonomic.partial.mantels.compR
taxonomic.partial.mantels.compR[taxonomic.partial.mantels.compR$gradient == "ele" &
                                 taxonomic.partial.mantels.compR$index == "sor", ]

#rare.partial.mantels.compR

endtime <- Sys.time()
endtime - starttime  # 45 mins





par(mfrow = c(3, 4))
metrics <- c("bsor", "bsim", "bnes")
for (k in unique(metrics)){
  for(m in unique(c("tax", "rar"))){
    if (m == "tax"){
      dissim <- comm.dissimR
      ele.list <- ele.dissimR
      spa.list <- spa.dissimR
      partial.mantels.comp <- taxonomic.partial.mantels.compR
    }
    else {
      dissim <- comm.dissimR
      ele.list <- ele.dissim.rR
      spa.list <- rare.dissimR
      partial.mantels.comp <- rare.partial.mantels.compR
    }
    
    
    if(k == "bnes" && m == "fun"){
      k <- "bsne"
    }
    else {
      k <- k
    }
    
    
    for (j in unique(commS$seasn)){
      ylabel <- paste("Community Dissimilarity Residuals:", k, sep = " ")
      plot(x = c(-500, 500), y = c(-0.25, 0.25), ylab = ylabel, xlab = "Elevational Distance Residuals", type = "n", 
           main = paste(k, m, j))
      
      for (i in unique(commS$year)){
        c.dist <- dissim[[paste(i, j, k, sep = "")]]
        e.dist <- ele.list[[paste(i, j, sep = "")]]
        
        
        if (m == "tax"){
          s.dist <- spa.list[[paste(i, j, sep = "")]]
        }
        else{
          s.dist <- spa.list[[paste(i, j, k, sep = "")]]
          
        }
        
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







par(mfrow = c(2, 7))
for (i in unique(comm$year)){
  for (j in unique(comm$seasn)){
    dat <- comm[comm$year == i & comm$seasn == j, ]
    #dat <- comm
    dat <- subset(dat, select = - c(s.id, Site, year, rep, seasn))
    present <- as.vector(names(dat[, - (which(colSums(dat) == 0))]))
    y2 <- y[which(rownames(y) %in% present), ]
    ran.temp <- rangesize(y2)
    hist(ran.temp, main = paste(i, j), xlim = c(0, 2200), breaks = 8)
  }
  
}

beta.pair(t(y2))
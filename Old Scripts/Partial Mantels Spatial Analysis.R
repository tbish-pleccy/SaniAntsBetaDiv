# Spatial Analysis with Partial Mantels

# Get packages
library(FD)
library(betapart)
library(vegan)
library(boot)
library(ecodist)
library(sp)

# Create empty objects

taxonomic.partial.mantels.comp <- NULL
functional.partial.mantels.comp <- NULL
comm.dissim <- list()
fun.dissim <- list()
ele.dissim <- list()
spa.dissim <- list()



# Start loop
starttime <- Sys.time()

for (i in unique(comm$year)){
  for (j in unique(comm$seasn)){
    c.ab <- comm[comm$year == i & comm$seasn == j, ]  # c.ab is abundance data
    
    c.pa <- subset(c.ab, select = - c(s.id, Site, rep, year, seasn))
    rownames(c.pa) <- c.ab$s.id
    c.pa[c.pa > 0] <- 1  # from here, c.pa is presence/absence data
    c.pa <- c.pa[- which(rowSums(c.pa) < 5), ]  # use only comms with >5 species
    absent.species <- which(colSums(c.pa) == 0)
    c.pa <- c.pa[, - absent.species]  # remove empty species
    c.pa <- c.pa[as.vector(which(rowSums(c.pa) > 0)), ]  # remove empty communities
    c.pa <- c.pa[, order(colnames(c.pa))]  # alphabetically order!
    present.species <- as.vector(names(c.pa))  # remeber present species
    
    
    # Create appropriate site info frame.
    site.info.temp <- site.info[site.info$year == i & site.info$seasn == j, ]
    rownames(site.info.temp) <- site.info.temp$s.id
    site.info.temp$s.id <- as.numeric(site.info.temp$s.id)
    site.info.temp <- site.info.temp[order(site.info.temp$s.id), ]
    
    rownames(c.pa) <- as.numeric(rownames(c.pa))
    site.info.temp <- site.info.temp[which(site.info.temp$s.id %in% 
                                             as.numeric(rownames(c.pa))), ]
    
    # Create altitudinal distance matrix
    altitudes <- data.frame(site.info.temp[, "altitude"])
    rownames(altitudes) <- rownames(site.info.temp)
    ele.dist <- dist(altitudes)
    
    ele.dissim[[paste(i, j, sep = "")]] <- ele.dist
    
    # create spatial distance matrix
    coord <- as.matrix(subset(site.info.temp, select = c(x, y)))
    spa.dist <- as.dist(spDists(coord, longlat = TRUE))
    
    spa.dissim[[paste(i, j, sep = "")]] <- spa.dist
    
    # Create community distance matrices
    comm.dist <- beta.pair(c.pa)
    comm.dist.bsor <- comm.dist$beta.sor
    comm.dist.bsim <- comm.dist$beta.sim
    comm.dist.bnes <- comm.dist$beta.sne
    
    comm.dissim[[paste(i, j, "bsor", sep = "")]] <- comm.dist.bsor
    comm.dissim[[paste(i, j, "bsim", sep = "")]] <- comm.dist.bsim
    comm.dissim[[paste(i, j, "bnes", sep = "")]] <- comm.dist.bnes    
    
    # Functional beta diversity
    traits.temp <- pcoa.axes[rownames(pcoa.axes) %in% present.species, ]
    fun.dist <- functional.beta.pair(x = c.pa, traits = traits.temp)
    fun.dist.bsor <- fun.dist$funct.beta.sor
    fun.dist.bsim <- fun.dist$funct.beta.sim
    fun.dist.bsne <- fun.dist$funct.beta.sne
    
    fun.dissim[[paste(i, j, "bsor", sep = "")]] <- fun.dist.bsor
    fun.dissim[[paste(i, j, "bsim", sep = "")]] <- fun.dist.bsim
    fun.dissim[[paste(i, j, "bsne", sep = "")]] <- fun.dist.bsne
    
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
    
    taxonomic.partial.mantels.comp <- rbind(taxonomic.partial.mantels.comp, 
                                            taxonomic.partial.mantels)
    rownames(taxonomic.partial.mantels.comp) <- 
      seq(1, nrow(taxonomic.partial.mantels.comp), 1)
    
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
    
    functional.partial.mantels.comp <- rbind(functional.partial.mantels.comp, 
                                            functional.partial.mantels)
    rownames(functional.partial.mantels.comp) <- 
      seq(1, nrow(functional.partial.mantels.comp), 1)
  }
  
  
}
  
  
taxonomic.partial.mantels.comp
functional.partial.mantels.comp

endtime <- Sys.time()
endtime - starttime  # 2.65















partial_plot <- function(i, j, k){
  c.dist <- comm.dissim[[paste(i, j, k, sep = "")]]
  e.dist <- ele.dissim[[paste(i, j, sep = "")]]
  s.dist <- spa.dissim[[paste(i, j, sep = "")]]
  
  residual_vectors <- cbind(c.dist, e.dist, s.dist)
  resid_vecs <- as.data.frame(residual_vectors)
  
  comm_residualsG <- resid(lm(c.dist ~ e.dist, data = resid_vecs))
  eco_residuals <- resid(lm(s.dist ~ e.dist, data = resid_vecs))
  comm_residualsE <- resid(lm(c.dist ~ s.dist, data = resid_vecs))
  geo_residuals <- resid(lm(e.dist ~ s.dist, data = resid_vecs))
  
  ind <- substring(k, 2, 4)
  
  elep <- taxonomic.partial.mantels.comp[taxonomic.partial.mantels.comp$index == ind & 
                                           taxonomic.partial.mantels.comp$gradient == "ele" & 
                                           taxonomic.partial.mantels.comp$year == i & 
                                           taxonomic.partial.mantels.comp$seasn == j, "pval3"]
  spap <- taxonomic.partial.mantels.comp[taxonomic.partial.mantels.comp$index == ind & 
                                           taxonomic.partial.mantels.comp$gradient == "spa" & 
                                           taxonomic.partial.mantels.comp$year == i & 
                                           taxonomic.partial.mantels.comp$seasn == j, "pval3"]
  par(mfrow = c(1, 2))
  ylabel <- paste("Community Dissimilarity Residuals:", k, sep = " ")
  plot(geo_residuals, comm_residualsE, ylab = ylabel, xlab = "Elevational Distance Residuals", 
       pch = 21, bg = "black", col = "white")
  eline <- ifelse(elep < 0.05, 1, 0)
  abline(lm(comm_residualsE~-1+geo_residuals), col = "red", lty = eline)
  plot(eco_residuals, comm_residualsG, ylab = ylabel, xlab = "Geographical Distance Residuals", 
       pch = 21, bg = "black", col = "white")
  sline <- ifelse(spap < 0.05, 1, 0)
  abline(lm(comm_residualsG~-1+eco_residuals), col = "red", lty = sline)
}

partial_plot(2008, "J", "bsim")
partial_plot(2012, "S", "bsim")
partial_plot(2012, "S", "bnes")





for (i in unique(comm$year)){
  for (j in unique(comm$seasn)){
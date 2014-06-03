
starttime <- Sys.time()
library(FD)
library(betapart)
library(vegan)
library(boot)
library(ecodist)
library(sp)

comm.dissim.time <- list()
fun.dissim.time <- list()
time.dissim <- list()

for (i in unique(comm$Site)){
  #for (j in unique(comm$seasn)){
    c.ab <- comm[comm$Site == i, ]  # c.ab is abundance data
    
    c.pa <- subset(c.ab, select = - c(s.id, Site, rep, year, seasn))
    rownames(c.pa) <- c.ab$s.id
    c.pa[c.pa > 0] <- 1  # from here, c.pa is presence/absence data
    c.pa <- c.pa[- which(rowSums(c.pa) < 5), ]  # use only comms with >5 species
    absent.species <- which(colSums(c.pa) == 0)
    c.pa <- c.pa[, - absent.species]  # remove empty species
    c.pa <- c.pa[as.vector(which(rowSums(c.pa) > 0)), ]  # remove empty communities
    c.pa <- c.pa[, order(colnames(c.pa))]  # alphabetically order!
    present.species <- as.vector(names(c.pa))  # remember present species
    
    
    # Create appropriate site info frame.
    site.info.temp <- site.info[site.info$Site == i, ]
    rownames(site.info.temp) <- site.info.temp$s.id
    site.info.temp$s.id <- as.numeric(site.info.temp$s.id)
    site.info.temp <- site.info.temp[order(site.info.temp$s.id), ]
    
    rownames(c.pa) <- as.numeric(rownames(c.pa))
    site.info.temp <- site.info.temp[which(site.info.temp$s.id %in% 
                                             as.numeric(rownames(c.pa))), ]
        
    # Create time distance matrix
    time <- data.frame(site.info.temp[, "period"])
    rownames(time) <- rownames(site.info.temp)
    time.dist <- dist(time)
    time.dissim[[paste(i, sep = "")]] <- time.dist
    
    
    
    
    # Create community distance matrices
    comm.dist <- beta.pair(c.pa)
    comm.dist.bsor <- comm.dist$beta.sor
    comm.dist.bsim <- comm.dist$beta.sim
    comm.dist.bnes <- comm.dist$beta.sne
    
    comm.dissim.time[[paste(i, "bsor", sep = "")]] <- comm.dist.bsor
    comm.dissim.time[[paste(i, "bsim", sep = "")]] <- comm.dist.bsim
    comm.dissim.time[[paste(i, "bnes", sep = "")]] <- comm.dist.bnes    
    
    # Functional beta diversity
    traits.temp <- pcoa.axes[rownames(pcoa.axes) %in% present.species, ]
    fun.dist <- functional.beta.pair(x = c.pa, traits = traits.temp)
    fun.dist.bsor <- fun.dist$funct.beta.sor
    fun.dist.bsim <- fun.dist$funct.beta.sim
    fun.dist.bsne <- fun.dist$funct.beta.sne
    
    fun.dissim.time[[paste(i, "bsor", sep = "")]] <- fun.dist.bsor
    fun.dissim.time[[paste(i, "bsim", sep = "")]] <- fun.dist.bsim
    fun.dissim.time[[paste(i, "bsne", sep = "")]] <- fun.dist.bsne
  #}
    print(paste("Progress: ",i, "/8", sep = ""))

}


endtime <- Sys.time()
endtime - starttime

Site <- 5
metric <- "bsim"
composition <- "taxonomic"
season <- "J"

betatime.plot <- function(Site, metric, composition, season){
  
  if (composition == "taxonomic"){
    comp.dist.list <- comm.dissim.time
  }
  else {
    comp.dist.list <- fun.dissim.time
  }
  
  comm.dist <- comp.dist.list[[paste(Site, metric, sep = "")]]
  grad.dist <- time.dissim[[paste(Site, sep = "")]]
  
  select <- site.info[site.info$Site == Site & site.info$seasn == season, "s.id"]
  
  tf <- labels(comm.dist) %in% select
  comm.dist <- as.dist(full(comm.dist)[tf, tf])
  grad.dist <- as.dist(full(grad.dist)[tf, tf])
  
  
  
  ylabel <- metric
  xlabel <- "Time Difference"
  mlabel <- paste(Site)
  
  
  plot(comm.dist ~ grad.dist, pch = 21, col = "white", bg = "black", 
       ylab = ylabel, xlab = xlabel, main = mlabel, ylim = c(0, 1))
  
  if (mantel(comm.dist ~ grad.dist)[4] < 0.05) {
    abline((lm(comm.dist ~ grad.dist)))
    
  }
  
  
  
}

par(mfrow = c(2, 4))
betatime.plot(1, "bsor", "functional", "J")
betatime.plot(2, "bsor", "functional", "J")
betatime.plot(3, "bsor", "functional", "J")
betatime.plot(4, "bsor", "functional", "J")
betatime.plot(5, "bsor", "functional", "J")
betatime.plot(6, "bsor", "functional", "J")
betatime.plot(7, "bsor", "functional", "J")
betatime.plot(8, "bsor", "functional", "J")


betatime.plot(1, "bsim", "functional", "J")
betatime.plot(2, "bsim", "functional", "J")
betatime.plot(3, "bsim", "functional", "J")
betatime.plot(4, "bsim", "functional", "J")
betatime.plot(5, "bsim", "functional", "J")
betatime.plot(6, "bsim", "functional", "J")
betatime.plot(7, "bsim", "functional", "J")
betatime.plot(8, "bsim", "functional", "J")


betatime.plot(1, "bsne", "functional", "J")
betatime.plot(2, "bsne", "functional", "J")
betatime.plot(3, "bsne", "functional", "J")
betatime.plot(4, "bsne", "functional", "J")
betatime.plot(5, "bsne", "functional", "J")
betatime.plot(6, "bsne", "functional", "J")
betatime.plot(7, "bsne", "functional", "J")
betatime.plot(8, "bsne", "functional", "J")


betatime.plot(1, "bsor", "taxonomic", "J")
betatime.plot(2, "bsor", "taxonomic", "J")
betatime.plot(3, "bsor", "taxonomic", "J")
betatime.plot(4, "bsor", "taxonomic", "J")
betatime.plot(5, "bsor", "taxonomic", "J")
betatime.plot(6, "bsor", "taxonomic", "J")
betatime.plot(7, "bsor", "taxonomic", "J")
betatime.plot(8, "bsor", "taxonomic", "J")



betatime.plot(1, "bsim", "taxonomic", "J")
betatime.plot(2, "bsim", "taxonomic", "J")
betatime.plot(3, "bsim", "taxonomic", "J")
betatime.plot(4, "bsim", "taxonomic", "J")
betatime.plot(5, "bsim", "taxonomic", "J")
betatime.plot(6, "bsim", "taxonomic", "J")
betatime.plot(7, "bsim", "taxonomic", "J")
betatime.plot(8, "bsim", "taxonomic", "J")


betatime.plot(1, "bnes", "taxonomic", "J")
betatime.plot(2, "bnes", "taxonomic", "J")
betatime.plot(3, "bnes", "taxonomic", "J")
betatime.plot(4, "bnes", "taxonomic", "J")
betatime.plot(5, "bnes", "taxonomic", "J")
betatime.plot(6, "bnes", "taxonomic", "J")
betatime.plot(7, "bnes", "taxonomic", "J")
betatime.plot(8, "bnes", "taxonomic", "J")





betatime.plot(1, "bsor", "functional", "S")
betatime.plot(2, "bsor", "functional", "S")
betatime.plot(3, "bsor", "functional", "S")
betatime.plot(4, "bsor", "functional", "S")
betatime.plot(5, "bsor", "functional", "S")
betatime.plot(6, "bsor", "functional", "S")
betatime.plot(7, "bsor", "functional", "S")
betatime.plot(8, "bsor", "functional", "S")


betatime.plot(1, "bsim", "functional", "S")
betatime.plot(2, "bsim", "functional", "S")
betatime.plot(3, "bsim", "functional", "S")
betatime.plot(4, "bsim", "functional", "S")
betatime.plot(5, "bsim", "functional", "S")
betatime.plot(6, "bsim", "functional", "S")
betatime.plot(7, "bsim", "functional", "S")
betatime.plot(8, "bsim", "functional", "S")


betatime.plot(1, "bsne", "functional", "S")
betatime.plot(2, "bsne", "functional", "S")
betatime.plot(3, "bsne", "functional", "S")
betatime.plot(4, "bsne", "functional", "S")
betatime.plot(5, "bsne", "functional", "S")
betatime.plot(6, "bsne", "functional", "S")
betatime.plot(7, "bsne", "functional", "S")
betatime.plot(8, "bsne", "functional", "S")


betatime.plot(1, "bsor", "taxonomic", "S")
betatime.plot(2, "bsor", "taxonomic", "S")
betatime.plot(3, "bsor", "taxonomic", "S")
betatime.plot(4, "bsor", "taxonomic", "S")
betatime.plot(5, "bsor", "taxonomic", "S")
betatime.plot(6, "bsor", "taxonomic", "S")
betatime.plot(7, "bsor", "taxonomic", "S")
betatime.plot(8, "bsor", "taxonomic", "S")



betatime.plot(1, "bsim", "taxonomic", "S")
betatime.plot(2, "bsim", "taxonomic", "S")
betatime.plot(3, "bsim", "taxonomic", "S")
betatime.plot(4, "bsim", "taxonomic", "S")
betatime.plot(5, "bsim", "taxonomic", "S")
betatime.plot(6, "bsim", "taxonomic", "S")
betatime.plot(7, "bsim", "taxonomic", "S")
betatime.plot(8, "bsim", "taxonomic", "S")


betatime.plot(1, "bnes", "taxonomic", "S")
betatime.plot(2, "bnes", "taxonomic", "S")
betatime.plot(3, "bnes", "taxonomic", "S")
betatime.plot(4, "bnes", "taxonomic", "S")
betatime.plot(5, "bnes", "taxonomic", "S")
betatime.plot(6, "bnes", "taxonomic", "S")
betatime.plot(7, "bnes", "taxonomic", "S")
betatime.plot(8, "bnes", "taxonomic", "S")








betatime.plot2 <- function(Site, metric, composition, season){
  
  if (composition == "taxonomic"){
    comp.dist.list <- comm.dissim.time
  }
  else {
    comp.dist.list <- fun.dissim.time
  }
  
  comm.dist <- comp.dist.list[[paste(Site, metric, sep = "")]]
  grad.dist <- time.dissim[[paste(Site, sep = "")]]
  
  #select <- site.info[site.info$Site == Site & site.info$seasn == season, "s.id"]
  
  #tf <- labels(comm.dist) %in% select
  #comm.dist <- as.dist(full(comm.dist)[tf, tf])
  #grad.dist <- as.dist(full(grad.dist)[tf, tf])
  
  
  
  ylabel <- metric
  xlabel <- "Time Difference"
  mlabel <- paste(Site)
  
  
  plot(comm.dist ~ grad.dist, pch = 21, col = "white", bg = "black", 
       ylab = ylabel, xlab = xlabel, main = mlabel, ylim = c(0, 1))
  
  if (mantel(comm.dist ~ grad.dist)[4] < 0.05) {
    abline((lm(comm.dist ~ grad.dist)))
    
  }
  
  
  
}




par(mfrow = c(2, 4))
betatime.plot2(1, "bsor", "functional", "J")
betatime.plot2(2, "bsor", "functional", "J")
betatime.plot2(3, "bsor", "functional", "J")
betatime.plot2(4, "bsor", "functional", "J")
betatime.plot2(5, "bsor", "functional", "J")
betatime.plot2(6, "bsor", "functional", "J")
betatime.plot2(7, "bsor", "functional", "J")
betatime.plot2(8, "bsor", "functional", "J")


betatime.plot2(1, "bsim", "functional", "J")
betatime.plot2(2, "bsim", "functional", "J")
betatime.plot2(3, "bsim", "functional", "J")
betatime.plot2(4, "bsim", "functional", "J")
betatime.plot2(5, "bsim", "functional", "J")
betatime.plot2(6, "bsim", "functional", "J")
betatime.plot2(7, "bsim", "functional", "J")
betatime.plot2(8, "bsim", "functional", "J")


betatime.plot2(1, "bsne", "functional", "J")
betatime.plot2(2, "bsne", "functional", "J")
betatime.plot2(3, "bsne", "functional", "J")
betatime.plot2(4, "bsne", "functional", "J")
betatime.plot2(5, "bsne", "functional", "J")
betatime.plot2(6, "bsne", "functional", "J")
betatime.plot2(7, "bsne", "functional", "J")
betatime.plot2(8, "bsne", "functional", "J")


betatime.plot2(1, "bsor", "taxonomic", "J")
betatime.plot2(2, "bsor", "taxonomic", "J")
betatime.plot2(3, "bsor", "taxonomic", "J")
betatime.plot2(4, "bsor", "taxonomic", "J")
betatime.plot2(5, "bsor", "taxonomic", "J")
betatime.plot2(6, "bsor", "taxonomic", "J")
betatime.plot2(7, "bsor", "taxonomic", "J")
betatime.plot2(8, "bsor", "taxonomic", "J")



betatime.plot2(1, "bsim", "taxonomic", "J")
betatime.plot2(2, "bsim", "taxonomic", "J")
betatime.plot2(3, "bsim", "taxonomic", "J")
betatime.plot2(4, "bsim", "taxonomic", "J")
betatime.plot2(5, "bsim", "taxonomic", "J")
betatime.plot2(6, "bsim", "taxonomic", "J")
betatime.plot2(7, "bsim", "taxonomic", "J")
betatime.plot2(8, "bsim", "taxonomic", "J")


betatime.plot2(1, "bnes", "taxonomic", "J")
betatime.plot2(2, "bnes", "taxonomic", "J")
betatime.plot2(3, "bnes", "taxonomic", "J")
betatime.plot2(4, "bnes", "taxonomic", "J")
betatime.plot2(5, "bnes", "taxonomic", "J")
betatime.plot2(6, "bnes", "taxonomic", "J")
betatime.plot2(7, "bnes", "taxonomic", "J")
betatime.plot2(8, "bnes", "taxonomic", "J")





betatime.plot2(1, "bsor", "functional", "S")
betatime.plot2(2, "bsor", "functional", "S")
betatime.plot2(3, "bsor", "functional", "S")
betatime.plot2(4, "bsor", "functional", "S")
betatime.plot2(5, "bsor", "functional", "S")
betatime.plot2(6, "bsor", "functional", "S")
betatime.plot2(7, "bsor", "functional", "S")
betatime.plot2(8, "bsor", "functional", "S")


betatime.plot2(1, "bsim", "functional", "S")
betatime.plot2(2, "bsim", "functional", "S")
betatime.plot2(3, "bsim", "functional", "S")
betatime.plot2(4, "bsim", "functional", "S")
betatime.plot2(5, "bsim", "functional", "S")
betatime.plot2(6, "bsim", "functional", "S")
betatime.plot2(7, "bsim", "functional", "S")
betatime.plot2(8, "bsim", "functional", "S")


betatime.plot2(1, "bsne", "functional", "S")
betatime.plot2(2, "bsne", "functional", "S")
betatime.plot2(3, "bsne", "functional", "S")
betatime.plot2(4, "bsne", "functional", "S")
betatime.plot2(5, "bsne", "functional", "S")
betatime.plot2(6, "bsne", "functional", "S")
betatime.plot2(7, "bsne", "functional", "S")
betatime.plot2(8, "bsne", "functional", "S")


betatime.plot2(1, "bsor", "taxonomic", "S")
betatime.plot2(2, "bsor", "taxonomic", "S")
betatime.plot2(3, "bsor", "taxonomic", "S")
betatime.plot2(4, "bsor", "taxonomic", "S")
betatime.plot2(5, "bsor", "taxonomic", "S")
betatime.plot2(6, "bsor", "taxonomic", "S")
betatime.plot2(7, "bsor", "taxonomic", "S")
betatime.plot2(8, "bsor", "taxonomic", "S")



betatime.plot2(1, "bsim", "taxonomic", "S")
betatime.plot2(2, "bsim", "taxonomic", "S")
betatime.plot2(3, "bsim", "taxonomic", "S")
betatime.plot2(4, "bsim", "taxonomic", "S")
betatime.plot2(5, "bsim", "taxonomic", "S")
betatime.plot2(6, "bsim", "taxonomic", "S")
betatime.plot2(7, "bsim", "taxonomic", "S")
betatime.plot2(8, "bsim", "taxonomic", "S")


betatime.plot2(1, "bnes", "taxonomic", "S")
betatime.plot2(2, "bnes", "taxonomic", "S")
betatime.plot2(3, "bnes", "taxonomic", "S")
betatime.plot2(4, "bnes", "taxonomic", "S")
betatime.plot2(5, "bnes", "taxonomic", "S")
betatime.plot2(6, "bnes", "taxonomic", "S")
betatime.plot2(7, "bnes", "taxonomic", "S")
betatime.plot2(8, "bnes", "taxonomic", "S")












i <- "a"
x <- "J"
y <- 5
z <- "bsim"

comp.frame.final <- NULL
for (i in c("a", "b", "c", "d")){
  for (x in levels(site.info$seasn)){
    for (y in unique(site.info$Site)){
      for (z in c("bsor", "bsim", "bnes")){
    
        data <- site.info[site.info$Site == y & site.info$seasn == x, ]
        
        rep.frame <- data[data$rep == paste(y, i, sep = ""), ]
        comp1 <- rep.frame[1:nrow(rep.frame)-1, "s.id"]
        comp2 <- rep.frame[2:nrow(rep.frame), "s.id"]
        comp.frame <- data.frame(comp1, comp2)
        
        c.name <- NULL
        for (j in rownames(comp.frame)){
          c.name.temp <- paste(rep.frame[rep.frame$s.id == comp.frame[j, "comp1"], "year"],
                               rep.frame[rep.frame$s.id == comp.frame[j, "comp2"], "year"], 
                               sep = " - ")
          c.name <- c(c.name, c.name.temp)    
        }
        comp.frame <- cbind(comp.frame, c.name)
        
#         if (composition == "taxonomic"){
#           comp.dist.list <- comm.dissim.time
#         }
#         else {
#           comp.dist.list <- fun.dissim.time
#         }
        
        comm.dist <- comp.dist.list[[paste(y, z, sep = "")]]
        
        distances <- NULL
        for (k in 1:nrow(comp.frame)){
          select <- c(as.character(comp.frame[k, 1]),as.character(comp.frame[k, 2]))
          tf <- labels(comm.dist) %in% select
          distances.temp <- as.dist(full(comm.dist)[tf, tf])
          if (length(distances.temp) == 0){
            distances.temp <- NA
          }
          else {
            distances.temp <- distances.temp
          }
          distances <- c(distances, distances.temp)
          
        }
        comp.frame <- cbind(comp.frame, distances)
        
        comp.frame <- cbind(comp.frame, z)
        comp.frame <- cbind(comp.frame, x)
        comp.frame <- cbind(comp.frame, y)
        comp.frame.final <- rbind(comp.frame.final, comp.frame)   
        comp.frame <- NULL
        
        
      }
      
      
      
      
    }
  }
  
    
  }
  


x <- "J"
y <- 1
z <- "bsor"

time.plot <- function(x, y, z){
  data <- comp.frame.final[comp.frame.final$x == x & comp.frame.final$y == y & 
                           comp.frame.final$z == z, ]
  plot(distances ~ c.name, data = data, ylab = z, main = unique(paste(x, y)))
  
}


time.plot("J", 1, "bsor")
time.plot("J", 2, "bsor")
time.plot("J", 3, "bsor")
time.plot("J", 4, "bsor")
time.plot("J", 5, "bsor")
time.plot("J", 6, "bsor")
time.plot("J", 7, "bsor")
time.plot("J", 8, "bsor")


time.plot("J", 1, "bsim")
time.plot("J", 2, "bsim")
time.plot("J", 3, "bsim")
time.plot("J", 4, "bsim")
time.plot("J", 5, "bsim")
time.plot("J", 6, "bsim")
time.plot("J", 7, "bsim")
time.plot("J", 8, "bsim")


time.plot("J", 1, "bnes")
time.plot("J", 2, "bnes")
time.plot("J", 3, "bnes")
time.plot("J", 4, "bnes")
time.plot("J", 5, "bnes")
time.plot("J", 6, "bnes")
time.plot("J", 7, "bnes")
time.plot("J", 8, "bnes")